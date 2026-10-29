"""Regression cases found by the 2026-10-19 scientific-code audit."""
from pathlib import Path
from types import SimpleNamespace
import json

import numpy as np
import pytest

from fv_energy_analysis.geometry import hcdr3_segments, paired_atoms, structural_rmsds, kabsch
from fv_energy_analysis.models import Atom, Residue, ResidueKey, Structure
from fv_energy_analysis.rotamers import CHI_ATOMS, angle_bin, compare_rotamers
from fv_energy_analysis.cli import _validated_pairs, _mutation_rows
from fv_energy_analysis.batch_workflow import _hbond_saltbridge_summary, _candidate_network_summary, _parent_cache_key
from fv_energy_analysis.validation import validate_ligand_against_params
from fv_energy_analysis.fingerprint_reports import _grouped_columns
from fv_energy_analysis.fast_workflow import _link_or_copy
from fv_energy_analysis.rosetta import read_residue_energy_breakdown
from fv_energy_analysis.structure import read_pdb
from fv_energy_analysis.estrogen_architecture import face_to_face_pi_geometry
from fv_energy_analysis.fast_workflow import validate_scoring_reuse
from fv_energy_analysis.hbonds import aggregate_hbond_csv


def structures():
    sequence = 'ARGYYGSSNAD'
    codes = {'A':'ALA','R':'ARG','G':'GLY','Y':'TYR','S':'SER','N':'ASN','D':'ASP'}
    p, m = [], []
    for i, aa in enumerate(sequence):
        key = ResidueKey('H', 97+i)
        xyz = np.array([i, np.sin(i), np.cos(i)])
        def residue(shift):
            r = Residue(key, codes[aa])
            r.atoms['CA'] = Atom(i+1,'CA','',codes[aa],key,xyz+shift,'C',1.0,i)
            return r
        p.append(residue(0))
        m.append(residue(np.array([0., 0., 0.5 if 3 <= i < 8 else 0.])))
    return Structure(Path('p'),[],p), Structure(Path('m'),[],m), {'H':(p,m)}, {r.key:'CDR3' for r in p}


def test_asn_and_gln_chi_definitions_match_rosetta():
    assert CHI_ATOMS['ASN'][1] == ('CA','CB','CG','OD1')
    assert CHI_ATOMS['GLN'][2] == ('CB','CG','CD','OE1')


@pytest.mark.parametrize("chain_type,region,whole,anchors,expected", [
    ("H", "CDR3", 5.0, 1.0, "performed"),
    ("H", "CDR3", 1.0, 4.0, "skipped_cdr_rmsd_above_cutoff"),
    ("H", "CDR3", 5.0, 3.0, "performed"),
    ("H", "CDR3", 1.0, 3.001, "skipped_cdr_rmsd_above_cutoff"),
    ("H", "CDR3", 1.0, None, "skipped_cdr_rmsd_unavailable"),
    ("H", "CDR3", 1.0, float("nan"), "skipped_cdr_rmsd_unavailable"),
    ("H", "CDR3", 1.0, float("inf"), "skipped_cdr_rmsd_unavailable"),
    ("H", "CDR1", 4.0, 1.0, "skipped_cdr_rmsd_above_cutoff"),
    ("H", "CDR2", 2.0, 4.0, "performed"),
    ("K", "CDR3", 4.0, 1.0, "skipped_cdr_rmsd_above_cutoff"),
    ("L", "CDR3", 2.0, 4.0, "performed"),
    ("H", "FR3", 5.0, 4.0, "performed"),
])
def test_hcdr3_rotamer_gate_uses_anchors_only(monkeypatch, chain_type, region, whole, anchors, expected):
    from unittest.mock import Mock
    from fv_energy_analysis import cli
    # Heavy-chain identity comes from annotation, not the literal chain ID H.
    key = ResidueKey("A", 104)
    parent, mutant = Residue(key, "ALA"), Residue(key, "VAL")
    ann = SimpleNamespace(chain_type=chain_type, region_by_key={key: region}, number_by_key={key: "104"})
    values = {f"A_{region}": whole}
    if anchors is not None:
        values["A_CDR3_ANCHORS"] = anchors
    compare = Mock(return_value={"rotamer_method": "canonical_chi_coarse_torsion_bins_not_dunbrack"})
    monkeypatch.setattr(cli, "compare_rotamers", compare)
    row, = _mutation_rows({"A": ([parent], [mutant])}, {"A": ann}, {"A": ann}, values, 3.0)
    assert row["rotamer_analysis_status"] == expected
    assert compare.call_count == int(expected == "performed")
    expected_region = "A_CDR3_ANCHORS" if chain_type == "H" and region == "CDR3" else f"A_{region}" if region.startswith("CDR") else ""
    assert row["rotamer_gate_region"] == expected_region
    assert row["cdr_global_fit_rmsd_A"] == (whole if region.startswith("CDR") else "")


@pytest.mark.parametrize('angle,expected', [(420,'p'),(-420,'m'),(180,'t'),(-180,'t'),(120,'t'),(-120,'m'),(0,'p')])
def test_normalized_torsion_boundaries(angle, expected):
    assert angle_bin(angle) == expected


def test_nonfinite_bin_is_not_silently_trans():
    with pytest.raises(ValueError):
        angle_bin(float('nan'))


def test_tip_and_anchors_are_a_disjoint_parent_position_partition():
    p, m, pairs, regions = structures()
    # Mutation inside the motif does not change its positional selection.
    m.residues[3].name = 'ALA'
    split = hcdr3_segments(pairs, regions)
    tip, anchors = split['H_CDR3_TIP'], split['H_CDR3_ANCHORS']
    assert len(tip['parent_keys']) == 5
    assert len(anchors['parent_keys']) == 6
    assert not set(tip['parent_keys']) & set(anchors['parent_keys'])
    assert set(tip['parent_keys']) | set(anchors['parent_keys']) == set(regions)
    assert tip['mutant_keys'] == tip['parent_keys']


def test_split_rmsds_use_one_transform_and_reconstruct_whole_sse():
    p, m, pairs, regions = structures()
    _, rows = structural_rmsds(p,m,pairs,regions,atom_names=('CA',))
    rows = {r['region']:r for r in rows}
    whole, anchors, tip = [rows[name] for name in ('H_CDR3','H_CDR3_ANCHORS','H_CDR3_TIP')]
    sse = lambda r:r['n_atoms']*r['rmsd_after_global_fit_A']**2
    assert sse(whole) == pytest.approx(sse(anchors)+sse(tip), abs=1e-12)
    assert tip['rmsd_after_global_fit_A'] > 0
    assert tip['local_fit_rmsd_A'] < 1e-10


def test_missing_tip_is_explicit_not_all_anchors_or_wrong_positions():
    p,m,pairs,regions = structures()
    _, rows = structural_rmsds(p,m,pairs,regions,atom_names=('CA',),hcdr3_tip_sequence='FFFFF')
    tip = next(r for r in rows if r['region']=='H_CDR3_TIP')
    assert tip['segment_status'] == 'tip_motif_not_found'
    assert tip['rmsd_after_global_fit_A'] == ''
    assert tip['n_residues'] == 0


def test_empty_coordinate_arrays_are_n_by_three():
    p,m,_,_ = structures()
    px,mx,n = paired_atoms(p,m,[],[])
    assert px.shape == mx.shape == (0,3)
    assert n == 0


def test_global_fit_rejects_nonfinite_coordinates():
    with pytest.raises(ValueError, match='finite'):
        kabsch(np.full((3,3),np.nan),np.zeros((3,3)))


def test_equal_length_but_different_pdb_positions_is_rejected():
    p,m,_,_ = structures()
    m.residues[0].key = ResidueKey('H',96)
    with pytest.raises(ValueError, match='PDB residue numbers'):
        _validated_pairs(p,m,['H'])


def test_ros_hbond_detector_swap_is_not_a_formed_and_broken_family(tmp_path):
    change = {'parent_mutated_residue':'H:99','mutant_mutated_residue':'H:99',
              'parent_partner':'L:55','mutant_partner':'L:55','pair_scope':'aa_to_aa',
              'parent_interaction_types':'hbond','mutant_interaction_types':'rosetta_hbond',
              'formed_interaction_types':'rosetta_hbond','broken_interaction_types':'hbond'}
    counts, rows = _hbond_saltbridge_summary('test',tmp_path,[change])
    assert counts['formed_inter_chain_hbond_pair_count']==0
    assert counts['broken_inter_chain_hbond_pair_count']==0
    assert rows == []


def test_energy_changed_retained_pair_is_not_omitted_from_changed_network(tmp_path):
    counts,_ = _candidate_network_summary('test',tmp_path,[{
        'parent_mutated_residue':'H:99','parent_partner':'L:55',
        'formed_interaction_types':'none','broken_interaction_types':'none',
        'network_status':'retained_changed','delta_total':'-1.25'}])
    assert counts['network_added_pair_type_events']==0
    assert counts['network_changed_unique_pairs']==1
    assert counts['network_changed_unique_pairs_energy_delta_REU']==-1.25


def test_empty_ligand_and_empty_params_are_invalid(tmp_path):
    pdb,params = tmp_path/'empty.pdb',tmp_path/'empty.params'
    pdb.write_text('END\n'); params.write_text('')
    assert not validate_ligand_against_params(pdb,params,'X')['valid']


def test_missing_hydrogens_cannot_silently_pass_with_nohydro(tmp_path):
    pdb,params = tmp_path/'lig.pdb',tmp_path/'lig.params'
    pdb.write_text('HETATM    1  O1  LIG X   1       0.000   0.000   0.000  1.00 20.00           O  \n')
    params.write_text('ATOM O1 OH X -0.5\nATOM H1 Hpol X 0.5\nBOND O1 H1\n')
    report=validate_ligand_against_params(pdb,params,'X')
    assert not report['valid']
    assert report['missing_hydrogens']==['H1']


def test_parent_cache_invalidates_parameter_changes_at_same_path(tmp_path):
    parent,params,template = [tmp_path/name for name in ('p.pdb','l.params','l.pdb')]
    for path in (parent,params,template):path.write_text('original')
    args=SimpleNamespace(fv_chains='H,L',ligand_chain='d',parent_ligand_chain='d',parent_chain_map={},
                         weights='ref2015',params=str(params),ligand_template_pdb=str(template),scheme='imgt',cdr_definition='imgt')
    original=_parent_cache_key(args,parent)
    params.write_text('different chemistry')
    assert _parent_cache_key(args,parent)!=original


def test_linking_a_file_to_itself_does_not_delete_it(tmp_path):
    source=tmp_path/'file';source.write_text('keep')
    _link_or_copy(source,source)
    assert source.read_text()=='keep'


def test_reuse_rejects_legacy_and_differently_parameterized_scores(tmp_path):
    params=tmp_path/'lig.params';params.write_text('new')
    (tmp_path/'workflow_manifest.json').write_text(json.dumps({'status':'complete'}))
    with pytest.raises(ValueError,match='legacy or incompatible'):
        validate_scoring_reuse(tmp_path,params,'ref2015')


def test_hbond_aggregation_does_not_confuse_raw_with_weighted_energy(tmp_path):
    path=tmp_path/'hb.csv'
    path.write_text('pair_id,hbond_energy,hbond_weighted_energy_REU\nH:1--L:1,-2,-1\nH:1--L:1,-1,-0.5\n')
    row=aggregate_hbond_csv(path)[0]
    assert row['sum_individual_hbond_energy']==-3
    assert row['sum_individual_hbond_weighted_energy_REU']==-1.5
    path.write_text('pair_id,hbond_energy\nH:1--L:1,-2\n')
    assert aggregate_hbond_csv(path)[0]['sum_individual_hbond_weighted_energy_REU']==''


def test_sdf_reference_is_taken_from_the_file_not_a_default_molecule(tmp_path):
    Chem=pytest.importorskip('rdkit.Chem')
    from rdkit.Chem import AllChem
    from fv_energy_analysis.ligand_chemistry import reference_smiles_from_source
    molecule=Chem.AddHs(Chem.MolFromSmiles('CCO'))
    AllChem.EmbedMolecule(molecule,randomSeed=7)
    path=tmp_path/'ethanol.sdf'
    writer=Chem.SDWriter(str(path));writer.write(molecule);writer.close()
    assert reference_smiles_from_source(path)=='CCO'


def test_rosetta_table_wrong_pose_mapping_is_rejected(tmp_path):
    source=Path('tests/fixtures/energy.sc').read_text().replace('1H ALA','9H ALA')
    table=tmp_path/'wrong.sc';table.write_text(source)
    with pytest.raises(ValueError,match='mapping mismatch'):
        read_residue_energy_breakdown(table,read_pdb('tests/fixtures/minimal_complex.pdb'))


def test_native_rosetta_insertion_code_pdb_ids_are_supported(tmp_path):
    table=tmp_path/'inserted.sc'
    table.write_text(Path('tests/fixtures/energy.sc').read_text().replace('1H ALA','1H:A ALA'))
    structure=read_pdb('tests/fixtures/minimal_complex.pdb')
    structure.residues[0].key=ResidueKey('H',1,'A')
    assert read_residue_energy_breakdown(table,structure)[0]['pair_id']=='H:1A--X:1'


def test_discovery_studio_resolves_recorded_prepared_structure(tmp_path):
    from fv_energy_analysis.discovery_studio import _prepared_pdb
    (tmp_path/'prepared').mkdir();(tmp_path/'final').mkdir()
    pdb=tmp_path/'prepared'/'model_input_rosetta_rosetta.pdb';pdb.write_text('END\n')
    (tmp_path/'final'/'summary.json').write_text(json.dumps({'mutant':str(pdb)}))
    assert _prepared_pdb(tmp_path)==pdb


def test_grouped_summary_preserves_overall_and_adds_two_split_metrics():
    keys={c['key'] for c in _grouped_columns()}
    assert {'h_cdr3_rmsd_A','h_cdr3_anchors_rmsd_A','h_cdr3_tip_rmsd_A'} <= keys


def test_opposite_face_contacts_are_not_automatically_parallel_pi_stacks():
    theta=np.arange(6)*np.pi/3
    ring=np.column_stack((np.cos(theta),np.sin(theta),np.zeros(6)))
    assert face_to_face_pi_geometry(ring,ring+[0,0,3.4])['passes']
    perpendicular=ring[:,[0,2,1]]+[0,0,3.4]
    assert not face_to_face_pi_geometry(ring,perpendicular)['passes']
    assert not face_to_face_pi_geometry(ring,ring+[3,0,3.4])['passes']
    assert not face_to_face_pi_geometry(ring,np.array([[0,0,3],[1,0,3],[2,0,3]]))['passes']


def test_estradiol_template_restores_aromaticity_and_checks_stereochemistry():
    Chem=pytest.importorskip('rdkit.Chem')
    from rdkit.Chem import AllChem
    from fv_energy_analysis.ligand_chemistry import assign_reference_chemistry, ESTRADIOL_SMILES, ESTRADIOL_INCHIKEY
    molecule=Chem.AddHs(Chem.MolFromSmiles(ESTRADIOL_SMILES))
    assert AllChem.EmbedMolecule(molecule,randomSeed=101)==0
    broken=Chem.RWMol(molecule)
    for bond in broken.GetBonds():
        bond.SetBondType(Chem.BondType.SINGLE);bond.SetIsAromatic(False)
    for atom in broken.GetAtoms():atom.SetIsAromatic(False)
    restored=assign_reference_chemistry(broken.GetMol(),ESTRADIOL_SMILES)
    assert Chem.MolToInchiKey(restored)==ESTRADIOL_INCHIKEY
    assert sum(b.GetIsAromatic() for b in restored.GetBonds())==6
    reflected=Chem.Mol(restored)
    conformer=reflected.GetConformer()
    for i in range(reflected.GetNumAtoms()):
        position=conformer.GetAtomPosition(i);conformer.SetAtomPosition(i,(-position.x,position.y,position.z))
    with pytest.raises(ValueError,match='stereochemistry'):
        assign_reference_chemistry(reflected,ESTRADIOL_SMILES)
