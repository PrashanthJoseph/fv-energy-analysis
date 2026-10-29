"""Developer-only numerical audit; not part of the installed runtime package."""
from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path
from types import SimpleNamespace
import xml.etree.ElementTree as ET

import numpy as np

from fv_energy_analysis.cli import compare, _validated_pairs
from fv_energy_analysis.geometry import BACKBONE, paired_atoms, structural_rmsds
from fv_energy_analysis.hbonds import write_pose_hbonds
from fv_energy_analysis.numbering import annotate_fv
from fv_energy_analysis.pyrosetta_score import _initialize, _scorefunction, _active_terms, _term_name, write_pair_energies, write_interface_metrics, complete_pair_energies_array
from fv_energy_analysis.report import write_json
from fv_energy_analysis.rosetta import read_scorefile_metrics
from fv_energy_analysis.rotamers import chi_angles, CHI_ATOMS
from fv_energy_analysis.structure import read_pdb
from fv_energy_analysis.validation import validate_ligand_against_params


def _csv(path):
    with Path(path).open(newline='', encoding='utf-8') as stream:
        return list(csv.DictReader(stream))


def run_sanity(result_dir: Path, params: Path, template: Path, output: Path, weights='ref2015'):
    """Check actual APIs, not only synthetic parsers, on one prepared mutant.

    Independent oracles: direct two-body score evaluation, onebody+pair total
    reconstruction, Cartesian separation in another direction, BioPython's fit,
    Rosetta's own chi angles, RDKit stereochemistry, and raw detector records.
    """
    if output.exists():
        raise ValueError('Sanity output must be a new directory; existing evidence is not overwritten')
    result_dir, output = result_dir.resolve(), output.resolve()
    if result_dir == output or result_dir in output.parents:
        raise ValueError('Choose an output outside the saved scientific result directory')
    output.mkdir(parents=True)
    old_summary = json.loads((result_dir / 'automatic_summary.json').read_text())
    metadata = json.loads((result_dir / 'final' / 'summary.json').read_text())
    parent = Path(metadata['parent'])
    mutant = Path(metadata['mutant'])
    if not parent.is_file() or not mutant.is_file():
        raise FileNotFoundError('Run this checker in the same filesystem environment as the saved paths')
    from fv_energy_analysis.estrogen_architecture import analyze_estrogen_architecture
    architecture = {label: analyze_estrogen_architecture(pdb, result_dir/f'arpeggio_{label}'/f'{label}_arpeggio.json',
                    params_path=params, fv_chains=set(metadata['fv_chains'])) for label,pdb in (('parent',parent),('mutant',mutant))}
    checks, warnings = [], []

    def check(name, actual, expected, tolerance=1e-6):
        if isinstance(actual, (float, int)) and not isinstance(actual, bool):
            okay = math.isfinite(actual) and abs(actual - expected) <= tolerance
        else:
            okay = actual == expected
        checks.append({'check': name, 'passed': bool(okay), 'actual': actual,
                       'expected': expected, 'absolute_tolerance': tolerance})

    pyrosetta = _initialize([str(params)])
    pose = pyrosetta.pose_from_file(str(mutant))
    sf = _scorefunction(pyrosetta, weights)
    total = float(sf(pose))
    terms, sf_weights = _active_terms(sf), sf.weights()
    pair_array = complete_pair_energies_array(pose, sf)
    pair_total = float(np.triu(pair_array['total_score'], 1).sum())
    onebody = sum(float(pose.energies().onebody_energies(i)[term]) * float(sf_weights[term])
                  for i in range(1, pose.size() + 1) for term in terms)
    check('full_score_equals_onebody_plus_unique_pair_energies', onebody + pair_total, total)
    write_pair_energies(pose, sf, output / 'mutant_pair_energies.sc', mutant.stem)
    write_interface_metrics(pose, sf, 'X', output / 'mutant_interface.sc', mutant.stem)
    write_pose_hbonds(pose, output / 'mutant_hbonds.csv', sf)
    weighted_hbonds = sum(float(row['hbond_weighted_energy_REU']) for row in _csv(output/'mutant_hbonds.csv'))
    scored_hbond_terms = sum(float(pose.energies().total_energies()[term])*float(sf_weights[term])
                             for term in terms if _term_name(term).startswith('hbond_'))
    check('individual_weighted_hbonds_reconstruct_scored_hbond_terms', weighted_hbonds, scored_hbond_terms, 1e-6)
    score_lines = (output/'mutant_pair_energies.sc').read_text().splitlines()
    header = score_lines[0].split()[1:]
    records = [dict(zip(header,line.split()[1:])) for line in score_lines[1:]]
    check('exported_pair_plus_self_table_reconstructs_bound_total', sum(float(r['total']) for r in records), total, 0.002)
    fresh = read_scorefile_metrics(output / 'mutant_interface.sc')
    old = read_scorefile_metrics(result_dir / 'mutant_interface.sc')
    check('fresh_interface_equals_saved_interface', fresh['interface_delta_X'], old['interface_delta_X'], 2e-5)
    check('fresh_bound_score_equals_saved_bound_score', total, old['bound_total'], 2e-5)
    term_delta = sum(value for key, value in fresh.items() if key.startswith('if_X_'))
    check('interface_equals_sum_of_weighted_term_differences', term_delta, fresh['interface_delta_X'], 2e-5)

    # Do not use the production separation method for this oracle.
    separated = pose.clone()
    from pyrosetta.rosetta.core.id import AtomID
    from pyrosetta.rosetta.numeric import xyzVector_double_t
    ligand_indices = [i for i in range(1, pose.size() + 1) if pose.pdb_info().chain(i) == 'X']
    for i in ligand_indices:
        for a in range(1, separated.residue(i).natoms() + 1):
            separated.set_xyz(AtomID(a, i), separated.residue(i).xyz(a) + xyzVector_double_t(0, 2000, 0))
    check('interface_equals_independent_cartesian_separation', total - float(sf(separated)), fresh['interface_delta_X'], 2e-5)

    from pyrosetta.rosetta.core.scoring import EMapVector
    maximum_error, evaluated = 0.0, 0
    short_array = pose.energies().residue_pair_energies_array()
    for i in range(1, pose.size() + 1):
        for j in range(i + 1, pose.size() + 1):
            if abs(float(pair_array[i-1, j-1]['total_score'])) < 1e-12:
                continue
            emap = EMapVector()
            sf.eval_ci_2b(pose.residue(i), pose.residue(j), pose, emap)
            sf.eval_cd_2b(pose.residue(i), pose.residue(j), pose, emap)
            for term in terms:
                name = _term_name(term)
                if name in short_array.dtype.names:
                    maximum_error = max(maximum_error, abs(float(emap[term]) - float(short_array[i-1, j-1][name])))
            evaluated += 1
    check('all_nonzero_pairs_match_direct_two_body_evaluation', maximum_error, 0.0)

    structure = read_pdb(mutant)
    check('coordinate_residue_count_equals_rosetta_pose', len(structure.residues), pose.size(), 0)
    maximum_chi_error, chi_count, atom_shift, rebuilt_hydrogens = 0.0, 0, 0.0, 0
    for i, residue in enumerate(structure.residues, 1):
        check(f'pose_mapping_{i}', f'{pose.pdb_info().chain(i)}:{pose.pdb_info().number(i)}', f'{residue.key.chain}:{residue.key.number}')
        for k, angle in enumerate(chi_angles(residue), 1):
            difference = (angle - float(pose.chi(k, i)) + 180) % 360 - 180
            maximum_chi_error = max(maximum_chi_error, abs(difference))
            chi_count += 1
        for a in range(1, pose.residue(i).natoms() + 1):
            atom = residue.atom(pose.residue(i).atom_name(a).strip())
            if atom is not None:
                xyz = pose.residue(i).xyz(a)
                atom_shift = max(atom_shift, float(np.linalg.norm(atom.coord - [xyz.x, xyz.y, xyz.z])))
            elif pose.residue(i).atom_is_hydrogen(a):
                rebuilt_hydrogens += 1
    check('all_available_canonical_chi_angles_match_rosetta', maximum_chi_error, 0.0, 1e-5)
    check('rosetta_loading_preserves_supplied_atom_coordinates', atom_shift, 0.0, 1e-5)
    if rebuilt_hydrogens:
        warnings.append(f'Rosetta rebuilt {rebuilt_hydrogens} missing hydrogens. Geometric tools reading the original prepared PDB do not see those added atoms.')
    ligand_validation = validate_ligand_against_params(mutant, params, 'X')
    check('ligand_atom_and_bond_validation', ligand_validation['valid'], True)
    params_charge = sum(float(line.split()[-1]) for line in params.read_text().splitlines() if line.startswith('ATOM '))
    check('estradiol_params_net_charge', params_charge, 0.0, 0.01)

    from rdkit import Chem
    from fv_energy_analysis.resources import resource_path
    import importlib.util
    spec = importlib.util.spec_from_file_location('_sanity_ligand_prepare', resource_path('scripts/prepare_ligand_complexes.py'))
    prep = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(prep)
    actual_molecule = prep.ligand_fragment(mutant, 'X')
    reference_molecule = Chem.MolFromPDBFile(str(template), removeHs=False)
    from fv_energy_analysis.ligand_chemistry import reference_smiles_from_source, assign_reference_chemistry
    chemical_reference = reference_smiles_from_source(params=params)
    if chemical_reference:
        actual_molecule = assign_reference_chemistry(actual_molecule, chemical_reference)
        reference_molecule = assign_reference_chemistry(reference_molecule, chemical_reference)
    stereo = []
    for molecule in (actual_molecule, reference_molecule):
        Chem.AssignAtomChiralTagsFromStructure(molecule, replaceExistingTags=True)
        Chem.AssignStereochemistry(molecule, cleanIt=True, force=True)
        stereo.append({a.GetPDBResidueInfo().GetName().strip(): a.GetProp('_CIPCode')
                       for a in molecule.GetAtoms() if a.HasProp('_CIPCode')})
    check('ligand_stereochemistry_matches_parameter_template', stereo[0], stereo[1])

    p_structure = read_pdb(parent)
    chains = metadata['fv_chains']
    pairs = _validated_pairs(p_structure, structure, chains)
    annotations = annotate_fv(p_structure, chains, metadata['numbering_scheme'], metadata['cdr_definition'])
    regions = {key: region for ann in annotations.values() for key, region in ann.region_by_key.items()}
    transform, rmsds = structural_rmsds(p_structure, structure, pairs, regions)
    all_p = [p.key for pair in pairs.values() for p in pair[0]]
    all_m = [m.key for pair in pairs.values() for m in pair[1]]
    px, mx, _ = paired_atoms(p_structure, structure, all_p, all_m)
    from Bio.SVDSuperimposer import SVDSuperimposer
    oracle = SVDSuperimposer()
    oracle.set(px, mx)
    oracle.run()
    rotation, translation = oracle.get_rotran()
    check('global_transform_matches_independent_biopython', float(np.max(np.abs(transform.apply(mx) - (mx @ rotation + translation)))), 0.0, 1e-8)
    by_region = {r['region']: r for r in rmsds}
    whole, anchors, tip = [by_region[name] for name in ('H_CDR3', 'H_CDR3_ANCHORS', 'H_CDR3_TIP')]
    check('unique_parent_YYGSS_tip_has_five_residues', tip['n_residues'], 5, 0)
    partition = sum(r['n_atoms'] * float(r['rmsd_after_global_fit_A']) ** 2 for r in (anchors, tip))
    check('overall_hcdr3_squared_displacement_partitions_into_anchors_and_tip', partition, whole['n_atoms'] * float(whole['rmsd_after_global_fit_A']) ** 2, 1e-8)

    compare(SimpleNamespace(parent=str(parent), mutant=str(mutant), fv_chains=chains, ligand_chains=['X'],
        scheme=metadata['numbering_scheme'], cdr_definition=metadata['cdr_definition'], fit_region='fv',
        rmsd_atoms='backbone', rotamer_rmsd_cutoff=3.0, hcdr3_tip_sequence='YYGSS',
        network_energy_presence_cutoff=0.1, network_energy_change_cutoff=0.25,
        parent_rosetta=str(result_dir/'parent_pair_energies.sc'), mutant_rosetta=str(output/'mutant_pair_energies.sc'),
        parent_interface=str(result_dir/'parent_interface.sc'), mutant_interface=str(output/'mutant_interface.sc'),
        parent_arpeggio=str(result_dir/'arpeggio_parent'/'parent_arpeggio.json'),
        mutant_arpeggio=str(result_dir/'arpeggio_mutant'/'mutant_arpeggio.json'),
        parent_hbonds=str(result_dir/'parent_hbonds.csv'), mutant_hbonds=str(output/'mutant_hbonds.csv'),
        output=str(output/'final')))
    pair_rows = _csv(output/'final'/'pair_comparison.csv')
    check('pair_energy_deltas_are_mutant_minus_parent', max(abs(float(r.get('delta_total') or 0) -
          (float(r.get('mutant_total') or 0) - float(r.get('parent_total') or 0))) for r in pair_rows), 0.0, 1e-12)
    from fv_energy_analysis.plip import summarize_plip_xml, INTERACTION_TAGS
    xml_path = result_dir/'plip_mutant'/'mutant_plip.xml'
    if xml_path.is_file():
        xml = ET.parse(xml_path).getroot()
        plip = summarize_plip_xml(xml_path)
        for group, tag in INTERACTION_TAGS.items():
            expected = sum(len(site.findall(f'interactions/{group}/{tag}')) for site in xml.findall('bindingsite') if site.findtext('identifiers/hetid') == 'E2L')
            check(f'plip_raw_xml_count_{group}', plip['interaction_counts'][group], expected, 0)
    else:
        warnings.append('Saved PLIP XML is not available; summary counts were not independently verified against XML.')
    from fv_energy_analysis.workflow import _write_automatic_summary
    (output/'ligand_atom_validation.json').write_text(json.dumps({'mutant': ligand_validation}) + '\n')
    _write_automatic_summary(output, architecture, old_summary.get('plip_ligand_interactions'))
    report = {'status': 'passed' if all(c['passed'] for c in checks) else 'failed',
              'candidate': str(mutant), 'checks': checks, 'warnings': warnings,
              'evaluated_pairs': evaluated, 'checked_chi_angles': chi_count, 'rebuilt_hydrogens': rebuilt_hydrogens,
              'interface_REU': fresh['interface_delta_X'], 'hcdr3_rmsds': [whole, anchors, tip],
              'scope': 'One prepared PDB plus its saved parent comparison; not experimental affinity or all-input validation'}
    write_json(output/'validation.json', report)
    lines = ['# Single-PDB independent validation', '', f"Status: {report['status']}", '', f'Candidate: `{mutant}`', '',
             '| Check | Result | Actual | Expected |', '|---|---|---|---|']
    lines.extend(f"| {c['check']} | {'PASS' if c['passed'] else 'FAIL'} | {c['actual']} | {c['expected']} |" for c in checks if not c['check'].startswith('pose_mapping_'))
    lines += ['', f'All {pose.size()} residue mappings checked; {evaluated} nonzero residue pairs; {chi_count} chi angles.', '', '## HCDR3 after the same global fit', '', '| Region | Residues | RMSD (A) |', '|---|---|---:|']
    lines.extend(f"| {r['region']} | {r['parent_residues']} | {float(r['rmsd_after_global_fit_A']):.6f} |" for r in (whole, anchors, tip))
    lines += ['', '## Limitations and warnings', '', *[f'- {w}' for w in warnings], '- This numerical validation does not establish affinity accuracy, protonation equilibria, or significance of static-pose score differences.']
    (output/'validation.md').write_text('\n'.join(lines) + '\n')
    return report


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--result-dir', type=Path, required=True)
    parser.add_argument('--params', type=Path, required=True)
    parser.add_argument('--ligand-template-pdb', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--weights', default='ref2015')
    args = parser.parse_args(argv)
    result = run_sanity(args.result_dir, args.params, args.ligand_template_pdb, args.output, args.weights)
    print(json.dumps({key: result[key] for key in ('status', 'evaluated_pairs', 'checked_chi_angles', 'interface_REU', 'warnings')}))
    return 0 if result['status'] == 'passed' else 1


if __name__ == '__main__':
    raise SystemExit(main())
