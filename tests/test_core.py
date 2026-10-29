from __future__ import annotations

import json
from types import SimpleNamespace

import numpy as np

from fv_energy_analysis.arpeggio import read_arpeggio
from fv_energy_analysis.batch_workflow import (
    _candidate_network_summary,
    _hbond_saltbridge_summary,
    _parent_cache_key,
    _result_names,
    _write_master_reports,
    build_parser as build_batch_parser,
)
from fv_energy_analysis.chemistry import substitution_chemistry
# Historical algorithm stays explicitly testable, never selected by the public
# DS-only API. Current feature tests live in test_ds_architecture.py.
from fv_energy_analysis.estrogen_architecture import _analyze_legacy_estrogen_architecture as analyze_estrogen_architecture
from fv_energy_analysis.discovery_studio import _normalize, _summary
from fv_energy_analysis.fast_workflow import _link_or_copy
from fv_energy_analysis.fingerprint_counts import _arpeggio_counts
from fv_energy_analysis.fingerprint_reports import _flatten
from fv_energy_analysis.geometry import kabsch, rmsd
from fv_energy_analysis.hbonds import aggregate_hbond_csv
from fv_energy_analysis.models import Residue, ResidueKey
from fv_energy_analysis.plip import summarize_plip_xml
from fv_energy_analysis.network import mutation_network_changes
from fv_energy_analysis.rotamers import angle_bin, compare_rotamers
from fv_energy_analysis.rosetta import read_residue_energy_breakdown
from fv_energy_analysis.arpeggio import run_arpeggio
from fv_energy_analysis.structure import read_pdb
from fv_energy_analysis.workflow import (
    _chain_map,
    _convert_structure_input,
    _mutation_interpretations,
    _residue_contact_present,
)
from fv_energy_analysis.validation import validate_ligand_against_params


def test_kabsch_recovers_rigid_transform():
    reference = np.array([[0., 0., 0.], [1., 0., 0.], [0., 2., 0.], [0., 0., 3.]])
    rotation = np.array([[0., -1., 0.], [1., 0., 0.], [0., 0., 1.]])
    mobile = reference @ rotation + np.array([5., -2., 1.])
    fit = kabsch(mobile, reference)
    assert rmsd(fit.apply(mobile), reference) < 1e-10


def test_plip_xml_summary(tmp_path):
    report = tmp_path / "report.xml"
    report.write_text(
        "<report><plipversion>3.0.1</plipversion><bindingsite><identifiers>"
        "<hetid>E2L</hetid></identifiers><interactions><hydrogen_bonds>"
        "<hydrogen_bond><reschain>H</reschain><resnr>35</resnr><restype>SER</restype>"
        "</hydrogen_bond></hydrogen_bonds><pi_stacks><pi_stack><reschain>L</reschain>"
        "<resnr>91</resnr><restype>TYR</restype></pi_stack></pi_stacks>"
        "</interactions></bindingsite></report>",
        encoding="utf-8",
    )
    summary = summarize_plip_xml(report)
    assert summary["binding_site_count"] == 1
    assert summary["interaction_counts"]["hydrogen_bonds"] == 1
    assert summary["interaction_counts"]["pi_stacks"] == 1
    assert summary["residues_by_type"]["pi_stacks"] == ["L:91:TYR"]


def test_pdb_chain_remap_preserves_conect(tmp_path):
    source = tmp_path / "input.pdb"
    destination = tmp_path / "mapped.pdb"
    source.write_text(
        "ATOM      1  CA  ALA H   1       0.000   0.000   0.000  1.00  0.00           C  \n"
        "HETATM    2  C1  LIG d   1       1.000   0.000   0.000  1.00  0.00           C  \n"
        "CONECT    2    3\n"
        "END\n",
        encoding="utf-8",
    )
    result = _convert_structure_input(source, destination, {"H": "L"}, "d")
    mapped = result.read_text(encoding="utf-8")
    assert mapped.splitlines()[0][21] == "L"
    assert mapped.splitlines()[1][21] == "d"
    assert "CONECT    2    3" in mapped


def test_chemical_classes():
    conservative = substitution_chemistry("L", "I")
    assert conservative["same_primary_chemical_class"] is True
    radical = substitution_chemistry("D", "K")
    assert radical["same_primary_chemical_class"] is False
    assert "negative" in radical["lost_properties"]


def test_rotamer_bins_and_noncomparable_full_bin():
    assert angle_bin(60) == "p"
    assert angle_bin(-60) == "m"
    assert angle_bin(179) == "t"
    ala = Residue(ResidueKey("H", 1), "ALA")
    gly = Residue(ResidueKey("H", 1), "GLY")
    result = compare_rotamers(ala, gly)
    assert result["same_full_rotamer_bin"] is None


def test_arpeggio_aggregates_contacts(tmp_path):
    endpoint_a = {"auth_asym_id": "H", "auth_seq_id": 31, "pdbx_PDB_ins_code": " ", "label_comp_id": "TYR"}
    endpoint_b = {"auth_asym_id": "X", "auth_seq_id": 1, "pdbx_PDB_ins_code": " ", "label_comp_id": "LIG"}
    data = [
        {"bgn": endpoint_a, "end": endpoint_b, "contact": ["hbond", "vdw"], "distance": 3.1, "type": "atom-atom"},
        {"bgn": endpoint_b, "end": endpoint_a, "contact": ["FT"], "distance": 4.8, "type": "plane-plane"},
    ]
    path = tmp_path / "contacts.json"
    path.write_text(json.dumps(data))
    rows = read_arpeggio(path, {"H", "L"})
    assert len(rows) == 1
    assert rows[0]["count_hbond"] == 1
    assert rows[0]["count_pi_stack"] == 1


def test_estrogen_architecture_detects_strict_pi_stack_sandwich(tmp_path):
    import pytest
    def pdb_atom(record, serial, name, resname, chain, number, xyz, element):
        x, y, z = xyz
        return (
            f"{record:<6}{serial:5d} {name:^4} {resname:>3} {chain}{number:4d}    "
            f"{x:8.3f}{y:8.3f}{z:8.3f}  1.00 20.00          {element:>2}\n"
        )

    lines = []
    serial = 1
    ring_names = ("C10", "C11", "C12", "C13", "C14", "C15")
    ring_xyz = ((1.4, 0, 0), (0.7, 1.2, 0), (-0.7, 1.2, 0),
                (-1.4, 0, 0), (-0.7, -1.2, 0), (0.7, -1.2, 0))
    for name, xyz in zip(ring_names, ring_xyz):
        lines.append(pdb_atom("HETATM", serial, name, "E2L", "X", 1, xyz, "C"))
        serial += 1
    for name, xyz, element in (("O1", (2.5, 0, 0), "O"), ("O2", (-3, 0, 0), "O"),
                               ("C5", (0, 2.5, 0), "C")):
        lines.append(pdb_atom("HETATM", serial, name, "E2L", "X", 1, xyz, element))
        serial += 1
    aromatic_names = ("CG", "CD1", "CD2", "CE1", "CE2", "CZ")
    for chain, number, z in (("H", 99, 3.4), ("L", 91, -3.4)):
        for index, name in enumerate(aromatic_names):
            lines.append(pdb_atom("ATOM", serial, name, "PHE", chain, number,
                                  (ring_xyz[index][0], ring_xyz[index][1], z), "C"))
            serial += 1
    pdb = tmp_path / "complex.pdb"
    pdb.write_text("".join(lines) + "END\n", encoding="utf-8")

    ligand = {
        "auth_asym_id": "X", "auth_seq_id": 1, "pdbx_PDB_ins_code": " ",
        "label_comp_id": "E2L",
    }
    def residue(chain, number, name, atoms):
        return {
            "auth_asym_id": chain, "auth_seq_id": number, "pdbx_PDB_ins_code": " ",
            "label_comp_id": name, "auth_atom_id": atoms,
        }
    interactions = [
        {"bgn": residue("H", 31, "TYR", "OH"), "end": {**ligand, "auth_atom_id": "O1"},
         "contact": ["hbond"], "type": "atom-atom"},
        {"bgn": residue("L", 35, "GLN", "NE2"), "end": {**ligand, "auth_atom_id": "O2"},
         "contact": ["hbond"], "type": "atom-atom"},
        {"bgn": residue("H", 50, "LEU", "CD1"), "end": {**ligand, "auth_atom_id": "C5"},
         "contact": ["hydrophobic"], "type": "atom-atom"},
        {"bgn": residue("H", 99, "PHE", ",".join(aromatic_names)),
         "end": {**ligand, "auth_atom_id": ",".join(ring_names)},
         "contact": ["FF"], "type": "plane-plane"},
        {"bgn": {**ligand, "auth_atom_id": ",".join(ring_names)},
         "end": residue("L", 91, "PHE", ",".join(aromatic_names)),
         "contact": ["FF"], "type": "plane-plane"},
    ]
    contacts = tmp_path / "contacts.json"
    contacts.write_text(json.dumps(interactions), encoding="utf-8")
    import csv
    ds = tmp_path / "ds.csv"
    with ds.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["interaction_type", "protein_atoms", "ligand_atoms"])
        writer.writeheader()
        for chain, number in (("H", 99), ("L", 91)):
            writer.writerow({"interaction_type": "pi_pi_stacked",
                             "protein_atoms": ";".join(f"{chain}:PHE{number}:{name}" for name in aromatic_names),
                             "ligand_atoms": ";".join(f"X:E2L1:{name}" for name in ring_names)})
    kwargs = dict(fv_chains={"H", "L"}, o3_atom="O1", o17_atom="O2",
                  c18_methyl_atom="C5", a_ring_atoms=ring_names)
    result = analyze_estrogen_architecture(pdb, contacts, discovery_studio_csv=ds, **kwargs)
    assert result["e2_conserved_architecture_score_0_to_4"] == 4
    assert result["e2_a_ring_pi_stack_count"] == 2
    assert result["e2_pi_stack_sandwich"] is True
    assert result["e2_strict_screen_pass"] is True
    # Arpeggio plane records are not required for DS-defined stacks.
    contacts.write_text(json.dumps(interactions[:3]))
    assert analyze_estrogen_architecture(pdb, contacts, discovery_studio_csv=ds, **kwargs)["e2_pi_stack_sandwich"] is True
    missing = analyze_estrogen_architecture(pdb, contacts, **kwargs)
    assert missing["e2_a_ring_pi_stack_count"] is None
    assert missing["e2_strict_screen_pass"] is None
    filtered = analyze_estrogen_architecture(pdb, contacts, discovery_studio_csv=ds, **{**kwargs, "fv_chains": {"H"}})
    assert filtered["e2_o17_hbond"] is False
    assert filtered["e2_a_ring_pi_stack_count"] == 1
    assert filtered["e2_pi_stack_sandwich"] is False
    with ds.open("a", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["pi_pi_stacked", ";".join(f"H:PHE99:{n}" for n in aromatic_names[:3]),
                         ";".join(f"X:E2L1:{n}" for n in ring_names)])
    multiple = analyze_estrogen_architecture(pdb, contacts, discovery_studio_csv=ds, **kwargs)
    assert multiple["e2_a_ring_pi_stack_ring_count"] == 3
    assert multiple["e2_a_ring_pi_stack_count"] == 2
    ds.write_text(ds.read_text().replace("pi_pi_stacked", "pi_pi_t_shaped"))
    tshaped = analyze_estrogen_architecture(pdb, contacts, discovery_studio_csv=ds, **kwargs)
    assert tshaped["e2_a_ring_pi_stack_count"] == 0
    assert tshaped["e2_ds_pi_t_stack_partner_count"] == 2
    assert tshaped["e2_pi_stack_sandwich"] is False
    assert tshaped["e2_a_ring_pi_contact_count"] == 2
    assert tshaped["e2_pi_contact_sandwich"] is True
    assert tshaped["e2_strict_screen_pass"] is True
    for ds_type, metric in (("pi_cation", "pi_cation"), ("pi_anion", "pi_anion")):
        charged = tmp_path / f"{ds_type}.csv"
        charged.write_text(ds.read_text().replace("pi_pi_t_shaped", ds_type))
        result = analyze_estrogen_architecture(pdb, contacts, discovery_studio_csv=charged, **kwargs)
        assert result[f"e2_ds_{metric}_partner_count"] == 2
        assert result["e2_a_ring_pi_stack_count"] == 0
        assert result["e2_a_ring_pi_contact_count"] == 2
        assert result["e2_pi_contact_sandwich"] is True
    # Charge partners need not be aromatic or contain enough atoms for a plane.
    pdb.write_text(pdb.read_text().replace("END\n", "")
                   + pdb_atom("ATOM", 900, "NZ", "LYS", "H", 55, (0, 0, 3.4), "N")
                   + pdb_atom("ATOM", 901, "OD1", "ASP", "L", 56, (0, 0, -3.4), "O") + "END\n")
    with ds.open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["interaction_type", "protein_atoms", "ligand_atoms"])
        ring_group = ";".join(f"X:E2L1:{name}" for name in ring_names)
        writer.writerow(["pi_cation", "H:LYS55:NZ", ring_group])
        writer.writerow(["pi_anion", "L:ASP56:OD1", ring_group])
        writer.writerow(["pi_cation", "H:LYS55:NZ", ring_group])
    mixed = analyze_estrogen_architecture(pdb, contacts, discovery_studio_csv=ds, **kwargs)
    assert mixed["e2_a_ring_pi_contact_count"] == 2
    assert mixed["e2_pi_contact_sandwich"] is True
    assert mixed["e2_strict_screen_pass"] is True
    with pytest.raises(ValueError, match="no legacy"):
        analyze_estrogen_architecture(pdb, contacts)


def test_discovery_studio_normalizes_detector_specific_interactions():
    rows = _normalize([
        {
            "structure_id": "m0001", "role": "mutant", "candidate": "x.pdb",
            "result_folder": "x", "interaction_type": "piAnionType",
            "interaction_name": "Pi-Anion", "interaction_categories": "allChargeNonbondTypes",
            "category_names": "Electrostatic", "distance_A": "4.326",
            "from_atoms": "H:ASP112:OD1",
            "to_atoms": "X:E2L1:C10;X:E2L1:C11;X:E2L1:C12",
            "from_role": "receptor", "to_role": "ligand",
        }
    ])
    assert rows[0]["interaction_type"] == "pi_anion"
    assert rows[0]["protein_residue"] == "H:112:ASP"
    assert rows[0]["residue_pair_id"] == "H:112:ASP<->X:1:E2L"
    summary = _summary(rows, "mutant")
    assert summary["mutant_ds_pi_anion_annotations"] == 1
    assert summary["mutant_ds_pi_anion_unique_residues"] == 1


def test_fingerprint_summary_is_flat_and_omits_unsupported_columns(tmp_path):
    arpeggio = tmp_path / "contacts.json"
    arpeggio.write_text(json.dumps([{
        "bgn": {"label_comp_id": "SER", "label_comp_type": "P"},
        "end": {"label_comp_id": "E2L", "label_comp_type": "B"},
        "contact": ["hbond", "vdw", "CATIONPI"],
    }]))
    counts = _arpeggio_counts(arpeggio, "E2L")
    assert counts["hbond"] == 1
    assert counts["vdw"] == 1
    assert counts["CATIONPI"] == 1
    assert counts["carbon_hbond"] == "NA"
    wide = {"candidate": "x.pdb", "result_folder": str(tmp_path), "tool": "Arpeggio", "count_unit": "annotations"}
    for name, value in counts.items():
        wide[f"parent_{name}"] = value
        wide[f"mutant_{name}"] = value
        wide[f"delta_{name}"] = 0 if isinstance(value, int) else "NA"
    flat = _flatten([wide])
    assert flat["parent_arpeggio_hbond"] == 1
    assert flat["mutant_arpeggio_CATIONPI"] == 1
    assert "parent_arpeggio_carbon_hbond" not in flat


def test_arpeggio_ring_contacts_are_not_aromatic_atom_labels(tmp_path):
    endpoints = {
        "bgn": {"label_comp_id": "PHE", "label_comp_type": "P"},
        "end": {"label_comp_id": "E2L", "label_comp_type": "B"},
    }
    path = tmp_path / "contacts.json"
    path.write_text(json.dumps([
        {**endpoints, "type": "atom-atom", "contact": ["aromatic"]},
        {**endpoints, "type": "atom-atom", "contact": ["aromatic", "CARBONPI"]},
        {**endpoints, "type": "plane-plane", "contact": ["FT", "OT"]},
    ]))
    counts = _arpeggio_counts(path, "E2L")
    assert counts["aromatic"] == 2
    assert counts["plane-plane"] == 1
    assert counts["CARBONPI"] == 1
    assert counts["pi_anion"] == "NA"


def test_rosetta_pair_mapping():
    structure = read_pdb("tests/fixtures/minimal_complex.pdb")
    rows = read_residue_energy_breakdown("tests/fixtures/energy.sc", structure)
    assert rows[0]["pair_id"] == "H:1--X:1"
    assert rows[0]["hbond_total"] == -0.6
    assert rows[0]["vdw_net"] == -1.3


def test_mutation_network_classifies_formed_broken_and_type_changes():
    substitutions = [{
        "substitution": "E1Y", "chain": "H", "parent_antibody_number": "H1",
        "parent_region": "CDR1", "parent_residue": "H:1", "mutant_residue": "H:1",
    }]
    comparison = [
        {
            "parent_residue1": "H:1", "parent_residue2": "L:2",
            "mutant_residue1": "H:1", "mutant_residue2": "L:2",
            "parent_interaction_types": "hbond;polar", "mutant_interaction_types": "hydrophobic",
            "parent_count_hbond": 1, "parent_count_polar": 1, "mutant_count_hydrophobic": 2,
            "parent_total": -1.5, "mutant_total": -0.5,
            "parent_pair_energy_strength": "moderate_favorable",
            "mutant_pair_energy_strength": "weak_favorable", "mutant_pair_scope": "aa_to_aa",
        },
        {
            "parent_residue1": "H:1", "parent_residue2": "L:3",
            "mutant_residue1": "", "mutant_residue2": "",
            "parent_interaction_types": "ionic", "mutant_interaction_types": "",
            "parent_count_ionic": 1, "parent_total": -2.0, "mutant_total": 0.0,
            "parent_pair_scope": "aa_to_aa",
        },
        {
            "parent_residue1": "", "parent_residue2": "",
            "mutant_residue1": "H:1", "mutant_residue2": "X:1",
            "parent_interaction_types": "", "mutant_interaction_types": "pi_stack",
            "mutant_count_pi_stack": 1, "parent_total": 0.0, "mutant_total": -1.2,
            "mutant_pair_scope": "aa_to_ligand",
        },
    ]
    details, summaries, type_energies = mutation_network_changes(comparison, substitutions)
    assert {row["network_status"] for row in details} == {"formed", "broken", "retained_changed"}
    changed = next(row for row in details if row["network_status"] == "retained_changed")
    assert changed["broken_interaction_types"] == "hbond;polar"
    assert changed["formed_interaction_types"] == "hydrophobic"
    assert summaries[0]["formed_pairs"] == 1
    assert summaries[0]["broken_pairs"] == 1
    formed = next(row for row in details if row["network_status"] == "formed")
    assert formed["event_basis"] == "typed_contact_or_hbond"
    pi_stack = next(row for row in type_energies if row["interaction_type"] == "pi_stack")
    assert pi_stack["contributing_pairs"] == "H:1<->X:1"


def test_mutation_interpretation_distinguishes_ligand_gain_from_local_cost():
    substitutions = [{
        "substitution": "Y1H", "parent_antibody_number": "L1", "parent_region": "CDR3",
        "parent_chemical_class": "aromatic", "mutant_chemical_class": "basic_or_titratable",
        "same_primary_chemical_class": "False", "rotamer_analysis_status": "performed",
        "same_full_rotamer_bin": "True", "same_common_chi_bins": "True",
        "n_comparable_chi": "2", "cdr_global_fit_rmsd_A": "0.7",
    }]
    effects = [{
        "substitution": "Y1H", "direct_ligand_delta_REU": -0.8,
        "direct_ligand_parent_REU": -1.0, "direct_ligand_mutant_REU": -1.8,
        "parent_contact_annotations": 3.0, "mutant_contact_annotations": 5.0,
    }]
    summaries = [{"mutation": "Y1H", "sum_pair_energy_delta_REU": "2.0"}]
    changes = [{
        "mutation": "Y1H", "mutant_partner": "X:1", "parent_partner": "X:1",
        "mutant_interaction_types": "hbond", "parent_interaction_types": "vdw",
        "delta_total": "-0.8",
    }, {
        "mutation": "Y1H", "mutant_partner": "L:2", "parent_partner": "L:2",
        "mutant_interaction_types": "energy_only", "parent_interaction_types": "pi_stack",
        "delta_total": "2.8",
    }]
    result = _mutation_interpretations(substitutions, effects, summaries, changes, -1.5)[0]
    assert result["assessment"] == "Mixed: ligand interaction improves, local network worsens"
    assert "X:1" in result["strongest_supporting_pair_changes"]
    assert "L:2" in result["strongest_opposing_pair_changes"]


def test_chain_map_parses_swapped_chains():
    assert _chain_map("H:L,L:H") == {"H": "L", "L": "H"}


def test_residue_contact_presence_collapses_atom_annotations():
    row = {
        "parent_total_contact_count": "14",
        "mutant_total_contact_count": "0",
        "parent_count_rosetta_hbond": "0",
        "mutant_count_rosetta_hbond": "1",
    }
    assert _residue_contact_present(row, "parent") is True
    assert _residue_contact_present(row, "mutant") is True
    row["mutant_count_rosetta_hbond"] = "0"
    assert _residue_contact_present(row, "mutant") is False


def test_mutation_interpretation_detects_local_gain_with_ligand_cost():
    substitutions = [{
        "substitution": "E1F", "parent_antibody_number": "H1", "parent_region": "CDR3",
        "parent_chemical_class": "acidic_negative", "mutant_chemical_class": "aromatic",
        "same_primary_chemical_class": "False", "rotamer_analysis_status": "skipped",
        "n_comparable_chi": "0", "cdr_global_fit_rmsd_A": "4.0",
    }]
    effects = [{
        "substitution": "E1F", "direct_ligand_delta_REU": 0.5,
        "direct_ligand_parent_REU": -3.0, "direct_ligand_mutant_REU": -2.5,
        "parent_contact_annotations": 10.0, "mutant_contact_annotations": 20.0,
    }]
    result = _mutation_interpretations(
        substitutions, effects,
        [{"mutation": "E1F", "sum_pair_energy_delta_REU": "-5.0"}], [], -2.0,
    )[0]
    assert result["assessment"] == "Mixed: local network improves, direct ligand interaction worsens"


def test_batch_report_ranks_completed_complexes(tmp_path):
    first = tmp_path / "first"
    second = tmp_path / "second"
    first.mkdir()
    second.mkdir()
    rows = [
        {
            "candidate": "first.pdb", "result_folder": str(first), "status": "complete",
            "substitutions": "A1V", "interface_delta_REU": 1.0, "bound_total_delta_REU": 2.0,
            "global_fv_rmsd_A": 1.0, "h_cdr3_rmsd_A": 1.0, "l_cdr3_rmsd_A": 1.0,
            "parent_rosetta_ligand_hbonds": 1, "mutant_rosetta_ligand_hbonds": 1,
            "parent_arpeggio_ligand_hbonds": 1, "mutant_arpeggio_ligand_hbonds": 1,
            "comparison_quality_warning": False,
        },
        {
            "candidate": "second.pdb", "result_folder": str(second), "status": "complete",
            "substitutions": "A1F", "interface_delta_REU": -2.0, "bound_total_delta_REU": -1.0,
            "global_fv_rmsd_A": 1.0, "h_cdr3_rmsd_A": 1.0, "l_cdr3_rmsd_A": 1.0,
            "parent_rosetta_ligand_hbonds": 1, "mutant_rosetta_ligand_hbonds": 2,
            "parent_arpeggio_ligand_hbonds": 1, "mutant_arpeggio_ligand_hbonds": 2,
            "comparison_quality_warning": False,
        },
    ]
    _write_master_reports(tmp_path, rows, [])
    assert rows[1]["rank_by_interface_delta"] == 1
    assert "second.pdb" in (tmp_path / "batch_summary.md").read_text()
    duplicate_paths = [tmp_path / "a" / "same.pdb", tmp_path / "b" / "same.pdb"]
    assert len(set(_result_names(duplicate_paths, tmp_path).values())) == 2


def test_optimized_batch_defaults_to_two_persistent_workers():
    args = build_batch_parser().parse_args([
        "--parent", "parent.pdb", "--input-folder", "inputs", "--output", "results",
    ])
    assert args.jobs == 2


def test_parent_cache_key_changes_with_scientific_settings(tmp_path):
    parent = tmp_path / "parent.pdb"
    params = tmp_path / "lig.params"
    template = tmp_path / "lig.pdb"
    parent.write_text("END\n")
    params.write_text("NAME LIG\n")
    template.write_text("END\n")
    base = dict(
        fv_chains="H,L", ligand_chain="d", parent_ligand_chain="d",
        parent_chain_map={}, weights="ref2015", params=str(params),
        ligand_template_pdb=str(template), scheme="imgt", cdr_definition="imgt",
    )
    first = _parent_cache_key(SimpleNamespace(**base), parent)
    second = _parent_cache_key(SimpleNamespace(**{**base, "weights": "beta_nov16"}), parent)
    assert first != second


def test_parent_artifact_link_or_copy_is_reusable(tmp_path):
    source = tmp_path / "source.sc"
    destination = tmp_path / "child" / "parent.sc"
    source.write_text("score\n")
    _link_or_copy(source, destination)
    assert destination.read_text() == "score\n"
    source.write_text("updated\n")
    _link_or_copy(source, destination)
    assert destination.read_text() == "updated\n"


def test_arpeggio_runs_once_with_fixed_seed(tmp_path):
    from unittest.mock import Mock
    runner = Mock()
    cif = tmp_path / "complex.cif"
    directory = tmp_path / "arpeggio"
    output = run_arpeggio(runner, "pdbe-arpeggio", directory, cif)
    runner.run.assert_called_once_with(
        ["pdbe-arpeggio", "-o", str(directory), "-m", str(cif)],
        env_overrides={"PYTHONHASHSEED": "0"},
    )
    assert output == directory / "complex.json"
    assert directory.is_dir()


def test_full_complex_contacts_keep_protein_pairs_but_not_in_ligand_fingerprints(tmp_path):
    def endpoint(chain, number, name, kind):
        return {"auth_asym_id": chain, "auth_seq_id": number,
                "label_comp_id": name, "label_comp_type": kind}
    path = tmp_path / "full_complex.json"
    path.write_text(json.dumps([
        {"bgn": endpoint("X", 1, "E2L", "B"), "end": endpoint("H", 35, "SER", "P"),
         "contact": ["hbond"], "type": "atom-atom"},
        {"bgn": endpoint("H", 100, "LYS", "P"), "end": endpoint("H", 108, "ASP", "P"),
         "contact": ["ionic"], "type": "atom-atom"},
        {"bgn": endpoint("H", 100, "LYS", "P"), "end": endpoint("L", 55, "GLU", "P"),
         "contact": ["hbond", "ionic"], "type": "atom-atom"},
    ]))
    pairs = read_arpeggio(path, {"H", "L"})
    assert len(pairs) == 3
    assert sum(row.get("count_ionic", 0) for row in pairs) == 2
    counts = _arpeggio_counts(path, "E2L")
    assert counts["hbond"] == 1
    assert counts["ionic"] == 0
    assert counts["salt_bridge"] == "NA"


def test_candidate_network_summary_deduplicates_shared_mutation_pair(tmp_path):
    shared = {
        "parent_mutated_residue": "H:99", "mutant_mutated_residue": "H:99",
        "parent_partner": "H:100", "mutant_partner": "H:100",
        "formed_interaction_types": "hbond;vdw", "broken_interaction_types": "ionic",
        "parent_total": "-1.0", "mutant_total": "-2.0", "delta_total": "-1.0",
    }
    duplicate = {**shared, "parent_mutated_residue": "H:100", "mutant_mutated_residue": "H:100",
                 "parent_partner": "H:99", "mutant_partner": "H:99"}
    columns, type_rows = _candidate_network_summary("candidate.pdb", tmp_path, [shared, duplicate])
    assert columns["network_added_pair_type_events"] == 2
    assert columns["network_removed_pair_type_events"] == 1
    assert columns["network_changed_unique_pairs"] == 1
    assert columns["network_changed_unique_pairs_energy_delta_REU"] == -1.0
    assert len(type_rows) == 3


def test_hbond_saltbridge_summary_splits_scope_and_deduplicates_sources(tmp_path):
    changes = [
        {
            "mutation": "E99W", "antibody_number": "H107",
            "parent_mutated_residue": "H:99", "mutant_mutated_residue": "H:99",
            "parent_partner": "L:55", "mutant_partner": "L:55", "pair_scope": "aa_to_aa",
            "formed_interaction_types": "hbond;rosetta_hbond", "broken_interaction_types": "ionic",
            "parent_total": "-1.0", "mutant_total": "-2.0", "delta_total": "-1.0",
        },
        {
            "mutation": "E99W", "antibody_number": "H107",
            "parent_mutated_residue": "H:99", "mutant_mutated_residue": "H:99",
            "parent_partner": "X:1", "mutant_partner": "X:1", "pair_scope": "aa_to_ligand",
            "formed_interaction_types": "rosetta_hbond", "broken_interaction_types": "none",
            "parent_total": "-2.0", "mutant_total": "-3.0", "delta_total": "-1.0",
        },
    ]
    columns, rows = _hbond_saltbridge_summary("candidate.pdb", tmp_path, changes)
    assert columns["formed_inter_chain_hbond_pair_count"] == 1
    assert columns["formed_protein_ligand_hbond_pair_count"] == 1
    assert columns["broken_inter_chain_salt_bridge_pair_count"] == 1
    assert columns["formed_inter_chain_hbond_energy_delta_REU"] == -1.0
    assert len(rows) == 3


def test_ligand_hydrogen_validation(tmp_path):
    import pytest
    pytest.importorskip("rdkit", reason="Ligand chemical validation requires the analysis extra")
    pdb = tmp_path / "ligand.pdb"
    params = tmp_path / "ligand.params"
    pdb.write_text(
        "HETATM    1  O1  LIG X   1       0.000   0.000   0.000  1.00 20.00           O  \n"
        "HETATM    2  H1  LIG X   1       0.990   0.000   0.000  1.00 20.00           H  \nEND\n"
    )
    params.write_text(
        "ATOM O1 OH X -0.5\nATOM H1 Hpol X 0.5\nBOND_TYPE O1 H1 1\n"
    )
    template = tmp_path / "custom_template.pdb"
    template.write_text(pdb.read_text())
    result = validate_ligand_against_params(pdb, params, "X", template=template)
    assert result["status"] == "complete_explicit_hydrogens"
    assert result["explicit_hydrogen_count"] == 1
    assert not result["bond_length_outliers"]
