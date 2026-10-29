from dataclasses import replace
import math

import numpy as np
import pytest

from fv_energy_analysis.e3_steric import (
    analyze_added_oxygen_contacts, analyze_e3_steric, orient_e3_hydroxyl,
    _potential_acceptor_for_hydroxyl_h,
)
from fv_energy_analysis.steric_sandwich import StericAtom


def atom(name, xyz, *, element="C", residue="ALA", number=1, chain="H"):
    return StericAtom(name, residue, chain, number, "", element, np.array(xyz, dtype=float))


def test_carbon_overlap_warning_counts_unique_residues_and_preserves_coordinates():
    protein = [atom("CB", (2.7, 0, 0)), atom("CA", (0, 2.8, 0)),
               atom("CB", (0, 0, 3.0), number=2)]
    before = [a.xyz.copy() for a in protein]
    result = analyze_added_oxygen_contacts(protein, (0, 0, 0))
    summary = result["summary"]
    assert summary["e3_steric_warning"] is True
    assert summary["e3_added_o_overlap_residue_count"] == 1
    assert summary["e3_added_o_nonpolar_overlap_residue_count"] == 1
    assert summary["e3_added_o_closest_distance_A"] == pytest.approx(2.7)
    assert summary["e3_added_o_deepest_overlap_A"] == pytest.approx(.52)
    assert result["contacts"][1]["is_backbone"] is True
    for a, original in zip(protein, before):
        np.testing.assert_array_equal(a.xyz, original)


def test_missing_polar_geometry_is_unassessed_not_false_or_hbond_exempt():
    result = analyze_added_oxygen_contacts(
        [atom("OD1", (2.5, 0, 0), element="O", residue="ASP")], (0, 0, 0))
    summary = result["summary"]
    assert summary["e3_steric_warning"] is None
    assert summary["e3_steric_assessment"] == "unassessed_polar_geometry"
    assert summary["e3_added_o_polar_close_residue_count"] == 1
    assert summary["e3_added_o_overlap_residue_count"] == 1
    assert result["contacts"][0]["hbond_assessed"] is False
    assert summary["e3_added_o_unassessed_overlap_residue_count"] == 1


def test_threshold_boundary_and_no_warning_outside_it():
    for distance, warning in ((2.82, True), (2.82001, False)):
        result = analyze_added_oxygen_contacts([atom("CB", (distance, 0, 0))], (0, 0, 0))
        assert result["summary"]["e3_steric_warning"] is warning
    assert result["summary"]["e3_steric_assessment"] == "no_large_overlap_in_fixed_pose"


def test_assessed_polar_contact_geometry_controls_exemption():
    oxygen = atom("OD1", (2.5, 0, 0), element="O", residue="ASP")
    for state, warning, assessment in (
        ("supported", False, "hbond_supported_overlaps_only"),
        ("not_supported", True, "steric_overlap_detected"),
        ("unassessed", None, "unassessed_polar_geometry"),
    ):
        result = analyze_added_oxygen_contacts(
            [oxygen], (0, 0, 0), polar_evaluations={oxygen.label: {"state": state}})
        assert result["summary"]["e3_steric_warning"] is warning
        assert result["summary"]["e3_steric_assessment"] == assessment
        assert result["contacts"][0]["hbond_geometry"]["state"] == state


def test_certain_clash_remains_true_when_other_polar_geometry_is_missing():
    result = analyze_added_oxygen_contacts(
        [atom("CB", (2.5, 0, 0)), atom("O", (0, 2.5, 0), element="O")], (0, 0, 0))
    assert result["summary"]["e3_steric_warning"] is True
    assert result["summary"]["e3_added_o_non_hbond_overlap_residue_count"] == 1
    assert result["summary"]["e3_added_o_unassessed_overlap_residue_count"] == 1


def test_no_neighbors_within_output_cutoff_still_reports_true_nearest_distance():
    result = analyze_added_oxygen_contacts([atom("CB", (10, 0, 0))], (0, 0, 0))
    assert not result["contacts"]
    assert result["summary"]["e3_added_o_closest_distance_A"] == 10
    assert result["summary"]["e3_added_o_deepest_overlap_A"] == 0


def test_rigid_rotation_and_translation_leave_metrics_unchanged():
    protein = [atom("CB", (2.6, 0, 0)), atom("O", (0, 3, 0), element="O")]
    result = analyze_added_oxygen_contacts(protein, (0, 0, 0))
    rot = np.array(((0, -1, 0), (1, 0, 0), (0, 0, 1)))
    shift = np.array((123, -7, 11))
    moved = [replace(a, xyz=rot @ a.xyz + shift) for a in protein]
    transformed = analyze_added_oxygen_contacts(moved, shift)
    for key, value in result["summary"].items():
        if isinstance(value, float):
            assert transformed["summary"][key] == pytest.approx(value)
        else:
            assert transformed["summary"][key] == value


def test_invalid_input_fails_explicitly():
    for xyz in ((math.nan, 0, 0), (1, 2)):
        with pytest.raises(ValueError, match="three finite"):
            analyze_added_oxygen_contacts([atom("CB", (3, 0, 0))], xyz)
    with pytest.raises(ValueError, match="No protein"):
        analyze_added_oxygen_contacts([], (0, 0, 0))
    with pytest.raises(ValueError, match="Unsupported"):
        analyze_added_oxygen_contacts([atom("ZN", (3, 0, 0), element="ZN")], (0, 0, 0))
    for setting in (0, -1, math.inf):
        with pytest.raises(ValueError, match="positive and finite"):
            analyze_added_oxygen_contacts([atom("CB", (3, 0, 0))], (0, 0, 0), overlap_warning_A=setting)


def test_public_probe_reader_and_coordinate_override(tmp_path):
    def line(serial, name, xyz, record, residue, chain, element):
        x, y, z = xyz
        return (f"{record:<6}{serial:5d} {name:>4} {residue:>3} {chain}   1    "
                f"{x:8.3f}{y:8.3f}{z:8.3f}{1:6.2f}{0:6.2f}          {element:>2}\n")
    text = line(1, "CB", (2.7, 0, 0), "ATOM", "ALA", "H", "C")
    text += line(2, "O3", (0, 0, 0), "HETATM", "E3L", "X", "O")
    text += "".join(line(i+2, f"C{i}", (20+i, 0, 0), "HETATM", "E3L", "X", "C")
                    for i in range(1, 21))
    path = tmp_path / "probe.pdb"
    path.write_text(text)
    before = path.read_bytes()
    result = analyze_e3_steric(path)
    assert result["summary"]["e3_steric_warning"] is True
    assert path.read_bytes() == before
    assert result["provenance"]["pdb_sha256"]
    source = tmp_path / "original.pdb"
    source.write_text(text.splitlines(keepends=True)[0])
    override = analyze_e3_steric(source, added_o_xyz_A=(0, 0, 0))
    assert override["summary"] == result["summary"]
    with pytest.raises(ValueError, match="twenty-one-heavy"):
        analyze_e3_steric(source)


def full_h_probe(tmp_path, *, blocker=None, name="probe.pdb"):
    def line(serial, atom_name, xyz, record, residue, chain, element):
        x, y, z = xyz
        return (f"{record:<6}{serial:5d} {atom_name:>4} {residue:>3} {chain}   1    "
                f"{x:8.3f}{y:8.3f}{z:8.3f}{1:6.2f}{0:6.2f}          {element:>2}\n")
    theta = math.radians(109.5)
    hydrogen = np.array((1.43-.96*math.cos(theta), .96*math.sin(theta), 0))
    entries = [("CA", (30, 20, 20), "ATOM", "ALA", "H", "C"),
               ("C12", (0, 0, 0), "HETATM", "E3L", "X", "C"),
               ("O3", (1.43, 0, 0), "HETATM", "E3L", "X", "O"),
               ("H25", hydrogen, "HETATM", "E3L", "X", "H")]
    entries += [(f"Q{i}", (20+i, 10, 0), "HETATM", "E3L", "X", "C") for i in range(19)]
    entries += [(f"H{i}", (20+i, 11, 0), "HETATM", "E3L", "X", "H") for i in range(1, 24)]
    if blocker:
        element, position = blocker
        entries += [("H1" if element == "H" else "O", position, "ATOM", "ALA", "H", element)]
    path = tmp_path / name
    path.write_text("".join(line(i, *entry) for i, entry in enumerate(entries, 1)))
    return path, hydrogen


def test_hydroxyl_sampling_only_moves_new_h_and_avoids_h_h_blocker(tmp_path):
    _, h0 = full_h_probe(tmp_path)
    path, _ = full_h_probe(tmp_path, blocker=("H", h0), name="blocked.pdb")
    before = path.read_text().splitlines()
    result = orient_e3_hydroxyl(path)
    after = path.read_text().splitlines()
    assert result["orientation_count"] == 12
    assert len(result["orientations"]) == 12
    assert result["selected_orientation_index"] != 0
    assert result["maximum_nonpolar_overlap_A"] > result["minimum_nonpolar_overlap_A"]
    assert result["minimum_nonpolar_overlap_A"] == pytest.approx(min(
        row["maximum_nonpolar_overlap_A"] for row in result["orientations"]))
    assert result["minimum_protein_nonpolar_overlap_A"] == pytest.approx(min(
        row["maximum_protein_nonpolar_overlap_A"] for row in result["orientations"]))
    for first, second in zip(before, after):
        if first[12:16].strip() != "H25":
            assert first == second
        else:
            assert first[:30] == second[:30] and first[54:] == second[54:]
    xyz = np.array(result["added_h_xyz_A"])
    assert np.linalg.norm(xyz-np.array((1.43, 0, 0))) == pytest.approx(.96, abs=.001)
    assert result["protocol"]["energies_assessed"] is False


def test_polar_hydrogen_acceptor_overlap_is_not_nonpolar_clash(tmp_path):
    _, h0 = full_h_probe(tmp_path)
    path, _ = full_h_probe(tmp_path, blocker=("O", h0+np.array((1.8, 0, 0))), name="polar.pdb")
    result = orient_e3_hydroxyl(path)
    contact = next(r for r in result["orientations"][0]["contacts"] if r["atom"] == "H:1:ALA:O")
    assert contact["potential_hbond_acceptor"] is True
    assert contact["vdw_overlap_A"] > .4
    assert contact["nonpolar_overlap_warning"] is False
    assert contact["short_polar_contact"] is False


def test_standard_donor_nitrogens_not_treated_as_acceptors():
    nd2 = atom("ND2", (0, 0, 0), element="N", residue="ASN")
    assert _potential_acceptor_for_hydroxyl_h(nd2, [nd2]) is False
    histidine = atom("NE2", (0, 0, 0), element="N", residue="HIS")
    assert _potential_acceptor_for_hydroxyl_h(histidine, [histidine]) is True
    hydrogen = atom("HE2", (1.0, 0, 0), element="H", residue="HIS")
    assert _potential_acceptor_for_hydroxyl_h(histidine, [histidine, hydrogen]) is False


def test_hydroxyl_sampler_rejects_missing_h_or_invalid_sampling(tmp_path):
    path, _ = full_h_probe(tmp_path)
    for n in (0, -1, True, 1.5):
        with pytest.raises(ValueError, match="positive integer"):
            orient_e3_hydroxyl(path, n_orientations=n)
    path.write_text("".join(line for line in path.read_text().splitlines(keepends=True)
                            if line[12:16].strip() != "H25"))
    with pytest.raises(ValueError, match="24 hydrogens"):
        orient_e3_hydroxyl(path)


def test_hbond_exemptions_must_coexist_in_one_saved_hydrogen_orientation(tmp_path, monkeypatch):
    path, h0 = full_h_probe(tmp_path)
    def protein_oxygen(serial, number, y):
        return (f"ATOM  {serial:5d}  OD1 ASP H{number:4d}    "
                f"{1.43:8.3f}{y:8.3f}{0:8.3f}{1:6.2f}{0:6.2f}           O\n")
    path.write_text(path.read_text()+protein_oxygen(100, 2, 2.5)+protein_oxygen(101, 3, -2.5))
    def classify(atom, protein, oxygen, added_h_positions_A):
        positive = added_h_positions_A[0][1] > 0
        supports = (atom.number == 2) == positive
        return {"status": "supported" if supports else "not_supported", "trials": []}
    monkeypatch.setattr("fv_energy_analysis.e3_hbond.classify_polar_overlap", classify)
    positions = [
        {"orientation_index": 0, "hydrogen_xyz_A": h0.tolist()},
        {"orientation_index": 1, "hydrogen_xyz_A": [h0[0], -h0[1], h0[2]]},
    ]
    before = path.read_bytes()
    result = analyze_e3_steric(path, hydrogen_positions=positions)
    assert result["summary"]["e3_steric_warning"] is True
    assert result["summary"]["e3_added_o_non_hbond_overlap_residue_count"] == 1
    assert result["summary"]["e3_added_o_hbond_supported_overlap_residue_count"] == 1
    assert result["selected_hydrogen_orientation_index"] == 0
    assert all(row["certain_clash_atom_count"] == 1
               for row in result["hydrogen_orientation_selection"]["per_orientation"])
    assert path.read_bytes() == before


def test_hbond_geometry_selection_does_not_prefer_a_saved_h_steric_blocker(tmp_path, monkeypatch):
    path, h0 = full_h_probe(tmp_path)
    line = (f"ATOM  {100:5d}  OD1 ASP H{2:4d}    "
            f"{1.43:8.3f}{2.5:8.3f}{0:8.3f}{1:6.2f}{0:6.2f}           O\n")
    path.write_text(path.read_text()+line)
    def classify(atom, protein, oxygen, added_h_positions_A):
        return {"status": "supported" if added_h_positions_A[0][1] > 0 else "not_supported", "trials": []}
    monkeypatch.setattr("fv_energy_analysis.e3_hbond.classify_polar_overlap", classify)
    positions = [
        {"orientation_index": 0, "hydrogen_xyz_A": h0.tolist(),
         "nonpolar_warning_atom_count": 1, "short_polar_contact_count": 0},
        {"orientation_index": 1, "hydrogen_xyz_A": [h0[0], -h0[1], h0[2]],
         "nonpolar_warning_atom_count": 0, "short_polar_contact_count": 0},
    ]
    result = analyze_e3_steric(path, hydrogen_positions=positions)
    assert result["selected_hydrogen_orientation_index"] == 1
    assert result["summary"]["e3_steric_warning"] is True
    assert result["hydrogen_orientation_selection"]["eligible_orientation_count"] == 1
