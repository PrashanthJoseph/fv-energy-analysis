"""Polar overlap exemptions require chemistry and explicit H geometry."""
import json
from pathlib import Path

import numpy as np
import pytest

from fv_energy_analysis.e3_hbond import classify_polar_overlap
from fv_energy_analysis.steric_sandwich import StericAtom


def atom(name, xyz, *, element="N", residue="TRP", number=35, chain="L"):
    return StericAtom(name, residue, chain, number, "", element, np.array(xyz, dtype=float))


def test_backbone_donor_valid_linear_geometry_is_supported():
    n, h = atom("N", (0, 0, 0)), atom("H", (1, 0, 0), element="H")
    result = classify_polar_overlap(n, [n, h], [2.8, 0, 0], added_h_positions_A=[[2.8, .96, 0]])
    assert result["status"] == "supported"
    assert result["trials"][0]["angle_DHA_deg"] == pytest.approx(180)
    assert result["trials"][0]["hydrogen_acceptor_distance_A"] == pytest.approx(1.8)


def test_backbone_n_is_not_an_acceptor_and_bad_donor_geometry_is_rejected():
    n, h = atom("N", (0, 0, 0)), atom("H", (1.02, 0, 0), element="H")
    result = classify_polar_overlap(n, [n, h], [.7, 2.41, 0], added_h_positions_A=[[.7, 1.45, 0]])
    assert result["status"] == "not_supported"
    assert result["chemistry"]["acceptor"] is False
    assert len(result["trials"]) == 1
    assert result["trials"][0]["angle_DHA_deg"] < 100


def test_missing_donor_h_is_unassessed_not_false_or_supported():
    n = atom("N", (0, 0, 0))
    result = classify_polar_overlap(n, [n], [2.5, 0, 0], added_h_positions_A=[[3.46, 0, 0]])
    assert result["status"] == "unassessed"
    assert "missing" in result["reason"].lower()


def test_nearest_wrong_named_h_is_not_assigned_to_backbone_n():
    n = atom("N", (0, 0, 0))
    wrong_h = atom("HA", (1, 0, 0), element="H")
    result = classify_polar_overlap(n, [n, wrong_h], [2.5, 0, 0], added_h_positions_A=[[3.46, 0, 0]])
    assert result["status"] == "unassessed"
    assert not result["trials"]


def test_ligand_donor_to_carbonyl_acceptor_uses_one_supplied_orientation():
    acceptor = atom("O", (2.8, 0, 0), element="O", residue="ALA")
    good = classify_polar_overlap(acceptor, [acceptor], [0, 0, 0], added_h_positions_A=[[.96, 0, 0]])
    bad = classify_polar_overlap(acceptor, [acceptor], [0, 0, 0], added_h_positions_A=[[-.96, 0, 0]])
    assert good["status"] == "supported"
    assert bad["status"] == "not_supported"
    assert good["trials"][0]["direction"] == "E3_donor_to_protein_acceptor"


def test_too_short_linear_contact_is_not_exempted():
    n, h = atom("N", (0, 0, 0)), atom("H", (1, 0, 0), element="H")
    result = classify_polar_overlap(n, [n, h], [2.0, 0, 0], added_h_positions_A=[[2.0, .96, 0]])
    assert result["status"] == "not_supported"
    assert result["trials"][0]["angle_DHA_deg"] == 180
    assert result["trials"][0]["criteria"]["hydrogen_acceptor_distance_pass"] is False


def test_sidechain_amide_n_donates_but_cannot_accept():
    n = atom("ND2", (0, 0, 0), residue="ASN")
    hs = [atom("HD21", (1, 0, 0), residue="ASN", element="H"),
          atom("HD22", (0, 1, 0), residue="ASN", element="H")]
    result = classify_polar_overlap(n, [n, *hs], [-2.5, 0, 0], added_h_positions_A=[[-1.54, 0, 0]])
    assert result["status"] == "not_supported"
    assert result["chemistry"]["acceptor"] is False
    assert all(t["direction"] == "protein_donor_to_E3_acceptor" for t in result["trials"])


def test_histidine_typing_uses_explicit_ring_protonation():
    nd = atom("ND1", (2.8, 0, 0), residue="HIS")
    ne = atom("NE2", (3.5, 1.0, 0), residue="HIS")
    he = atom("HE2", (4.5, 1.0, 0), residue="HIS", element="H")
    result = classify_polar_overlap(nd, [nd, ne, he], [0, 0, 0], added_h_positions_A=[[.96, 0, 0]])
    assert result["status"] == "supported"
    assert result["chemistry"]["acceptor"] is True
    assert result["chemistry"]["donor"] is False
    assert classify_polar_overlap(nd, [nd, ne], [0, 0, 0], added_h_positions_A=[[.96, 0, 0]])["status"] == "unassessed"


def test_internal_proline_is_neither_donor_nor_acceptor():
    n = atom("N", (0, 0, 0), residue="PRO")
    c = atom("C", (-1.3, 0, 0), residue="ALA", number=34, element="C")
    result = classify_polar_overlap(n, [n, c], [2.5, 0, 0], added_h_positions_A=[[1.54, 0, 0]])
    assert result["status"] == "not_supported"
    assert not result["chemistry"]["donor"] and not result["chemistry"]["acceptor"]


def test_unknown_atom_typing_is_explicitly_unassessed():
    n = atom("NZ", (0, 0, 0), residue="UNK")
    result = classify_polar_overlap(n, [n], [2.5, 0, 0], added_h_positions_A=[[1.54, 0, 0]])
    assert result["status"] == "unassessed"


def test_criteria_are_strict_for_plip_gates_and_inclusive_short_contact_guard():
    n, h = atom("N", (0, 0, 0)), atom("H", (1, 0, 0), element="H")
    result = classify_polar_overlap(n, [n, h], [4.1, 0, 0], added_h_positions_A=[[4.1, .96, 0]])
    assert result["status"] == "not_supported"
    assert result["thresholds"]["donor_acceptor_distance_A"] == {"operator": "<", "value": 4.1}
    result = classify_polar_overlap(n, [n, h], [2.2, 0, 0], added_h_positions_A=[[2.2, .96, 0]])
    assert result["status"] == "supported"


@pytest.mark.parametrize("xyz", [[float("nan"), 0, 0], [1, 2]])
def test_invalid_coordinates_fail(xyz):
    n = atom("N", (0, 0, 0))
    with pytest.raises(ValueError):
        classify_polar_overlap(n, [n], xyz)


def test_actual_rank1_backbone_contact_is_not_hbond_supported():
    folder = Path(__file__).resolve().parents[1] / "results/full_2436_e3_steric_20261027/cases/set_A_3mut_type4_contiguous_0_anchored_0_E1Y_R2T_P5D_estradiol_model_0_relaxed"
    path = folder / "e3_probe.pdb"
    if not path.is_file():
        pytest.skip("Local audited rank1 E3 probe absent")
    atoms = []
    for line in path.read_text().splitlines():
        if line.startswith(("ATOM  ", "HETATM")):
            atoms.append(StericAtom(line[12:16].strip(), line[17:20].strip(), line[21], int(line[22:26]),
                                    line[26].strip(), line[76:78].strip(),
                                    np.array([float(line[i:i+8]) for i in (30, 38, 46)]), record=line[:6].strip()))
    n = next(a for a in atoms if a.chain == "L" and a.number == 35 and a.name == "N")
    oxygen = next(a for a in atoms if a.residue == "E3L" and a.name == "O3")
    hydrogen = next(a for a in atoms if a.residue == "E3L" and a.name == "H25")
    result = classify_polar_overlap(n, [a for a in atoms if a.record == "ATOM"], oxygen.xyz,
                                    added_h_positions_A=[hydrogen.xyz])
    assert result["status"] == "not_supported"
    trial = result["trials"][0]
    assert trial["angle_DHA_deg"] == pytest.approx(83.1665752, abs=1e-6)
    assert trial["donor_acceptor_distance_A"] == pytest.approx(2.5062111, abs=1e-6)
    assert trial["hydrogen_acceptor_distance_A"] == pytest.approx(2.4138362, abs=1e-6)
    json.dumps(result, allow_nan=False)
