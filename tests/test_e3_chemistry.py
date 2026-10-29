"""E3 probe graph/stereochemistry and exact framework preservation tests."""
from pathlib import Path

import numpy as np
import pytest

Chem = pytest.importorskip("rdkit.Chem")
from rdkit.Chem import AllChem, rdMolDescriptors

from fv_energy_analysis.e3_chemistry import (
    ESTRIOL_FORMULA, ESTRIOL_INCHIKEY, _c16_site, _graft_e3, build_e3_probe,
)
from fv_energy_analysis.ligand_chemistry import ESTRADIOL_SMILES, _params_molecule
from fv_energy_analysis.steric_sandwich import read_steric_atoms


@pytest.fixture
def chemical_input(tmp_path):
    molecule = Chem.AddHs(Chem.MolFromSmiles(ESTRADIOL_SMILES))
    assert AllChem.EmbedMolecule(molecule, randomSeed=9123) == 0
    numbers = {6: 0, 8: 0, 1: 0}
    for atom in molecule.GetAtoms():
        number = atom.GetAtomicNum()
        numbers[number] += 1
        name = atom.GetSymbol() + str(numbers[number])
        info = Chem.AtomPDBResidueInfo()
        info.SetName(f"{name:>4}")
        info.SetResidueName("E2L")
        info.SetResidueNumber(1)
        info.SetChainId("X")
        info.SetIsHeteroAtom(True)
        info.SetOccupancy(1)
        atom.SetMonomerInfo(info)
    template = tmp_path / "E2L_0001.pdb"
    template.write_text(Chem.MolToPDBBlock(molecule))
    names = [atom.GetPDBResidueInfo().GetName().strip() for atom in molecule.GetAtoms()]
    lines = [f"ATOM {name} X X 0.0" for name in names]
    for bond in molecule.GetBonds():
        order = 4 if bond.GetIsAromatic() else int(bond.GetBondTypeAsDouble())
        lines.append(f"BOND_TYPE {names[bond.GetBeginAtomIdx()]} {names[bond.GetEndAtomIdx()]} {order}")
    params = tmp_path / "E2L.params"
    params.write_text("\n".join(lines) + "\n")
    source = tmp_path / "complex.pdb"
    protein = "ATOM   1001  CA  ALA H   1      20.000  20.000  20.000  1.00  0.00           C"
    ligand_lines = [line for line in template.read_text().splitlines() if line.startswith("HETATM")]
    source.write_text(protein + "\n" + "\n".join(ligand_lines) + "\nEND\n")
    return source, params, template


def test_only_correct_c16_stereoisomer_is_accepted(chemical_input):
    _, params, template = chemical_input
    molecule, _ = _params_molecule(params, template)
    original = molecule.GetConformer().GetPositions().copy()
    probe, site = _graft_e3(molecule)
    assert Chem.MolToInchiKey(probe) == ESTRIOL_INCHIKEY
    assert rdMolDescriptors.CalcMolFormula(probe) == ESTRIOL_FORMULA
    assert len(site["tested_directions"]) == 2
    assert sum(item["inchikey"] == ESTRIOL_INCHIKEY for item in site["tested_directions"]) == 1
    assert len(Chem.FindMolChiralCenters(probe, includeUnassigned=True)) == 6
    for atom in molecule.GetAtoms():
        if atom.GetAtomicNum() > 1:
            np.testing.assert_array_equal(original[atom.GetIdx()], probe.GetConformer().GetPositions()[atom.GetIdx()])
    np.testing.assert_array_equal(original, molecule.GetConformer().GetPositions())
    new = probe.GetConformer().GetPositions()
    assert np.linalg.norm(new[site["added_o_index"]]-new[site["c16_index"]]) == pytest.approx(1.43)


def test_export_preserves_source_heavy_coordinate_fields_and_protein(chemical_input, tmp_path):
    source, params, template = chemical_input
    before = source.read_bytes()
    output = tmp_path / "probe.pdb"
    result = build_e3_probe(source, params, template, output)
    assert source.read_bytes() == before
    assert result["identity_validated"] is True
    assert result["shared_heavy_atom_count"] == 20
    assert result["energetics_assessed"] is False
    original = read_steric_atoms(source)
    probe = read_steric_atoms(output)
    before_atoms = {a.name: a for a in original if a.chain == "X"}
    ligand = {a.name: a for a in probe if a.chain == "X"}
    assert len(ligand) == 21
    assert {a.residue for a in ligand.values()} == {"E3L"}
    for name, atom in before_atoms.items():
        np.testing.assert_array_equal(atom.xyz, ligand[name].xyz)
    np.testing.assert_allclose(ligand["O3"].xyz, result["added_o_xyz_A"], atol=.0005)
    protein_before = [line for line in before.decode().splitlines() if line.startswith("ATOM")]
    protein_after = [line for line in output.read_text().splitlines() if line.startswith("ATOM")]
    assert protein_before == protein_after
    all_ligand = {line[12:16].strip(): line for line in output.read_text().splitlines()
                  if line.startswith("HETATM") and line[21:22] == "X"}
    assert len(all_ligand) == 45
    assert sum(line[76:78].strip() == "H" for line in all_ligand.values()) == 24
    assert "H25" in all_ligand and result["replaced_h_atom"] not in all_ligand
    assert result["formal_charge"] == 0
    source_ligand = {line[12:16].strip(): line for line in before.decode().splitlines()
                     if line.startswith("HETATM")}
    for name, line in source_ligand.items():
        if name != result["replaced_h_atom"]:
            assert line[30:54] == all_ligand[name][30:54]
    assert result["ligand_hydrogens_exported"] is True
    assert result["ligand_hydrogen_count"] == 24
    assert set(result["added_h_bond_exclusions"]) == {"H25", "O3", result["c16_atom"]}


def test_added_oh_hydrogen_can_rotate_without_moving_other_atoms(chemical_input, tmp_path):
    source, params, template = chemical_input
    original = build_e3_probe(source, params, template)
    center = np.array(original["added_o_xyz_A"])
    axis = center - np.array(original["c16_xyz_A"])
    axis /= np.linalg.norm(axis)
    h_vector = np.array(original["added_h_xyz_A"]) - center
    # Rotate by90degrees around the C-O axis using Rodrigues' formula.
    rotated_h = center + np.cross(axis, h_vector) + axis * np.dot(axis, h_vector)
    updated = build_e3_probe(source, params, template, tmp_path/"rotated.pdb", added_oh_h_xyz_A=rotated_h)
    np.testing.assert_allclose(updated["added_h_xyz_A"], rotated_h, atol=1e-12)
    np.testing.assert_array_equal(updated["added_o_xyz_A"], original["added_o_xyz_A"])
    assert updated["inchikey"] == ESTRIOL_INCHIKEY
    assert updated["formal_charge"] == 0
    with pytest.raises(ValueError, match="0.96 A"):
        build_e3_probe(source, params, template, added_oh_h_xyz_A=center+2*axis)


def test_direction_rotates_and_translates_with_steroid(chemical_input):
    _, params, template = chemical_input
    molecule, _ = _params_molecule(params, template)
    before, first_site = _graft_e3(molecule)
    rotation = np.array([[0, -1, 0], [1, 0, 0], [0, 0, 1]], dtype=float)
    shift = np.array([10, -5, 2])
    modified = Chem.Mol(molecule)
    for i, xyz in enumerate(molecule.GetConformer().GetPositions()):
        modified.GetConformer().SetAtomPosition(i, rotation @ xyz + shift)
    after, next_site = _graft_e3(modified)
    old_o = before.GetConformer().GetPositions()[first_site["added_o_index"]]
    new_o = after.GetConformer().GetPositions()[next_site["added_o_index"]]
    np.testing.assert_allclose(new_o, rotation @ old_o + shift, atol=1e-12)
    assert next_site["replaced_h_atom"] == first_site["replaced_h_atom"]


def test_wrong_source_stereo_is_not_silently_repaired(chemical_input):
    _, params, template = chemical_input
    molecule, _ = _params_molecule(params, template)
    for i, xyz in enumerate(molecule.GetConformer().GetPositions()):
        molecule.GetConformer().SetAtomPosition(i, [-xyz[0], xyz[1], xyz[2]])
    with pytest.raises(ValueError, match="not canonical E2"):
        _graft_e3(molecule)


def test_missing_hydrogens_and_source_overwrite_are_rejected(chemical_input):
    source, params, template = chemical_input
    with pytest.raises(ValueError, match="differ from the source"):
        build_e3_probe(source, params, template, source)
    lines = [line for line in source.read_text().splitlines()
             if not (line.startswith("HETATM") and line[76:78].strip() == "H")]
    source.write_text("\n".join(lines)+"\n")
    with pytest.raises(ValueError, match="complete explicit-H"):
        build_e3_probe(source, params, template)
