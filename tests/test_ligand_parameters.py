"""Explicit charge-model handling in ligand parameter generation."""

from unittest.mock import Mock

import pytest

from fv_energy_analysis import ligand_parameters


def mol2_text(model="USER_CHARGES", charge="-0.83456"):
    return (
        f"@<TRIPOS>MOLECULE\nwater\n3 2 1 0 0\nSMALL\n{model}\n\n"
        "@<TRIPOS>ATOM\n"
        f"1 O1 0.0000 0.0000 0.0000 O.3 1 LIG {charge}\n"
        "2 H1 0.9600 0.0000 0.0000 H 1 LIG 0.41728\n"
        "3 H2 -0.2400 0.9300 0.0000 H 1 LIG 0.41728\n"
        "@<TRIPOS>BOND\n1 1 2 1\n2 1 3 1\n"
    )


def test_mol2_reader_preserves_declared_charge_values(tmp_path):
    path = tmp_path / "lig.mol2"
    path.write_text(mol2_text())
    atoms, model = ligand_parameters._mol2_atoms(path)
    assert model == "USER_CHARGES"
    assert [float(atom[8]) for atom in atoms] == [-0.83456, 0.41728, 0.41728]


@pytest.mark.parametrize("model", ["NO_CHARGES", "NONE", ""])
def test_mol2_without_a_charge_model_is_rejected(tmp_path, model):
    path = tmp_path / "lig.mol2"
    path.write_text(mol2_text(model=model))
    with pytest.raises(ValueError, match="no declared partial-charge model"):
        ligand_parameters._mol2_atoms(path)


@pytest.mark.parametrize("charge", ["nan", "inf", "-inf", "not-a-number"])
def test_bad_charge_values_are_rejected(tmp_path, charge):
    path = tmp_path / "lig.mol2"
    path.write_text(mol2_text(charge=charge))
    with pytest.raises(ValueError, match="partial charge"):
        ligand_parameters._mol2_atoms(path)


def test_missing_atom_charge_column_is_rejected(tmp_path):
    path = tmp_path / "lig.mol2"
    path.write_text(mol2_text().replace(" LIG -0.83456", " LIG"))
    with pytest.raises(ValueError, match="partial-charge column"):
        ligand_parameters._mol2_atoms(path)


def test_mol2_does_not_require_openbabel(tmp_path, monkeypatch):
    pytest.importorskip("rdkit.Chem")
    path = tmp_path / "lig.mol2"
    path.write_text(mol2_text())
    lookup = Mock(side_effect=AssertionError("Open Babel must not be queried for MOL2 input"))
    monkeypatch.setattr(ligand_parameters.shutil, "which", lookup)
    generator = Mock(side_effect=RuntimeError("reached-generator"))
    monkeypatch.setattr(ligand_parameters, "molfile_to_params_path", generator)
    with pytest.raises(RuntimeError, match="reached-generator"):
        ligand_parameters.prepare_parameters(path, tmp_path / "out")
    lookup.assert_not_called()
    generator.assert_called_once()
    assert (tmp_path / "out" / "ligand.mol2").read_bytes() == path.read_bytes()


def test_sdf_conversion_explicitly_requests_gasteiger(tmp_path, monkeypatch):
    Chem = pytest.importorskip("rdkit.Chem")
    from rdkit.Chem import AllChem

    molecule = Chem.AddHs(Chem.MolFromSmiles("O"))
    AllChem.EmbedMolecule(molecule, randomSeed=11)
    source = tmp_path / "lig.sdf"
    writer = Chem.SDWriter(str(source))
    writer.write(molecule)
    writer.close()
    monkeypatch.setattr(ligand_parameters.shutil, "which", lambda _: "obabel")
    conversion = Mock(side_effect=RuntimeError("conversion-command-captured"))
    monkeypatch.setattr(ligand_parameters.subprocess, "run", conversion)
    with pytest.raises(RuntimeError, match="conversion-command-captured"):
        ligand_parameters.prepare_parameters(source, tmp_path / "out")
    command = conversion.call_args.args[0]
    assert command[-2:] == ["--partialcharge", "gasteiger"]
