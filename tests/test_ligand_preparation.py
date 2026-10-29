"""Exact ligand preparation without changes to coordinates or input chemistry."""

import importlib.util
import json
import sys

import pytest

Chem = pytest.importorskip("rdkit.Chem")
from rdkit.Chem import AllChem

from fv_energy_analysis.resources import resource_path


@pytest.fixture
def preparation():
    spec = importlib.util.spec_from_file_location("ligand_preparation", resource_path("scripts/prepare_ligand_complexes.py"))
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def ligand_files(tmp_path, smiles="C[C@H](O)c1ccc(O)cc1", reorder=False, rename=False):
    molecule = Chem.AddHs(Chem.MolFromSmiles(smiles))
    assert AllChem.EmbedMolecule(molecule, randomSeed=37) == 0
    for index, atom in enumerate(molecule.GetAtoms(), 1):
        info = Chem.AtomPDBResidueInfo()
        info.SetName(f"{atom.GetSymbol()}{index}".ljust(4))
        info.SetResidueName("E2L")
        info.SetResidueNumber(1)
        info.SetChainId("X")
        info.SetIsHeteroAtom(True)
        atom.SetMonomerInfo(info)
    names = [atom.GetPDBResidueInfo().GetName().strip() for atom in molecule.GetAtoms()]
    lines = [f"ATOM {name} X X 0.0" for name in names]
    lines += [
        f"BOND_TYPE {names[bond.GetBeginAtomIdx()]} {names[bond.GetEndAtomIdx()]} "
        f"{4 if bond.GetIsAromatic() else int(bond.GetBondTypeAsDouble())}"
        for bond in molecule.GetBonds()
    ]
    lines += [f"CHARGE {names[atom.GetIdx()]} FORMAL {atom.GetFormalCharge()}"
              for atom in molecule.GetAtoms() if atom.GetFormalCharge()]
    params = tmp_path / "E2L.params"
    params.write_text("\n".join(lines) + "\n")
    template = Chem.Mol(molecule)
    conformer = template.GetConformer()
    for index, xyz in enumerate(conformer.GetPositions()):
        conformer.SetAtomPosition(index, xyz + [50, 50, 50])
    template_pdb = tmp_path / "E2L_0001.pdb"
    # Deliberately omit template CONECT: bond orders must come from params.
    template_pdb.write_text("\n".join(line for line in Chem.MolToPDBBlock(template).splitlines()
                                     if not line.startswith("CONECT")) + "\n")
    source = Chem.Mol(molecule)
    for index, atom in enumerate(source.GetAtoms(), 1):
        atom.GetPDBResidueInfo().SetChainId("d")
        if rename:
            atom.GetPDBResidueInfo().SetName(f"Q{index}".ljust(4))
    if reorder:
        source = Chem.RenumberAtoms(source, list(reversed(range(source.GetNumAtoms()))))
    pdb = tmp_path / "bound.pdb"
    protein = "ATOM   9999  CA  ALA H   1      20.000  20.000  20.000  1.00 20.00           C  \n"
    pdb.write_text(protein + Chem.MolToPDBBlock(source))
    return pdb, template_pdb, params, molecule, protein


@pytest.mark.parametrize("reorder,rename", [(False, False), (True, False), (True, True)])
def test_exact_mapping_preserves_all_bound_coordinates(tmp_path, preparation, reorder, rename):
    source, template, params, molecule, protein = ligand_files(tmp_path, reorder=reorder, rename=rename)
    destination = tmp_path / "prepared.pdb"
    before_params, before_template = params.read_bytes(), template.read_bytes()
    preparation.rename_ligand_atoms(source, destination, "d", "X", template)
    before = [line for line in source.read_text().splitlines() if line.startswith("HETATM")]
    after = [line for line in destination.read_text().splitlines() if line.startswith("HETATM")]
    assert [line[30:54] for line in after] == [line[30:54] for line in before]
    assert all(line[21] == "X" for line in after)
    assert destination.read_text().startswith(protein)
    assert params.read_bytes() == before_params
    assert template.read_bytes() == before_template
    output = Chem.MolFromPDBFile(str(destination), removeHs=False)
    ligand = [part for part in Chem.GetMolFrags(output, asMols=True) if part.GetNumAtoms() > 1][0]
    assert sum(atom.GetIsAromatic() for atom in ligand.GetAtoms()) == 6
    assert Chem.MolToInchiKey(ligand) == Chem.MolToInchiKey(molecule)


def test_matching_names_cannot_bypass_stereochemistry(tmp_path, preparation):
    source, template, params, molecule, _ = ligand_files(tmp_path)
    lines = source.read_text().splitlines(keepends=True)
    for index, line in enumerate(lines):
        if line.startswith("HETATM"):
            lines[index] = line[:30] + f"{-float(line[30:38]):8.3f}" + line[38:]
    source.write_text("".join(lines))
    destination = tmp_path / "must_not_be_written.pdb"
    with pytest.raises(ValueError, match="stereochemistry"):
        preparation.rename_ligand_atoms(source, destination, "d", "X", template, params)
    assert not destination.exists()


def test_params_formal_charges_are_used(tmp_path, preparation):
    source, template, params, molecule, _ = ligand_files(tmp_path, smiles="C[NH3+]")
    loaded = preparation.ligand_template(template, params)
    assert Chem.GetFormalCharge(loaded) == 1
    assert Chem.MolToInchiKey(loaded) == Chem.MolToInchiKey(molecule)


def test_missing_explicit_hydrogens_are_not_silently_dropped(tmp_path, preparation):
    source, template, params, _, _ = ligand_files(tmp_path)
    lines = source.read_text().splitlines(keepends=True)
    # Strip all hydrogens and CONECT records to make the PDB consistently heavy-only.
    source.write_text("".join(line for line in lines
                              if not line.startswith("CONECT")
                              and not (line.startswith("HETATM") and line[76:78].strip() == "H")))
    with pytest.raises(ValueError):
        preparation.rename_ligand_atoms(source, tmp_path / "out.pdb", "d", "X", template, params)


def test_legacy_smiles_and_sidecar_do_not_override_params_template(tmp_path, preparation, monkeypatch):
    source, template, params, molecule, _ = ligand_files(tmp_path)
    (tmp_path / "ligand_definition.json").write_text(json.dumps({"reference_smiles": "CC"}))
    output = tmp_path / "out"
    monkeypatch.setattr(sys, "argv", [
        "prepare_ligand_complexes.py", "--complex", str(source),
        "--generated-ligand-pdb", str(template), "--params", str(params),
        "--ligand-chain", "d", "--output", str(output),
        "--reference", str(source), "--ligand-smiles", "CC",
    ])
    preparation.main()
    assert (output / "bound_rosetta.pdb").is_file()


def test_explicit_params_can_live_elsewhere(tmp_path, preparation):
    source, template, params, _, _ = ligand_files(tmp_path)
    alternate = tmp_path / "custom_parameters.params"
    params.rename(alternate)
    preparation.rename_ligand_atoms(source, tmp_path / "out.pdb", "d", "X", template, alternate)


def test_chemistry_resolver_prioritizes_params_template_over_source_and_sidecar(tmp_path):
    from fv_energy_analysis.ligand_chemistry import reference_smiles_from_source

    _, template, params, molecule, _ = ligand_files(tmp_path)
    (tmp_path / "ligand_definition.json").write_text(json.dumps({"reference_smiles": "CC"}))
    before = params.read_bytes(), template.read_bytes()
    smiles = reference_smiles_from_source(tmp_path / "ignored.sdf", params, template)
    assert Chem.MolToInchiKey(Chem.MolFromSmiles(smiles)) == Chem.MolToInchiKey(molecule)
    assert (params.read_bytes(), template.read_bytes()) == before


@pytest.mark.parametrize("smiles,charge", [("C[NH3+]", 1), ("CC(=O)[O-]", -1)])
def test_chemistry_reader_uses_formal_and_partial_charges_separately(tmp_path, smiles, charge):
    from fv_energy_analysis.ligand_chemistry import inspect_params_chemistry, reference_smiles_from_source

    _, template, params, molecule, _ = ligand_files(tmp_path, smiles=smiles)
    formal = {atom.GetPDBResidueInfo().GetName().strip(): atom.GetFormalCharge()
              for atom in molecule.GetAtoms()}
    lines = []
    for line in params.read_text().splitlines():
        fields = line.split()
        if fields[0] == "ATOM":
            fields[-1] = str(formal[fields[1]])
        lines.append(" ".join(fields))
    params.write_text("\n".join(lines) + "\n")
    report = inspect_params_chemistry(params, smiles, template)
    assert report["valid"], report["errors"]
    assert report["formal_charge"] == charge
    assert report["partial_charge_sum"] == charge
    recovered = reference_smiles_from_source(params=params)  # Standard companion template.
    assert Chem.MolToInchiKey(Chem.MolFromSmiles(recovered)) == Chem.MolToInchiKey(molecule)


def test_chemistry_reader_takes_elements_from_template_not_atom_names(tmp_path):
    from fv_energy_analysis.ligand_chemistry import inspect_params_chemistry

    _, template, params, molecule, _ = ligand_files(tmp_path, smiles="CCl")
    mapping = {atom.GetPDBResidueInfo().GetName().strip(): f"Q{atom.GetIdx() + 1}"
               for atom in molecule.GetAtoms()}
    params.write_text("\n".join(" ".join(mapping.get(value, value) for value in line.split())
                                for line in params.read_text().splitlines()) + "\n")
    template.write_text("\n".join(
        line[:12] + mapping[line[12:16].strip()].ljust(4) + line[16:]
        if line.startswith("HETATM") else line
        for line in template.read_text().splitlines()
    ) + "\n")
    report = inspect_params_chemistry(params, "CCl", template)
    assert report["valid"], report["errors"]
    assert report["atom_elements"] == {
        f"Q{atom.GetIdx() + 1}": atom.GetSymbol() for atom in molecule.GetAtoms()
    }


def test_template_is_not_silently_used_without_params(tmp_path):
    from fv_energy_analysis.ligand_chemistry import reference_smiles_from_source

    with pytest.raises(ValueError, match="accompanied by params"):
        reference_smiles_from_source(template=tmp_path / "lig.pdb")


def test_unused_charge_redistribution_function_is_removed():
    from fv_energy_analysis import ligand_chemistry

    assert not hasattr(ligand_chemistry, "correct_charge_rounding")
