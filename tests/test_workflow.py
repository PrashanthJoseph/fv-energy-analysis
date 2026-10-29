"""Parent-cache execution and the two supported ligand-input modes."""

import json
from pathlib import Path
from unittest.mock import Mock

import pytest

from fv_energy_analysis import workflow
from fv_energy_analysis.cli import build_parser as cli_parser, score_pyrosetta


@pytest.mark.parametrize("same_complex", [True, False])
def test_workflow_runs_each_distinct_complex_once(tmp_path, monkeypatch, same_complex):
    parent = tmp_path / "parent.pdb"
    parent.write_text("END\n")
    mutant = parent if same_complex else tmp_path / "mutant" / "parent.pdb"
    mutant.parent.mkdir(exist_ok=True)
    mutant.write_text("END\n")
    params, template = tmp_path / "lig.params", tmp_path / "lig.pdb"
    params.write_text("NAME LIG\n")
    template.write_text("END\n")
    output = tmp_path / "results"
    commands = []

    def run(runner, command, **kwargs):
        command = [str(value) for value in command]
        commands.append(command)
        runner.commands.append(command)

        def value(flag):
            return command[command.index(flag) + 1]

        if "--generated-ligand-pdb" in command:
            destination = Path(value("--output")) / f"{Path(value('--complex')).stem}_rosetta.pdb"
            destination.write_text("END\n")
        elif "score-pyrosetta" in command:
            for flag in ("--pair-output", "--interface-output", "--hbond-output"):
                Path(value(flag)).write_text("saved result\n")
        elif command[0] == "plip":
            directory, name = Path(value("-o")), value("--name")
            (directory / f"{name}.xml").write_text("<report/>\n")
            (directory / f"{name}.txt").write_text("saved result\n")
        elif Path(command[1]).name == "pdb_to_arpeggio_cif.py":
            Path(command[-1]).write_text("data_test\n")

    def arpeggio(runner, executable, directory, cif):
        path = directory / f"{cif.stem}.json"
        path.write_text("[]\n")
        return path

    monkeypatch.setattr(workflow.Runner, "run", run)
    monkeypatch.setattr(workflow.shutil, "which", lambda name: name)
    detectors = Mock(side_effect=arpeggio)
    monkeypatch.setattr(workflow, "run_arpeggio", detectors)
    atoms = Mock(return_value={"valid": True})
    monkeypatch.setattr(workflow, "validate_ligand_against_params", atoms)
    plip = Mock(return_value={"interaction_counts": {}})
    monkeypatch.setattr(workflow, "summarize_plip_xml", plip)
    architecture = Mock(return_value={})
    monkeypatch.setattr(workflow, "analyze_estrogen_architecture", architecture)
    monkeypatch.setattr(workflow, "_write_automatic_summary", lambda *args: None)
    chemistry = Mock(return_value="CCO")
    monkeypatch.setattr("fv_energy_analysis.ligand_chemistry.reference_smiles_from_source", chemistry)

    args = workflow.build_parser().parse_args([
        "--parent", str(parent), "--mutant", str(mutant), "--output", str(output),
        "--params", str(params), "--ligand-template-pdb", str(template),
        "--arpeggio-scope", "full",
    ])
    workflow.run_workflow(args)

    count = 1 if same_complex else 2
    assert len([command for command in commands if "--generated-ligand-pdb" in command]) == count
    scores = [command for command in commands if "score-pyrosetta" in command]
    assert len(scores) == count
    assert all("--hbond-output" in command for command in scores)
    preparations = [command for command in commands if "--generated-ligand-pdb" in command]
    assert all(command[command.index("--params") + 1] == str(params.resolve()) for command in preparations)
    assert all(call.kwargs["template"] == template.resolve() for call in atoms.call_args_list)
    if not same_complex:
        assert scores[0][scores[0].index("--pdb") + 1] != scores[1][scores[1].index("--pdb") + 1]
    assert not any("extract-hbonds" in command for command in commands)
    assert len([command for command in commands if command[0] == "plip"]) == count
    assert detectors.call_count == atoms.call_count == plip.call_count == architecture.call_count == count
    chemistry.assert_called_once_with(None, params.resolve(), template.resolve())
    manifest = json.loads((output / "workflow_manifest.json").read_text())
    assert manifest["status"] == "complete"
    assert manifest["parent_artifacts_shared"] is same_complex
    if same_complex:
        for suffix in ("pair_energies.sc", "interface.sc", "hbonds.csv"):
            assert (output / f"parent_{suffix}").read_bytes() == (output / f"mutant_{suffix}").read_bytes()
        assert (output / "arpeggio_mutant" / "mutant_arpeggio.json").is_file()
        assert (output / "plip_mutant" / "mutant_plip.xml").is_file()


def test_workflow_preserves_existing_output(tmp_path):
    saved = tmp_path / "workflow.log"
    saved.write_text("existing results\n")
    args = workflow.build_parser().parse_args([
        "--parent", "parent.pdb", "--mutant", "mutant.pdb", "--output", str(tmp_path),
    ])
    with pytest.raises(FileExistsError, match="not empty"):
        workflow.run_workflow(args)
    assert saved.read_text() == "existing results\n"


def test_score_cli_extracts_hbonds_without_reloading_the_pose(monkeypatch):
    combined = Mock()
    separate = Mock()
    monkeypatch.setattr("fv_energy_analysis.cli.score_and_extract_hbonds", combined)
    monkeypatch.setattr("fv_energy_analysis.cli.score_with_pyrosetta", separate)
    args = cli_parser().parse_args([
        "score-pyrosetta", "--pdb", "p.pdb", "--pair-output", "pairs.sc",
        "--interface-output", "interface.sc", "--hbond-output", "hbonds.csv",
        "--ligand-chain", "X", "--params", "lig.params",
    ])
    assert score_pyrosetta(args) == 0
    combined.assert_called_once_with(
        "p.pdb", "pairs.sc", "interface.sc", "hbonds.csv", "X",
        params=["lig.params"], weights="ref2015",
    )
    separate.assert_not_called()


def test_existing_params_and_template_need_no_sdf_or_provenance_file(tmp_path):
    Chem = pytest.importorskip("rdkit.Chem")
    from rdkit.Chem import AllChem
    from fv_energy_analysis.ligand_chemistry import ESTRADIOL_SMILES, reference_smiles_from_source

    molecule = Chem.AddHs(Chem.MolFromSmiles(ESTRADIOL_SMILES))
    assert AllChem.EmbedMolecule(molecule, randomSeed=71) == 0
    for atom in molecule.GetAtoms():
        info = Chem.AtomPDBResidueInfo()
        info.SetName(f"{atom.GetSymbol()}{atom.GetIdx() + 1}".ljust(4))
        info.SetResidueName("E2L")
        info.SetResidueNumber(1)
        info.SetChainId("X")
        info.SetIsHeteroAtom(True)
        atom.SetMonomerInfo(info)
    template, params = tmp_path / "E2L.pdb", tmp_path / "E2L.params"
    # No CONECT records: bond orders must come from params, not PDB inference.
    template.write_text("\n".join(
        line for line in Chem.MolToPDBBlock(molecule).splitlines()
        if not line.startswith("CONECT")
    ) + "\n")
    names = [atom.GetPDBResidueInfo().GetName().strip() for atom in molecule.GetAtoms()]
    params.write_text("\n".join(f"ATOM {name} X X 0.0" for name in names) + "\n" + "\n".join(
        f"BOND_TYPE {names[bond.GetBeginAtomIdx()]} {names[bond.GetEndAtomIdx()]} "
        f"{4 if bond.GetIsAromatic() else int(bond.GetBondTypeAsDouble())}"
        for bond in molecule.GetBonds()
    ) + "\n")
    recovered = reference_smiles_from_source(params=params, template=template)
    assert Chem.MolToInchiKey(Chem.MolFromSmiles(recovered)) == Chem.MolToInchiKey(molecule)
    assert not (tmp_path / "ligand_definition.json").exists()


@pytest.mark.parametrize("filename,relocated", [
    ("candidate_rosetta.pdb", False),
    ("candidate_input_rosetta.pdb", False),
    ("different_recorded_name.pdb", False),
    ("candidate_rosetta.pdb", True),
])
def test_rebase_uses_recorded_mutant_not_a_filename_pattern(tmp_path, monkeypatch, filename, relocated):
    from fv_energy_analysis import fast_workflow

    source, cache, output = [tmp_path / name for name in ("saved", "parent_cache", "rebased")]
    for directory in (source / "prepared", source / "final", source / "arpeggio_mutant", cache):
        directory.mkdir(parents=True)
    (source / "workflow_manifest.json").write_text("{}")
    saved_mutant = source / "prepared" / filename
    saved_mutant.write_text("saved mutant coordinates\n")
    recorded = tmp_path / "previous_location" / filename if relocated else saved_mutant
    (source / "final" / "summary.json").write_text(json.dumps({"mutant": str(recorded)}))
    # This tempting filename must not influence the selection.
    (source / "prepared" / "candidate_decoy_input_rosetta.pdb").write_text("wrong coordinates\n")
    for name in ("mutant_pair_energies.sc", "mutant_interface.sc", "mutant_hbonds.csv"):
        (source / name).write_text("saved scores\n")
    (source / "arpeggio_mutant" / "mutant_arpeggio.json").write_text("[]\n")
    for name in ("plip_ligand_interactions.json", "ligand_atom_validation.json", "estrogen_binding_architecture.json"):
        (source / name).write_text(json.dumps({"mutant": {}}))
    cached = {}
    for name in ("prepared", "pairs", "interface", "hbonds", "arpeggio", "validation", "architecture", "plip"):
        path = cache / ("parent_rosetta.pdb" if name == "prepared" else name)
        path.write_text(json.dumps({"parent": {}}))
        cached[name] = path
    monkeypatch.setattr(fast_workflow, "_cache_paths", lambda _: cached)
    monkeypatch.setattr(fast_workflow, "validate_scoring_reuse", lambda *args: None)
    monkeypatch.setattr(fast_workflow, "validate_arpeggio_coverage", lambda *args: None)
    monkeypatch.setattr(fast_workflow, "_write_automatic_summary", lambda *args: None)
    comparison, scoring = Mock(), Mock()
    monkeypatch.setattr(fast_workflow, "compare", comparison)
    monkeypatch.setattr(fast_workflow, "score_and_extract_hbonds", scoring)
    config = dict(
        candidate=str(tmp_path / "candidate.pdb"), source_result=str(source),
        output=str(output), parent_cache=str(cache), parent="parent.pdb", params="lig.params",
        weights="ref2015", fv_chains="H,L", rosetta_ligand_chain="X",
        scheme="imgt", cdr_definition="imgt", fit_region="fv", rmsd_atoms="backbone",
        rotamer_rmsd_cutoff=3.0, network_energy_presence_cutoff=0.1, network_energy_change_cutoff=0.25,
    )
    # Reading params for the output manifest does not trigger any scoring.
    config["params"] = str(tmp_path / "lig.params")
    Path(config["params"]).write_text("saved parameters\n")
    result = fast_workflow.run_rebased_candidate(config)

    assert result["status"] == "complete", result.get("error")
    comparison.assert_called_once()
    chosen = Path(comparison.call_args.args[0].mutant)
    assert chosen.name == filename
    assert chosen.read_bytes() == saved_mutant.read_bytes()
    assert (output / "mutant_pair_energies.sc").read_bytes() == (source / "mutant_pair_energies.sc").read_bytes()
    scoring.assert_not_called()


@pytest.mark.parametrize("relocated", [False, True])
def test_parent_cache_uses_recorded_parent_not_first_pdb(tmp_path, relocated):
    from fv_energy_analysis.fast_workflow import _cache_paths

    (tmp_path / "prepared").mkdir()
    (tmp_path / "final").mkdir()
    parent = tmp_path / "prepared" / "z_actual_parent_rosetta.pdb"
    parent.write_text("actual parent\n")
    (tmp_path / "prepared" / "a_wrong_parent_rosetta.pdb").write_text("wrong parent\n")
    recorded = tmp_path / "old_location" / parent.name if relocated else parent
    (tmp_path / "final" / "summary.json").write_text(json.dumps({"parent": str(recorded)}))
    (tmp_path / "workflow_manifest.json").write_text(json.dumps({"parent_analysis_input": str(parent)}))

    cached = _cache_paths(tmp_path)
    assert cached["prepared"] == parent
    assert cached["prepared"].read_text() == "actual parent\n"


def test_missing_recorded_parent_does_not_fall_back_to_another_pdb(tmp_path):
    from fv_energy_analysis.fast_workflow import _cache_paths

    (tmp_path / "prepared").mkdir()
    (tmp_path / "final").mkdir()
    (tmp_path / "prepared" / "unrelated_rosetta.pdb").write_text("wrong parent\n")
    missing = tmp_path / "prepared" / "missing_parent.pdb"
    (tmp_path / "final" / "summary.json").write_text(json.dumps({"parent": str(missing)}))
    (tmp_path / "workflow_manifest.json").write_text("{}\n")

    with pytest.raises(FileNotFoundError, match="Recorded prepared parent PDB is missing"):
        _cache_paths(tmp_path)
