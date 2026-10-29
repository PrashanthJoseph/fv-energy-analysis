"""Optional tools are truly omitted and their absence is explicit provenance."""
import json
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest

from fv_energy_analysis import cli, workflow
from fv_energy_analysis.analysis_stages import analysis_status, validate_analysis_reuse, validate_comparison_reuse
from fv_energy_analysis.batch_workflow import build_parser, _parent_cache_key
from fv_energy_analysis.fast_workflow import _cache_paths
from fv_energy_analysis.geometry import Transform
from fv_energy_analysis.models import Residue, ResidueKey


def test_individual_workflow_can_omit_both_optional_detectors(tmp_path, monkeypatch):
    source, params, template = [tmp_path / name for name in ("parent.pdb", "lig.params", "lig.pdb")]
    source.write_text("END\n")
    params.write_text("NAME LIG\n")
    template.write_text("END\n")
    output = tmp_path / "result"
    commands = []

    def run(runner, command, **kwargs):
        command = [str(x) for x in command]
        commands.append(command)
        runner.commands.append(command)
        value = lambda flag: command[command.index(flag) + 1]
        if "--generated-ligand-pdb" in command:
            (Path(value("--output")) / f"{Path(value('--complex')).stem}_rosetta.pdb").write_text("END\n")
        elif "score-pyrosetta" in command:
            for flag in ("--pair-output", "--interface-output", "--hbond-output"):
                Path(value(flag)).write_text("saved\n")
        elif "compare" in command:
            final = Path(value("--output"))
            final.mkdir()
            (final / "summary.json").write_text("{}")

    def forbidden(*args, **kwargs):
        pytest.fail("Disabled detector was looked up or executed")

    monkeypatch.setattr(workflow.Runner, "run", run)
    monkeypatch.setattr(workflow.shutil, "which", forbidden)
    monkeypatch.setattr(workflow, "run_arpeggio", forbidden)
    monkeypatch.setattr(workflow, "summarize_plip_xml", forbidden)
    monkeypatch.setattr(workflow, "validate_ligand_against_params", lambda *a, **k: {"valid": True})
    monkeypatch.setattr(workflow, "analyze_estrogen_architecture", lambda *a, **k: {})
    monkeypatch.setattr(workflow, "_write_automatic_summary", lambda *a: None)
    monkeypatch.setattr("fv_energy_analysis.ligand_chemistry.reference_smiles_from_source", lambda *a: "CCO")
    flags = ["--skip-arpeggio", "--skip-plip", "--skip-network-analysis", "--skip-rotamer-analysis", "--skip-e2-architecture"]
    args = workflow.build_parser().parse_args([
        "--parent", str(source), "--mutant", str(source), "--output", str(output),
        "--params", str(params), "--ligand-template-pdb", str(template), *flags,
    ])
    workflow.run_workflow(args)
    saved = json.loads((output / "workflow_manifest.json").read_text())
    assert saved["status"] == "complete"
    assert set(saved["analysis_status"].values()) == {"not_run"}
    assert saved["arpeggio_selection"] == {"scope": "not_run"}
    assert not (output / "plip_ligand_interactions.json").exists()
    assert not list(output.glob("arpeggio_*/*.json"))
    comparison = next(command for command in commands if "compare" in command)
    assert "--parent-arpeggio" not in comparison
    assert all(flag in comparison for flag in flags)
    assert len([command for command in commands if "score-pyrosetta" in command]) == 1
    assert json.loads((output / "final/summary.json").read_text())["analysis_status"] == saved["analysis_status"]


def test_disabled_detectors_are_not_required_parent_cache_artifacts(tmp_path):
    prepared = tmp_path / "parent.pdb"
    prepared.write_text("END\n")
    (tmp_path / "final").mkdir()
    (tmp_path / "final/summary.json").write_text(json.dumps({"parent": str(prepared)}))
    (tmp_path / "workflow_manifest.json").write_text(json.dumps({
        "parent_analysis_input": str(prepared),
        "analysis_status": {"arpeggio": "not_run", "plip": "not_run"},
    }))
    paths = _cache_paths(tmp_path)
    assert "arpeggio" not in paths and "plip" not in paths
    assert {"pairs", "interface", "hbonds", "validation", "architecture"} <= paths.keys()


def test_reuse_requires_matching_analysis_selection(tmp_path):
    (tmp_path / "workflow_manifest.json").write_text(json.dumps({
        "analysis_status": {"arpeggio": "not_run", "plip": "complete"},
    }))
    validate_analysis_reuse(tmp_path, {"skip_arpeggio": True})
    with pytest.raises(ValueError, match="selection differs"):
        validate_analysis_reuse(tmp_path, {})


def test_detector_selection_changes_parent_cache_key(tmp_path):
    parent = tmp_path / "parent.pdb"
    parent.write_text("END\n")
    args = build_parser().parse_args(["--parent", str(parent), "--input-folder", str(tmp_path),
                                     "--output", str(tmp_path / "results")])
    ordinary = _parent_cache_key(args, parent)
    args.skip_arpeggio = True
    assert _parent_cache_key(args, parent) != ordinary


@pytest.mark.parametrize("setting,value", [("fit_region", "framework"), ("rmsd_atoms", "ca"),
                                           ("rotamer_rmsd_cutoff", 2.0), ("hcdr3_tip_sequence", "GSS"),
                                           ("network_energy_change_cutoff", 0.5)])
def test_changed_comparison_settings_cannot_reuse_old_derived_metrics(tmp_path, setting, value):
    (tmp_path / "final").mkdir()
    (tmp_path / "final/summary.json").write_text("{}")
    validate_comparison_reuse(tmp_path, {})
    with pytest.raises(ValueError, match="comparison setting differs"):
        validate_comparison_reuse(tmp_path, {setting: value})


def test_changed_same_path_mutant_content_cannot_reuse_old_result(tmp_path):
    from hashlib import sha256
    source = tmp_path / "mutant.pdb"
    source.write_text("original coordinates")
    (tmp_path / "workflow_manifest.json").write_text(json.dumps({
        "mutant": str(source), "mutant_input_sha256": sha256(source.read_bytes()).hexdigest(),
    }))
    validate_comparison_reuse(tmp_path, {}, source)
    source.write_text("changed coordinates")
    with pytest.raises(ValueError, match="input content differs"):
        validate_comparison_reuse(tmp_path, {}, source)


def test_explicitly_skipped_rotamers_do_not_calculate_torsions(monkeypatch):
    key = ResidueKey("H", 99)
    parent, mutant = Residue(key, "ALA"), Residue(key, "SER")
    annotation = SimpleNamespace(region_by_key={key: "CDR3"}, number_by_key={key: "H105"}, chain_type="H")
    monkeypatch.setattr(cli, "compare_rotamers", lambda *a: pytest.fail("Rotamers were calculated"))
    rows = cli._mutation_rows({"H": ([parent], [mutant])}, {"H": annotation}, {"H": annotation},
                              {"H_CDR3": 1, "H_CDR3_ANCHORS": 1}, 3, perform_rotamers=False)
    assert rows[0]["rotamer_analysis_status"] == "not_run"
    assert all(rows[0][field] == "" for field in cli.ROTAMER_FIELDS)


def test_explicitly_skipped_network_does_not_run_network_calculation(tmp_path, monkeypatch):
    args = cli.build_parser().parse_args([
        "compare", "--parent", "parent.pdb", "--mutant", "mutant.pdb", "--fv-chains", "H,L",
        "--ligand-chains", "X", "--parent-rosetta", "p.sc", "--mutant-rosetta", "m.sc",
        "--skip-network-analysis", "--output", str(tmp_path),
    ])
    monkeypatch.setattr(cli, "read_pdb", lambda *a: object())
    monkeypatch.setattr(cli, "_validated_pairs", lambda *a: {})
    monkeypatch.setattr(cli, "annotate_fv", lambda *a: {})
    monkeypatch.setattr(cli, "structural_rmsds", lambda *a, **k: (Transform(np.eye(3), np.zeros(3)), []))
    monkeypatch.setattr(cli, "write_transformed_pdb", lambda *a: None)
    monkeypatch.setattr(cli, "read_residue_energy_breakdown", lambda *a: [{"residue1": "H:1", "residue2": "X:1"}])
    monkeypatch.setattr(cli, "merge_by_pair", lambda *a: [{"residue1": "H:1", "residue2": "X:1"}])
    monkeypatch.setattr(cli, "outer_compare", lambda *a: [{"parent_total": -1}])
    monkeypatch.setattr(cli, "mutation_network_changes", lambda *a, **k: pytest.fail("Network calculation executed"))
    cli.compare(args)
    assert (tmp_path / "pair_comparison.csv").is_file()
    assert not (tmp_path / "mutation_network_changes.csv").exists()
    assert json.loads((tmp_path / "summary.json").read_text())["analysis_status"]["network"] == "not_run"


def test_persistent_worker_omits_disabled_tools_and_their_cache_artifacts(tmp_path, monkeypatch):
    from fv_energy_analysis import fast_workflow

    candidate, params, template = [tmp_path / name for name in ("mutant.pdb", "lig.params", "lig.pdb")]
    candidate.write_text("END\n")
    params.write_text("validated params")
    template.write_text("END\n")
    cache = tmp_path / "cache"
    (cache / "final").mkdir(parents=True)
    prepared = cache / "parent.pdb"
    prepared.write_text("END\n")
    omitted = {"skip_arpeggio": True, "skip_plip": True, "skip_network_analysis": True,
               "skip_rotamer_analysis": True, "skip_e2_architecture": True}
    (cache / "workflow_manifest.json").write_text(json.dumps({
        "parent_analysis_input": str(prepared), "analysis_status": analysis_status(omitted),
    }))
    (cache / "final/summary.json").write_text(json.dumps({"parent": str(prepared)}))
    for name in ("parent_pair_energies.sc", "parent_interface.sc", "parent_hbonds.csv"):
        (cache / name).write_text("saved parent evidence")
    (cache / "ligand_atom_validation.json").write_text(json.dumps({"parent": {"valid": True}}))
    (cache / "estrogen_binding_architecture.json").write_text(json.dumps({"parent": {"e2_architecture_status": "not_run"}}))
    output = tmp_path / "output"
    config = dict(parent=str(prepared), candidate=str(candidate), output=str(output), parent_cache=str(cache),
                  params=str(params), ligand_template_pdb=str(template), ligand_resname="LIG", ligand_chain="X",
                  rosetta_ligand_chain="X", fv_chains="H,L", scheme="imgt", cdr_definition="imgt",
                  fit_region="fv", rmsd_atoms="backbone", rotamer_rmsd_cutoff=3,
                  network_energy_presence_cutoff=.1, network_energy_change_cutoff=.25, weights="ligand", **omitted)

    def run(runner, command, **kwargs):
        runner.commands.append(command)
        assert "--generated-ligand-pdb" in command, "Unexpected detector/conversion subprocess"
        folder = Path(command[command.index("--output") + 1])
        (folder / "mutant_input_rosetta.pdb").write_text("END\n")

    def score(pdb, pairs, interface, hbonds, *args, **kwargs):
        for path in (pairs, interface, hbonds):
            path.write_text("new mutant scoring evidence")

    def compare(args):
        assert args.parent_arpeggio is None and args.mutant_arpeggio is None
        folder = Path(args.output)
        folder.mkdir()
        (folder / "summary.json").write_text("{}")

    def forbidden(*args, **kwargs):
        pytest.fail("Disabled scientific stage executed")

    monkeypatch.setattr(fast_workflow.Runner, "run", run)
    monkeypatch.setattr(fast_workflow, "score_and_extract_hbonds", score)
    monkeypatch.setattr(fast_workflow, "compare", compare)
    monkeypatch.setattr(fast_workflow, "validate_ligand_against_params", lambda *a, **k: {"valid": True})
    monkeypatch.setattr(fast_workflow, "_write_automatic_summary", lambda *a: None)
    monkeypatch.setattr(fast_workflow.shutil, "which", forbidden)
    monkeypatch.setattr(fast_workflow, "run_arpeggio", forbidden)
    monkeypatch.setattr(fast_workflow, "summarize_plip_xml", forbidden)
    monkeypatch.setattr(fast_workflow, "analyze_estrogen_architecture", forbidden)
    result = fast_workflow.run_cached_candidate(config)
    assert result["status"] == "complete", result.get("error")
    manifest = json.loads((output / "workflow_manifest.json").read_text())
    assert set(manifest["analysis_status"].values()) == {"not_run"}
    assert not (output / "plip_ligand_interactions.json").exists()
    assert not list(output.glob("arpeggio_*/*.json"))
