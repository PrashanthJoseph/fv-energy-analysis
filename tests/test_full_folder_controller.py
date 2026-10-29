"""Controller configuration checks without launching scientific work or Excel."""
from __future__ import annotations

import json
from pathlib import Path
import subprocess
import sys

import pytest


from fv_energy_analysis import full_pipeline as controller

MODULE_PATH = Path(controller.__file__).resolve()


def command_config(**updates):
    config = {
        "workspace": "D:/analysis", "input_folder": "D:/inputs", "parent": "D:/parent.pdb",
        "output": "D:/results", "params": "D:/E2L.params",
        "ligand_template_pdb": "D:/E2L_0001.pdb", "jobs": 3,
        "wsl_python": "/home/user/venv/bin/python", "distro": "Ubuntu-24.04",
    }
    return {**config, **updates}


def test_ligand_weights_with_canonical_chains_do_not_swap():
    command = controller.batch_command(command_config(weights="ligand", mutant_chain_map=""))
    script = command[-1]
    assert "--weights ligand" in script
    assert "--jobs 3" in script
    assert "--mutant-chain-map" not in script
    assert "--parent-chain-map" not in script
    assert "--arpeggio-scope neighbourhood" in script
    assert "OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1" in script


def test_controller_imports_package_without_parent_pythonpath():
    script = (
        "import importlib.util; "
        f"spec=importlib.util.spec_from_file_location('controller',{str(MODULE_PATH)!r}); "
        "module=importlib.util.module_from_spec(spec); spec.loader.exec_module(module); "
        f"module.configure_workspace_imports({{'workspace':{str(MODULE_PATH.parents[2])!r}}}); "
        "from fv_energy_analysis.report import write_json; print('package import succeeded')"
    )
    result = subprocess.run([sys.executable, "-I", "-c", script], capture_output=True, text=True, timeout=30)
    assert result.returncode == 0, result.stderr
    assert "package import succeeded" in result.stdout


def test_fresh_folder_opth_is_one_pass_before_analysis():
    script = controller.hydrogen_command(command_config(
        hydrogen_parent="D:/raw_parent.pdb", hydrogen_output="D:/optimized",
        reuse_parent_result="D:/cached_parent/result.json", expected_count=845,
    ))[-1]
    assert "fv_energy_analysis.hydrogen_batch" in script
    assert "--parent /mnt/d/raw_parent.pdb" in script
    assert "--reuse-parent-result /mnt/d/cached_parent/result.json" in script
    assert "--template /mnt/d/E2L_0001.pdb" in script
    assert "--expected-count 845" in script
    assert "--seed 20261022" in script
    assert "--ligand-chain d" in script


def test_opth_retains_explicit_protein_chain_maps_but_normalizes_ligand_chain():
    config = controller.optimized_analysis_config(command_config(
        hydrogen_output="D:/optimized", mutant_chain_map="A:H,B:L", parent_chain_map="C:H,D:L",
        ligand_chain="d", mutant_ligand_chain="d", parent_ligand_chain="d"))
    script = controller.batch_command(config)[-1]
    assert "--mutant-chain-map A:H,B:L" in script
    assert "--parent-chain-map C:H,D:L" in script
    assert "--parent-ligand-chain X" in script and "--mutant-ligand-chain X" in script


def test_combination_refreshes_only_saved_architecture_before_merging():
    config = command_config(
        architecture_refresh={"source": "original1397", "output": "refreshed1397",
                              "opth_results": "savedOptH", "expected_count": 1397},
        c18_calibration="calibration.json", combine_sources=["refreshed1397", "completed845", "remainder194"],
        combined_output="combined2436", combined_expected_count=2436,
        combined_expected_candidates="all_candidates.json",
        node="node.exe", marker_script="marker.mjs", export_builder="export_full.mjs",
    )
    commands = list(controller.finalization_commands(config))
    assert [name for name, _ in commands] == ["refresh_previous_architecture", "combine_full_csv"]
    assert "fv_energy_analysis.architecture_batch" in commands[0][1]
    assert "--opth-results" in commands[0][1]
    assert "--xlsx" not in commands[0][1]
    assert commands[1][1].count("--source") == 3
    assert "--expected-count" in commands[1][1] and "2436" in commands[1][1]
    assert "--expected-candidates" in commands[1][1]


def test_legacy_configuration_requires_explicit_weights_and_chain_swap():
    script = controller.batch_command(command_config(weights="ref2015", mutant_chain_map="H:L,L:H"))[-1]
    assert "--weights ref2015" in script
    assert "--mutant-chain-map H:L,L:H" in script


def test_production_defaults_are_ligand_and_no_chain_swap_or_recursive_scan():
    script = controller.batch_command(command_config())[-1]
    assert "--weights ligand" in script
    assert "--mutant-chain-map" not in script
    assert "--recursive" not in script


def test_explicit_chain_maps_and_ligand_chains_are_forwarded():
    script = controller.batch_command(command_config(
        mutant_chain_map="A:H,B:L", parent_chain_map="C:H,D:L", ligand_chain="d",
        parent_ligand_chain="P", mutant_ligand_chain="M", rosetta_ligand_chain="X",
    ))[-1]
    for flag, value in (("--mutant-chain-map", "A:H,B:L"), ("--parent-chain-map", "C:H,D:L"),
                        ("--ligand-chain", "d"), ("--parent-ligand-chain", "P"),
                        ("--mutant-ligand-chain", "M"), ("--rosetta-ligand-chain", "X")):
        assert f"{flag} {value}" in script


def test_missing_optional_metadata_does_not_forward_old_dataset_paths():
    command = controller.fingerprint_command(command_config())
    assert command[-2:] == ["--batch-results", "D:/results"]
    assert "--cif-folder" not in command
    assert "--parent-cif" not in command
    assert "--scores-csv" not in command


def test_supplied_optional_metadata_is_forwarded():
    command = controller.fingerprint_command(command_config(
        cif_folder="D:/cifs", parent_cif="D:/parent.cif", scores_csv="D:/scores.csv",
    ))
    assert command[-6:] == ["--cif-folder", "D:/cifs", "--parent-cif", "D:/parent.cif",
                           "--scores-csv", "D:/scores.csv"]


def test_explicit_empty_metadata_sections_are_forwarded():
    command = controller.fingerprint_command(command_config(include_empty_metadata_sections=True))
    assert command[-1] == "--include-empty-metadata-sections"


def saved_config(tmp_path, **updates):
    values = {}
    for key, kind in controller.PATH_KINDS.items():
        path = tmp_path / key
        if kind == "dir":
            path.mkdir()
        else:
            path.write_text("fixture", encoding="utf-8")
        values[key] = str(path)
    values.update(output=str(tmp_path / "output"), distro="Ubuntu-24.04", skip_discovery_studio=True,
                  wsl_python="/home/user/venv/bin/python")
    values.update(updates)
    path = tmp_path / "config.json"
    path.write_text(json.dumps(values), encoding="utf-8")
    return path


def test_read_config_accepts_absent_optional_metadata(tmp_path):
    config = controller.read_config(saved_config(tmp_path, weights="ligand", mutant_chain_map=""))
    assert config["weights"] == "ligand"
    assert config["mutant_chain_map"] == ""
    assert "cif_folder" not in config
    assert "parent_cif" not in config
    assert "scores_csv" not in config


def test_config_requires_no_node_or_temporary_exporter_and_defaults_to_opth(tmp_path):
    config = controller.read_config(saved_config(tmp_path))
    assert config["weights"] == "ligand"
    assert config["hydrogen_parent"] == config["parent"]
    assert Path(config["hydrogen_output"]).name == "output_opth"
    assert not {"node", "export_builder", "marker_script"} & config.keys()


def test_native_linux_science_uses_configured_python_without_wsl_or_remapping():
    config = command_config(science_backend="native", science_python="/venv/bin/python",
                            parent="/inputs/parent.pdb", input_folder="/inputs", output="/outputs",
                            params="/chem/E2L.params", ligand_template_pdb="/chem/E2L.pdb")
    command = controller.batch_command(config)
    assert command[:4] == ["/venv/bin/python", "-u", "-m", "fv_energy_analysis.batch_workflow"]
    assert command[command.index("--parent") + 1] == "/inputs/parent.pdb"
    assert "--mutant-chain-map" not in command and "wsl.exe" not in command


def test_optional_e3_reuses_existing_modules_and_persisted_prepared_coordinates():
    config = command_config(e3_discrimination=True)
    commands = list(controller.finalization_commands(config))
    assert [name for name, _ in commands] == ["e3_source_inventory", "e3_steric_discrimination"]
    assert "fv_energy_analysis.combine_batches" in commands[0][1]
    assert "--include-failures" in commands[0][1]
    assert "fv_energy_analysis.e3_discrimination" in commands[1][1][-1]
    assert "--resume" in commands[1][1][-1]
    assert not any("hydrogen_batch" in str(command) for _, command in commands)


def test_native_linux_requires_explicit_partial_run_when_discovery_studio_unavailable(tmp_path):
    if sys.platform == "win32":
        pytest.skip("Linux-specific native DS gate")
    with pytest.raises(ValueError, match="requires a Windows controller"):
        controller.read_config(saved_config(tmp_path, skip_discovery_studio=False))


def test_read_config_checks_optional_paths_when_supplied(tmp_path):
    with pytest.raises(FileNotFoundError, match="cif_folder"):
        controller.read_config(saved_config(tmp_path, cif_folder=str(tmp_path / "missing_cifs")))


def test_parent_cif_requires_mutant_cif_folder(tmp_path):
    parent_cif = tmp_path / "parent.cif"
    parent_cif.write_text("data_parent", encoding="utf-8")
    with pytest.raises(ValueError, match="parent_cif requires cif_folder"):
        controller.read_config(saved_config(tmp_path, parent_cif=str(parent_cif)))


def test_empty_scoring_function_is_rejected(tmp_path):
    with pytest.raises(ValueError, match="weights must name"):
        controller.read_config(saved_config(tmp_path, weights="  "))


def test_preflight_checks_actual_science_and_reporting_interpreters():
    commands = list(controller.preflight_commands(command_config()))
    assert [name for name, _ in commands] == ["science_dependencies", "reporting_dependencies"]
    assert "/home/user/venv/bin/python" in commands[0][1][-1]
    assert "--profile analysis" in commands[0][1][-1]
    assert commands[1][1][-2:] == ["--profile", "reporting"]


def combination_fixture(tmp_path):
    output = tmp_path / "combined"
    output.mkdir()
    manifest = {"status": "complete", "candidate_count": 1, "completed_count": 1,
                "failure_count": 0, "columns": 1, "scientific_protocol": {"weights": "ligand"}}
    (output / "combined_manifest.json").write_text(json.dumps(manifest))
    for name, contents in (("batch_summary_grouped_fingerprints.json", '{"rows": [{"candidate": "a.pdb"}]}'),
                           ("batch_summary_grouped_fingerprints.csv", "candidate\na.pdb\n"),
                           ("candidate_sources.csv", "candidate,source\na.pdb,source\n")):
        (output / name).write_text(contents)
    identity = {"sources": [{"payload_sha256": "verified source", "scientific_protocol": {"weights": "ligand"}}]}
    (output / "full_pipeline_combination_identity.json").write_text(json.dumps({
        "identity": identity, "artifacts": controller.combination_artifacts(output)}))
    return output, identity


def test_completed_combination_reused_only_when_sources_and_saved_reports_match(tmp_path):
    output, identity = combination_fixture(tmp_path)
    assert controller.reuse_combination(output, identity) is True
    # Later workbook export may add presentation metadata, but not scientific values.
    manifest_path = output / "combined_manifest.json"
    manifest = json.loads(manifest_path.read_text())
    manifest["interaction_fingerprints"] = {"grouped_summary_xlsx": "report.xlsx"}
    manifest_path.write_text(json.dumps(manifest))
    assert controller.reuse_combination(output, identity) is True


def test_changed_combination_source_or_saved_values_require_new_output(tmp_path):
    output, identity = combination_fixture(tmp_path)
    with pytest.raises(ValueError, match="has changed"):
        controller.reuse_combination(output, {"sources": []})
    (output / "batch_summary_grouped_fingerprints.csv").write_text("candidate\nchanged.pdb\n")
    with pytest.raises(ValueError, match="has changed"):
        controller.reuse_combination(output, identity)


def test_unbound_existing_combination_is_not_silently_overwritten(tmp_path):
    fresh = tmp_path / "fresh"
    assert controller.reuse_combination(fresh, {}) is False
    fresh.mkdir()
    (fresh / "unrelated.txt").write_text("keep")
    with pytest.raises(FileExistsError, match="no verified controller provenance"):
        controller.reuse_combination(fresh, {})


def test_explicit_linux_partial_run_records_unavailable_discovery_studio(tmp_path, monkeypatch):
    output = tmp_path / "output"
    output.mkdir()
    (output / "batch_manifest.json").write_text(json.dumps({
        "status": "complete", "candidate_count": 1, "completed_count": 1, "failure_count": 0}))
    params = tmp_path / "params"
    params.write_text("params")
    commands = []
    monkeypatch.setattr(controller.subprocess, "run", lambda command, **kwargs:
                        (commands.append(command) or type("Completed", (), {"returncode": 0})()))
    controller.run({"output": str(output), "params": str(params), "jobs": 1,
                    "postprocessing_only": True, "skip_discovery_studio": True, "export_xlsx": False})
    status = json.loads((output / "full_pipeline_status.json").read_text())
    assert status["status"] == "complete_with_unavailable_tools"
    assert status["unavailable_tools"] == ["Discovery Studio"]
    assert not any("fv_energy_analysis.discovery_studio" in command for command in commands)
    assert json.loads((output / "batch_manifest.json").read_text())["discovery_studio"]["status"] == "unavailable"


def test_preparation_failure_kept_in_master_counts_and_rows(tmp_path):
    (tmp_path / "batch_manifest.json").write_text(json.dumps({
        "candidate_count": 1, "completed_count": 1, "failure_count": 0,
        "candidates": [{"candidate": "valid.pdb", "status": "complete"}],
    }))
    (tmp_path / "batch_summary.csv").write_text("candidate,status\nvalid.pdb,complete\n")
    hydrogen = {"candidate_count": 2, "input_folder": "raw", "cases": [
        {"candidate": "valid.pdb", "status": "complete"},
        {"candidate": "rejected.pdb", "status": "failed", "error_type": "LigandIdentityError",
         "error": "stereochemistry mismatch", "diagnostic_folder": "diagnostics"},
    ]}
    controller.merge_preparation_failures({"output": str(tmp_path)}, hydrogen)
    manifest = json.loads((tmp_path / "batch_manifest.json").read_text())
    assert manifest["candidate_count"] == 2
    assert manifest["completed_count"] == 1
    assert manifest["failure_count"] == 1
    assert manifest["candidates"][-1]["stage"] == "hydrogen_optimization"
    assert "stereochemistry mismatch" in (tmp_path / "batch_summary.csv").read_text()
    assert (tmp_path / "preparation_failures.json").is_file()
    original = (tmp_path / "batch_summary.csv").read_bytes()
    controller.merge_preparation_failures({"output": str(tmp_path)}, hydrogen)
    assert (tmp_path / "batch_summary.csv").read_bytes() == original
    assert json.loads((tmp_path / "batch_manifest.json").read_text())["failure_count"] == 1


def test_all_rejected_inputs_keep_parent_but_no_mutant_scores(tmp_path):
    reference, output = tmp_path / "reference", tmp_path / "output"
    reference.mkdir()
    output.mkdir()
    parent = tmp_path / "parent.pdb"
    parent.write_text("same optimized parent")
    (reference / "batch_manifest.json").write_text(json.dumps({"parent": str(parent)}))
    payload = {"columns": [{"key": "candidate"}, {"key": "energy"}],
               "parent": {"candidate": "Parent", "energy": -10},
               "rows": [{"candidate": "previous.pdb", "energy": -12}]}
    (reference / "batch_summary_grouped_fingerprints.json").write_text(json.dumps(payload))
    hydrogen = {"candidate_count": 1, "parent_optimized_pdb": str(parent), "input_folder": "raw",
                "cases": [{"candidate": "invalid.pdb", "status": "failed", "error_type": "LigandIdentityError",
                           "error": "stereochemistry mismatch", "diagnostic_folder": "diagnostics"}]}
    controller.initialize_rejected_cohort({"output": str(output), "failure_reference_batch": str(reference),
        "params": "E2L.params", "ligand_template_pdb": "template.pdb", "hydrogen_output": "opth"}, hydrogen)
    manifest = json.loads((output / "batch_manifest.json").read_text())
    result = json.loads((output / "batch_summary_grouped_fingerprints.json").read_text())
    assert manifest["preparation_only"] is True
    assert manifest["completed_count"] == 0 and manifest["failure_count"] == 1
    assert result["rows"] == []
    assert result["parent"] == payload["parent"]
    assert "mutant scientific scores" in result["summary_subtitle"]
