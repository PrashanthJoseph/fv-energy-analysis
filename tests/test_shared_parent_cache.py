"""Parent reuse must avoid rescoring without relaxing protocol validation."""
import io
import json
from hashlib import sha256
from pathlib import Path

import pytest

from fv_energy_analysis.batch_workflow import build_parser, _ensure_parent_cache, _parent_cache_key
from fv_energy_analysis.pyrosetta_score import SCORING_PROTOCOL_VERSION


def write_json(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value), encoding="utf-8")


@pytest.fixture
def shared_cache(tmp_path):
    original = tmp_path / "old_parent.pdb"
    original.write_text("ATOM      1  CA  ALA H   1       1.000   2.000   3.000  1.00  0.00           C\nEND\n")
    current = tmp_path / "new_parent.pdb"
    current.write_bytes(original.read_bytes())
    params, template = tmp_path / "E2L.params", tmp_path / "E2L_0001.pdb"
    params.write_text("validated params")
    template.write_text("validated template")
    output = tmp_path / "new_output"
    args = build_parser().parse_args([
        "--parent", str(current), "--input-folder", str(tmp_path), "--output", str(output),
        "--params", str(params), "--ligand-template-pdb", str(template),
        "--weights", "ligand", "--ligand-chain", "X", "--parent-ligand-chain", "X",
    ])
    cache = tmp_path / ("_parent_cache_" + _parent_cache_key(args, original))
    cache.mkdir()
    prepared = cache / "prepared.pdb"
    prepared.write_bytes(original.read_bytes())
    write_json(cache / "workflow_manifest.json", {
        "status": "complete", "parent": str(original), "parent_analysis_input": str(original),
        "scoring_protocol": {"version": SCORING_PROTOCOL_VERSION, "weights": "ligand",
                             "params_sha256": sha256(params.read_bytes()).hexdigest()},
    })
    write_json(cache / "automatic_summary.json", {})
    write_json(cache / "final/summary.json", {"parent": str(prepared)})
    for name in ("parent_pair_energies.sc", "parent_interface.sc", "parent_hbonds.csv",
                 "arpeggio_parent/parent_arpeggio.json", "plip_ligand_interactions.json",
                 "ligand_atom_validation.json", "estrogen_binding_architecture.json"):
        target = cache / name
        target.parent.mkdir(parents=True, exist_ok=True)
        target.write_text("saved artifact")
    args.reuse_parent_cache = str(cache)
    return args, current, output, cache


def test_reuses_existing_parent_directly_without_new_cache(shared_cache, monkeypatch):
    import fv_energy_analysis.batch_workflow as module
    args, current, output, cache = shared_cache
    monkeypatch.setattr(module.subprocess, "run", lambda *a, **k: pytest.fail("Parent was rescored"))
    result = _ensure_parent_cache(args, current, output, io.StringIO())
    assert result == cache
    assert not list(output.glob("_parent_cache_*"))
    provenance = json.loads((output / "shared_parent_cache_reuse.json").read_text())
    assert provenance["parent_scoring_repeated"] is False
    assert provenance["parent_sha256"] == sha256(current.read_bytes()).hexdigest()


def test_rejects_different_parent(shared_cache):
    args, current, output, _ = shared_cache
    current.write_text(current.read_text().replace("1.000", "4.000"))
    with pytest.raises(ValueError, match="different parent"):
        _ensure_parent_cache(args, current, output, io.StringIO())


def test_rejects_changed_parameters_or_scoring_settings(shared_cache):
    args, current, output, _ = shared_cache
    args.weights = "ref2015"
    with pytest.raises(ValueError, match="settings/code differ"):
        _ensure_parent_cache(args, current, output, io.StringIO())


def test_rejects_missing_required_parent_artifacts(shared_cache):
    args, current, output, cache = shared_cache
    (cache / "parent_hbonds.csv").unlink()
    with pytest.raises(FileNotFoundError, match="artifacts are missing"):
        _ensure_parent_cache(args, current, output, io.StringIO())


def test_rejects_changed_prepared_coordinates(shared_cache):
    args, current, output, cache = shared_cache
    prepared = cache / "prepared.pdb"
    prepared.write_text(prepared.read_text().replace("1.000", "4.000"))
    with pytest.raises(ValueError, match="prepared parent differs"):
        _ensure_parent_cache(args, current, output, io.StringIO())


def controller_module():
    from fv_energy_analysis import full_pipeline as module
    return module


def test_controller_forwards_shared_parent_cache():
    controller = controller_module()
    config = {"workspace": "D:/analysis", "input_folder": "D:/inputs", "parent": "D:/parent.pdb",
              "output": "D:/out", "params": "D:/E2L.params", "ligand_template_pdb": "D:/E2L_0001.pdb",
              "jobs": 3, "distro": "Ubuntu", "wsl_python": "/home/user/bin/python",
              "reuse_parent_cache": "D:/saved/_parent_cache_123", "weights": "ligand"}
    assert "--reuse-parent-cache /mnt/d/saved/_parent_cache_123" in controller.batch_command(config)[-1]


def test_controller_can_skip_intermediate_xlsx(tmp_path, monkeypatch):
    controller = controller_module()
    output = tmp_path / "output"
    output.mkdir()
    write_json(output / "batch_manifest.json", {"status": "complete", "candidate_count": 1,
                                               "completed_count": 1, "failure_count": 0})
    params = tmp_path / "params"
    params.write_text("params")
    commands = []
    monkeypatch.setattr(controller.subprocess, "run", lambda command, **kwargs:
                        (commands.append(command) or type("Completed", (), {"returncode": 0})()))
    config = {"workspace": str(Path(__file__).resolve().parents[1]), "output": str(output),
              "params": str(params), "postprocessing_only": True, "export_xlsx": False, "jobs": 3,
              "node": "node", "marker_script": "marker", "export_builder": "builder"}
    controller.run(config)
    assert not any("marker" in command or "builder" in command for command in commands)
    status = json.loads((output / "full_pipeline_status.json").read_text())
    assert status["status"] == "complete"
    assert status["workbook"] is None and status["expected_xlsx"] is None
