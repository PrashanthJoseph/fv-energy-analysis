"""A saved-batch rescore must reuse only validated prepared structures."""

import csv
import hashlib
import json
from pathlib import Path
import sys
from types import SimpleNamespace
from unittest.mock import Mock

import pytest

from fv_energy_analysis.rescore_saved import collect_inputs


def write_json(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value), encoding="utf-8")


def saved_batch(tmp_path, statuses=("complete", "complete", "failed")):
    batch = tmp_path / "saved_batch"
    batch.mkdir()
    params = tmp_path / "E2L.params"
    params.write_text("NAME E2L\nIO_STRING E2L X\n", encoding="utf-8")
    digest = hashlib.sha256(params.read_bytes()).hexdigest()
    rows, cases = [], []
    for index, status in enumerate(statuses):
        folder = batch / f"candidate_{index}"
        parent, mutant = folder / "prepared" / "parent.pdb", folder / "prepared" / "mutant.pdb"
        parent.parent.mkdir(parents=True)
        parent.write_text("REMARK same prepared parent\nEND\n", encoding="utf-8")
        mutant.write_text(f"REMARK prepared candidate {index}\nEND\n", encoding="utf-8")
        write_json(folder / "final" / "summary.json", {
            "parent": str(parent), "mutant": str(mutant),
            "fv_chains": ["H", "L"], "ligand_chains": ["X"],
        })
        write_json(folder / "workflow_manifest.json", {
            "status": status, "params": str(params), "ligand_resname": "E2L",
            "scoring_protocol": {"version": 2, "weights": "ref2015", "params_sha256": digest},
        })
        rows.append({
            "candidate": f"candidate_{index}.pdb", "result_folder": str(folder),
            "candidate_path": str(tmp_path / "nonexistent_raw_input.pdb"),
            "status": status, "substitutions": "E99W",
        })
        cases.append({"folder": folder, "parent": parent, "pdb": mutant})
    with (batch / "batch_summary.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    return batch, params, cases


def test_collects_complete_prepared_structures_not_raw_inputs(tmp_path):
    batch, params, cases = saved_batch(tmp_path)
    result = collect_inputs(batch)
    assert result["params"] == params.resolve()
    assert result["parent"] == cases[0]["parent"].resolve()
    assert result["ligand_chain"] == "X"
    assert result["source_count"] == 3
    assert result["excluded_count"] == 1
    assert [item["candidate"] for item in result["candidates"]] == ["candidate_0.pdb", "candidate_1.pdb"]
    assert [item["pdb"] for item in result["candidates"]] == [case["pdb"].resolve() for case in cases[:2]]


def test_rejected_cases_need_no_prepared_artifacts(tmp_path):
    batch, _, cases = saved_batch(tmp_path, ("complete", "failed"))
    cases[1]["parent"].unlink()
    cases[1]["pdb"].unlink()
    (cases[1]["folder"] / "final" / "summary.json").unlink()
    result = collect_inputs(batch)
    assert len(result["candidates"]) == 1


def test_rejects_different_prepared_parent_contents(tmp_path):
    batch, _, cases = saved_batch(tmp_path)
    cases[1]["parent"].write_text("REMARK different prepared parent\nEND\n", encoding="utf-8")
    with pytest.raises(ValueError, match="(?i)parent"):
        collect_inputs(batch)


def test_rejects_changed_parameter_contents(tmp_path):
    batch, params, _ = saved_batch(tmp_path)
    params.write_text("NAME CHANGED\n", encoding="utf-8")
    with pytest.raises(ValueError, match="(?i)(param|hash|chemistry)"):
        collect_inputs(batch)


def test_rejects_candidate_with_different_recorded_parameter_hash(tmp_path):
    batch, _, cases = saved_batch(tmp_path)
    manifest_path = cases[1]["folder"] / "workflow_manifest.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    manifest["scoring_protocol"]["params_sha256"] = "0" * 64
    write_json(manifest_path, manifest)
    with pytest.raises(ValueError, match="(?i)(param|hash|chemistry)"):
        collect_inputs(batch)


def test_missing_prepared_mutant_fails_before_scoring(tmp_path):
    batch, _, cases = saved_batch(tmp_path)
    cases[0]["pdb"].unlink()
    with pytest.raises((FileNotFoundError, ValueError)):
        collect_inputs(batch)


def test_rejects_batch_with_no_completed_candidates(tmp_path):
    batch, _, _ = saved_batch(tmp_path, ("failed",))
    with pytest.raises(ValueError, match="(?i)(complete|candidate|success)"):
        collect_inputs(batch)


@pytest.mark.parametrize("weights,legacy", [
    ("ligand", True), ("ligand.wts", True), ("ligand_soft_rep", True),
    ("/opt/rosetta/scoring/weights/ligand.wts", True),
    ("ref2015", False), ("beta_nov16", False),
])
def test_historical_weight_files_require_legacy_initialization(weights, legacy):
    from fv_energy_analysis.pyrosetta_score import _legacy_weights
    assert _legacy_weights(weights) is legacy


def test_initialize_sets_legacy_flag_once_without_real_pyrosetta(tmp_path, monkeypatch):
    from fv_energy_analysis import pyrosetta_score as scoring
    params = tmp_path / "lig.params"
    params.write_text("NAME LIG\n", encoding="utf-8")
    fake = SimpleNamespace(init=Mock())
    monkeypatch.setitem(sys.modules, "pyrosetta", fake)
    monkeypatch.setattr(scoring, "_PYROSETTA", None)
    monkeypatch.setattr(scoring, "_INITIALIZED_PARAMS", None)
    monkeypatch.setattr(scoring, "_INITIALIZED_PARAM_DIGESTS", None)
    monkeypatch.setattr(scoring, "_INITIALIZED_LEGACY", None)
    assert scoring._initialize([str(params)], legacy=True) is fake
    assert scoring._initialize([str(params)], legacy=True) is fake
    fake.init.assert_called_once()
    options = fake.init.call_args.args[0]
    assert "-restore_pre_talaris_2013_behavior true" in options
    assert "-extra_res_fa" in options
    with pytest.raises(RuntimeError, match="separate PyRosetta processes"):
        scoring._initialize([str(params)], legacy=False)


def test_persistent_legacy_worker_rejects_changed_params(tmp_path, monkeypatch):
    from fv_energy_analysis import pyrosetta_score as scoring
    params = tmp_path / "lig.params"
    params.write_text("NAME LIG\n", encoding="utf-8")
    monkeypatch.setattr(scoring, "_PYROSETTA", object())
    monkeypatch.setattr(scoring, "_INITIALIZED_PARAMS", (str(params.resolve()),))
    monkeypatch.setattr(scoring, "_INITIALIZED_PARAM_DIGESTS", (hashlib.sha256(params.read_bytes()).hexdigest(),))
    monkeypatch.setattr(scoring, "_INITIALIZED_LEGACY", True)
    params.write_text("NAME ALTERED\n", encoding="utf-8")
    with pytest.raises(RuntimeError, match="cannot change ligand params"):
        scoring._initialize([str(params)], legacy=True)
