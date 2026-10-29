"""Scoring choices are explicit, hash-bound and do not require an engine to inspect."""
from hashlib import sha256
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import MagicMock

import pytest

from fv_energy_analysis import hydrogen_batch, hydrogen_rescore, pyrosetta_score as scoring


@pytest.mark.parametrize("name,expected", [
    ("ligand", True), ("ligand.wts", True), ("ligand_soft_rep", True),
    ("ref2015", False), ("beta_nov16", False), ("custom.wts", False),
])
def test_historical_defaults_are_automatic_only_for_known_names(name, expected):
    assert scoring.resolve_legacy_pre_talaris(name) is expected


def test_explicit_historical_override_is_independent_of_custom_filename():
    assert scoring.resolve_legacy_pre_talaris("my_custom.wts", True) is True
    assert scoring.resolve_legacy_pre_talaris("ligand", False) is False
    with pytest.raises(ValueError, match="true, false or null"):
        scoring.resolve_legacy_pre_talaris("ref2015", "true")


def test_custom_file_identity_detects_in_place_weight_changes(tmp_path):
    weights = tmp_path / "my_custom.wts"
    weights.write_text("fa_atr 1.0\n")
    before = scoring.scoring_weights_identity(weights, True)
    assert before == {"weights": str(weights.resolve()), "legacy_pre_talaris": True,
                      "weights_sha256": sha256(weights.read_bytes()).hexdigest()}
    weights.write_text("fa_atr 0.8\n")
    assert scoring.scoring_weights_identity(weights, True)["weights_sha256"] != before["weights_sha256"]


def test_named_weights_can_be_inspected_without_pyrosetta():
    assert scoring.scoring_weights_identity("ref2015") == {
        "weights": "ref2015", "legacy_pre_talaris": False, "weights_sha256": None,
    }
    with pytest.raises(ValueError, match="weights must name"):
        scoring.scoring_weights_identity(" ")


def test_database_weight_identity_resolves_builtin_and_rejects_missing(tmp_path):
    weights = tmp_path / "scoring/weights"
    weights.mkdir(parents=True)
    (weights / "ligand.wts").write_text("fa_atr 0.8\n")
    result = scoring.scoring_weights_identity("ligand", database=tmp_path)
    assert result["weights"] == "ligand"
    assert result["legacy_pre_talaris"] is True
    assert result["weights_sha256"] == sha256((weights / "ligand.wts").read_bytes()).hexdigest()
    with pytest.raises(FileNotFoundError, match="were not found"):
        scoring.scoring_weights_identity("unknown_weights", database=tmp_path)
    with pytest.raises(FileNotFoundError, match="Custom Rosetta weights file"):
        scoring.scoring_weights_identity(tmp_path / "missing.wts")


def test_combined_scoring_forwards_explicit_legacy_override(monkeypatch, tmp_path):
    calls = []
    engine = SimpleNamespace(pose_from_file=lambda path: object())
    monkeypatch.setattr(scoring, "_initialize", lambda params, legacy: calls.append((params, legacy)) or engine)
    monkeypatch.setattr(scoring, "_scorefunction", lambda engine, weights: calls.append(weights) or object())
    monkeypatch.setattr(scoring, "write_pair_energies", lambda *args: None)
    monkeypatch.setattr(scoring, "write_interface_metrics", lambda *args: None)
    from fv_energy_analysis import hbonds
    monkeypatch.setattr(hbonds, "write_pose_hbonds", lambda *args, **kwargs: None)
    scoring.score_and_extract_hbonds("pose.pdb", "pairs.sc", "interface.sc", "hbonds.csv", "X",
                                     params=["E2L.params"], weights="custom.wts", legacy_pre_talaris=True)
    assert calls == [(["E2L.params"], True), "custom.wts"]


def test_persistent_scorefunction_rejects_in_place_custom_weight_changes(tmp_path, monkeypatch):
    weights = tmp_path / "custom.wts"
    weights.write_text("fa_atr 1.0\n")
    engine, score = MagicMock(), MagicMock()
    engine.create_score_function.return_value = score
    monkeypatch.setattr(scoring, "_SCORE_WEIGHTS_ENGINE", None)
    monkeypatch.setattr(scoring, "_SCORE_WEIGHTS_SIGNATURE", None)
    assert scoring._scorefunction(engine, str(weights)) is score
    assert scoring._scorefunction(engine, str(weights)) is score
    weights.write_text("fa_atr 0.8\n")
    with pytest.raises(RuntimeError, match="cannot change scoring weights"):
        scoring._scorefunction(engine, str(weights))
    assert engine.create_score_function.call_count == 2


def test_opth_worker_uses_selected_scorefunction_and_chemical_state_guards(monkeypatch):
    calls, guards = [], []
    engine = SimpleNamespace(rosetta=SimpleNamespace(basic=SimpleNamespace(options=SimpleNamespace(
        set_boolean_option=lambda name, value: guards.append((name, value))))))
    monkeypatch.setattr(scoring, "_initialize", lambda params, legacy: calls.append((params, legacy)) or engine)
    score = object()
    monkeypatch.setattr(scoring, "_scorefunction", lambda engine, weights: calls.append(weights) or score)
    monkeypatch.setattr(hydrogen_rescore, "_ENGINE", None)
    monkeypatch.setattr(hydrogen_rescore, "_SCORE", None)
    hydrogen_rescore._init_worker("E2L.params", "ref2015", False)
    assert calls == [(["E2L.params"], False), "ref2015"]
    assert guards == [("packing:no_optH", True), ("packing:flip_HNQ", False)]
    assert hydrogen_rescore._SCORE is score


def test_opth_protocol_records_actual_weights_and_explicit_legacy_choice(tmp_path, monkeypatch):
    module = tmp_path / "pyrosetta/__init__.py"
    module.parent.mkdir()
    module.write_text("")
    weights = module.parent / "database/scoring/weights"
    weights.mkdir(parents=True)
    (weights / "ligand.wts").write_text("fa_atr 0.8\n")
    (weights / "ref2015.wts").write_text("fa_atr 1.0\n")
    params, template = tmp_path / "E2L.params", tmp_path / "E2L.pdb"
    params.write_text("params")
    template.write_text("template")
    monkeypatch.setattr(hydrogen_batch.util, "find_spec", lambda name: SimpleNamespace(origin=str(module)))
    monkeypatch.setattr(hydrogen_batch.metadata, "version", lambda name: "test-version")
    default = hydrogen_batch._protocol(params, template, 20261022)
    assert default["source_protocol"]["weights"] == "ligand"
    assert default["source_protocol"]["legacy_pre_talaris"] is True
    selected = hydrogen_batch._protocol(params, template, 42, "ref2015")
    assert selected["source_protocol"]["weights"] == "ref2015"
    assert selected["source_protocol"]["legacy_pre_talaris"] is False
    assert selected["passes_per_complex"] == 1
    custom = tmp_path / "custom.wts"
    custom.write_text("fa_atr 0.9\n")
    historical = hydrogen_batch._protocol(params, template, 42, str(custom), True)
    assert historical["source_protocol"]["legacy_pre_talaris"] is True
    assert historical["source_protocol"]["weights_sha256"] == sha256(custom.read_bytes()).hexdigest()
