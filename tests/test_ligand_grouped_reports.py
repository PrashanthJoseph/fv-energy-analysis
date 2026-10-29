"""Legacy-weight grouped reports reuse saved data and keep provenance explicit."""

import csv
import json

import pytest

from fv_energy_analysis import fingerprint_reports as reports
from fv_energy_analysis.confidence import confidence_columns
from fv_energy_analysis.fingerprint_counts import (
    _empty, _rosetta_energies, _row,
)
from fv_energy_analysis.model_scores import score_columns
from fv_energy_analysis.report import write_csv


LEGACY_KEYS = {"rosetta_fa_pair_REU", "rosetta_hbond_bb_sc_REU"}


def _saved_fingerprints(candidate, result_dir, _parent_ds, summary=None):
    """Representative saved detector values; no structure files are required."""
    result = []
    for tool in reports.ESSENTIAL_COUNTS:
        parent, mutant = _empty(), _empty()
        for _, name in reports.ESSENTIAL_COUNTS[tool]:
            parent[name], mutant[name] = 1, 2
        energies = {}
        if tool == "Rosetta":
            for label, pair, bb_sc in (("parent", -0.1, -0.2), ("mutant", -0.4, -0.5)):
                energies.update({f"{label}_{name}": 0 for name in (
                    "fa_atr_REU", "fa_rep_REU", "fa_elec_REU", "solvation_REU", "hbond_sc_REU",
                )})
                energies.update({f"{label}_fa_pair_REU": pair, f"{label}_hbond_bb_sc_REU": bb_sc})
        result.append(_row(candidate, result_dir, tool, "saved records", parent, mutant, energies))
    return result


def _batch(tmp_path, monkeypatch, weights=("ligand",)):
    """Create only report inputs and intercept structure/detector operations."""
    bases = []
    for index, profile in enumerate(weights):
        directory = tmp_path / f"case_{index}"
        directory.mkdir()
        (directory / "automatic_summary.json").write_text("{}")
        (directory / "automatic_summary.md").write_text("# Saved analysis\n")
        manifest = {} if profile is None else {"scoring_protocol": {"weights": profile}}
        (directory / "workflow_manifest.json").write_text(json.dumps(manifest))
        bases.append({
            "candidate": f"case_{index}.pdb", "status": "complete", "result_folder": str(directory),
            "substitutions": "S89A", "rank_by_interface_delta": index + 1,
            "interface_parent_REU": -10, "interface_mutant_REU": -11,
            "interface_delta_REU": -1,
        })
    write_csv(tmp_path / "batch_summary.csv", bases)
    (tmp_path / "batch_manifest.json").write_text(json.dumps({"candidate_count": len(bases)}))
    monkeypatch.setattr(reports, "candidate_fingerprints", _saved_fingerprints)
    monkeypatch.setattr(reports, "_raw_count_rows", lambda *args, **kwargs: [])
    monkeypatch.setattr(reports, "_saved_architecture", lambda *args, **kwargs: {
        "e2_architecture_status": "unavailable: fixture has no geometry",
    })
    from fv_energy_analysis import arpeggio, pyrosetta_score

    def forbidden(*args, **kwargs):
        pytest.fail("Saved reporting must not run scoring or interaction detectors")

    monkeypatch.setattr(arpeggio, "run_arpeggio", forbidden)
    monkeypatch.setattr(pyrosetta_score, "_initialize", forbidden)
    return tmp_path


def _payload(batch):
    return json.loads((batch / "batch_summary_grouped_fingerprints.json").read_text())


def test_legacy_energy_terms_use_recorded_chain_and_preserve_missing_values():
    summary = {"ligand_chain": "Z", "interface_terms": [
        {"metric": "if_Z_fa_pair", "parent": -0.1, "mutant": -0.4,
         "delta_mutant_minus_parent": -0.3},
        {"metric": "if_Z_hbond_bb_sc", "parent": -0.2, "mutant": -0.5,
         "delta_mutant_minus_parent": -0.3},
    ]}
    values = _rosetta_energies(summary)
    assert values["parent_fa_pair_REU"] == -0.1
    assert values["mutant_hbond_bb_sc_REU"] == -0.5
    assert values["delta_fa_pair_REU"] == -0.3
    assert _rosetta_energies({})["parent_fa_pair_REU"] is None
    assert _rosetta_energies({})["mutant_hbond_bb_sc_REU"] is None


def test_custom_weights_honor_recorded_legacy_options(tmp_path, monkeypatch):
    batch = _batch(tmp_path, monkeypatch, weights=("custom_ligand.wts",))
    manifest = batch / "case_0/workflow_manifest.json"
    saved = json.loads(manifest.read_text())
    saved["scoring_protocol"]["legacy_pre_talaris"] = True
    manifest.write_text(json.dumps(saved))
    reports.generate(batch)
    payload = _payload(batch)
    assert payload["scoring"]["weights"] == "custom_ligand.wts"
    assert payload["scoring"]["legacy_pre_talaris"] is True
    assert LEGACY_KEYS <= {column["key"] for column in payload["columns"]}


def test_disabled_metadata_ignores_saved_sources_and_e2_is_not_computed(tmp_path, monkeypatch):
    batch = _batch(tmp_path, monkeypatch)
    (batch / "report_sources.json").write_text(json.dumps({
        "scores_csv": str(batch / "missing_scores.csv"), "scores_csv_sha256": "old_hash",
        "cif_folder": str(batch / "missing_cifs"), "parent_cif": str(batch / "missing_parent.cif")}))
    monkeypatch.setattr(reports, "_saved_architecture", lambda *a, **k: pytest.fail("E2 analysis was disabled"))
    reports.generate(batch, include_empty_metadata_sections=True, analysis_options={
        "e2_architecture": False, "hcdr3_plddt": False, "boltz2_scores": False, "gnina_scores": False})
    payload = _payload(batch)
    assert payload["rows"][0]["e2_architecture_status"] == "not_run"
    assert all(payload["rows"][0][column["key"]] is None for column in confidence_columns() + score_columns())
    assert set(payload["metadata_availability"].values()) == {"not_run"}


def test_legacy_schema_adds_only_active_extra_terms_and_preserves_default():
    default = reports._grouped_columns()
    ligand = reports._grouped_columns(weights="ligand")
    assert not LEGACY_KEYS & {column["key"] for column in default}
    assert LEGACY_KEYS <= {column["key"] for column in ligand}
    assert len(ligand) == len(default) + 2
    assert {column["section"] for column in ligand} == {
        "Identity", "RosettaLigand", "Arpeggio", "PLIP", "Discovery Studio",
        "Structure & network", "E2 architecture",
    }
    unchanged = [column for column in ligand if column["section"] != "RosettaLigand"]
    assert unchanged == [column for column in default if column["section"] != "Rosetta"]
    assert reports._grouped_columns() == default  # No in-place mutation of shared definitions.


@pytest.mark.parametrize("parent, expected", [(True, (-0.1, -0.2)), (False, (-0.4, -0.5))])
def test_grouped_row_exposes_legacy_terms_only_for_legacy_profile(tmp_path, parent, expected):
    base = {"candidate": "mutant.pdb", "rank_by_interface_delta": 1,
            "interface_parent_REU": -10, "interface_mutant_REU": -11, "interface_delta_REU": -1}
    fingerprints = _saved_fingerprints("mutant.pdb", tmp_path, None)
    legacy = reports._grouped_row(base, fingerprints, parent=parent, weights="ligand")
    assert (legacy["rosetta_fa_pair_REU"], legacy["rosetta_hbond_bb_sc_REU"]) == expected
    original = reports._grouped_row(base, fingerprints, parent=parent)
    assert not LEGACY_KEYS & original.keys()


def test_fresh_ligand_report_retains_full_layout_with_explicitly_missing_metadata(tmp_path, monkeypatch):
    batch = _batch(tmp_path, monkeypatch)
    reports.generate(batch, include_empty_metadata_sections=True)
    payload = _payload(batch)
    assert payload["scoring_profile"] == "rosettaligand_legacy"
    assert len(payload["columns"]) == 110
    sections = list(dict.fromkeys(column["section"] for column in payload["columns"]))
    assert sections == ["Identity", "RosettaLigand", "Arpeggio", "PLIP", "Discovery Studio",
                        "Structure & network", "E2 architecture", "HCDR3 pLDDT", "Boltz-2", "GNINA"]
    # The full layout retains both HCDR3 partitions and the current enclosure
    # measurements, including when no model-confidence metadata is supplied.
    keys = {column["key"] for column in payload["columns"]}
    assert {"h_cdr3_anchors_rmsd_A", "h_cdr3_tip_rmsd_A",
            "e2_steric_opposing_residue_pair_count",
            "e2_steric_overlapping_atom_pair_count"} <= keys
    metadata = confidence_columns() + score_columns()
    for row in [payload["parent"], *payload["rows"]]:
        assert all(row[column["key"]] is None for column in metadata)
    assert payload["rows"][0]["rosetta_fa_pair_REU"] == -0.4
    assert payload["parent"]["rosetta_hbond_bb_sc_REU"] == -0.2
    assert payload["metadata_availability"]
    with (batch / "batch_summary_grouped_fingerprints.csv").open(newline="") as handle:
        rows = list(csv.DictReader(handle))
    assert len(rows) == 1
    assert rows[0]["hcdr3_plddt_mean"] == ""
    assert rows[0]["gnina_rank1_CNNscore"] == ""


def test_default_reporting_keeps_original_ref2015_schema(tmp_path, monkeypatch):
    batch = _batch(tmp_path, monkeypatch, weights=(None,))
    reports.generate(batch)
    payload = _payload(batch)
    assert payload["columns"] == reports._grouped_columns()
    assert not LEGACY_KEYS & payload["rows"][0].keys()
    assert "hcdr3_plddt_mean" not in payload["rows"][0]
    assert "gnina_rank1_CNNscore" not in payload["rows"][0]


def test_mixed_scorefunctions_cannot_be_combined_in_one_grouped_report(tmp_path, monkeypatch):
    batch = _batch(tmp_path, monkeypatch, weights=("ligand", "ref2015"))
    with pytest.raises(ValueError, match="(?i)(mixed|weights|score|profile)"):
        reports.generate(batch)


def test_late_native_boltz_source_is_joined_without_restarting_analysis(tmp_path, monkeypatch):
    from hashlib import sha256

    batch = _batch(tmp_path, monkeypatch)
    affinity = tmp_path / "boltz_affinity.csv"
    fields = ("affinity_pred_value", "affinity_probability_binary", "affinity_pred_value1",
              "affinity_probability_binary1", "affinity_pred_value2", "affinity_probability_binary2")
    write_csv(affinity, [{"job_name": "case_0", "result_status": "success",
                          **{field: index / 10 for index, field in enumerate(fields)}}])
    (batch / "report_sources.json").write_text(json.dumps({
        "scores_csv": str(affinity), "scores_csv_sha256": sha256(affinity.read_bytes()).hexdigest(),
    }))
    reports.generate(batch, include_empty_metadata_sections=True)
    payload = _payload(batch)
    assert payload["rows"][0]["boltz2_affinity_pred_value"] == 0.0
    assert payload["rows"][0]["rosetta_fa_pair_REU"] == -0.4
    assert payload["rows"][0]["gnina_rank1_CNNscore"] is None
    assert payload["rows"][0]["boltz2_complex_plddt"] is None
    assert payload["metadata_availability"]["Boltz-2"] == "partial: 6/15 source fields supplied"
    assert payload["metadata_availability"]["GNINA"].startswith("unavailable:")
    assert len(payload["columns"]) == 110
    assert str(affinity.name) in payload["model_scores_note"]
    assert payload["columns"][-44:] == confidence_columns() + score_columns()


def test_late_source_hash_is_checked_before_reporting(tmp_path, monkeypatch):
    batch = _batch(tmp_path, monkeypatch)
    affinity = tmp_path / "boltz_affinity.csv"
    affinity.write_text("job_name\nchanged\n")
    (batch / "report_sources.json").write_text(json.dumps({
        "scores_csv": str(affinity), "scores_csv_sha256": "not-the-validated-hash",
    }))
    with pytest.raises(ValueError, match="source changed"):
        reports.generate(batch, include_empty_metadata_sections=True)


def test_late_cif_sources_are_forwarded_to_confidence_reporting(tmp_path, monkeypatch):
    from fv_energy_analysis import confidence

    batch = _batch(tmp_path, monkeypatch)
    originals = batch / "original_cifs"
    originals.mkdir()
    parent_cif = batch / "m7_parent.cif"
    parent_cif.write_text("data_parent\n")
    (batch / "report_sources.json").write_text(json.dumps({
        "cif_folder": str(originals), "parent_cif": str(parent_cif),
    }))
    calls = []

    def saved_confidence_only(actual_batch, actual_folder, actual_parent):
        calls.append((actual_batch, actual_folder, actual_parent))
        return actual_batch / "batch_summary_grouped_fingerprints.csv"

    monkeypatch.setattr(confidence, "add_confidence", saved_confidence_only)
    reports.generate(batch, include_empty_metadata_sections=True)
    assert calls == [(batch, originals, parent_cif)]
    assert _payload(batch)["metadata_availability"]["HCDR3 pLDDT"] == "available"
