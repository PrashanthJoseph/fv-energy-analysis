"""E3 reports never infer missing values or change existing screening/ranks."""
import copy
import csv
import json

import pytest

from fv_energy_analysis.e3_reports import (
    SECTION, append_e3_report, append_e3_to_payload, e3_columns,
    flatten_e3_result, load_saved_e3, render_e3_markdown,
)


def saved(candidate="mutant.pdb"):
    return {"schema_version": 1, "status": "complete", "candidate": candidate,
            "summary": {"e3_analysis_scope": "fixed_pocket_steric_only",
                        "e3_added_o_closest_distance_A": 2.5,
                        "e3_added_o_closest_residue": "H:112:ASP",
                        "e3_added_o_deepest_overlap_A": .72,
                        "e3_added_o_overlap_residue_count": 2,
                        "e3_added_o_non_hbond_overlap_residue_count": 1,
                        "e3_added_o_hbond_supported_overlap_residue_count": 1,
                        "e3_added_o_unassessed_overlap_residue_count": 0,
                        "e3_added_o_overlap_residues": "H:112:ASP;L:35:SER",
                        "e3_steric_assessment": "steric_overlap_detected",
                        "e3_steric_warning": True}}


def payload():
    return {"columns": [{"key": "rank", "section": "Identity"},
                         {"key": "e2_score", "section": "E2 architecture"},
                         {"key": "plddt", "section": "HCDR3 pLDDT"}],
            "rows": [{"candidate": "mutant.pdb", "rank": 9, "e2_score": 4,
                      "rosetta_interface_delta_REU": -2.3},
                     {"candidate": "unassessed.pdb", "rank": 10, "e2_score": 3},
                     {"candidate": "rejected.pdb", "status": "failed", "rank": None}],
            "parent": {"candidate": "Parent reference", "e2_score": 4}}


def test_saved_measurements_are_typed_and_missing_are_not_zero():
    values = flatten_e3_result(saved())
    assert values["e3_steric_warning"] is True
    assert values["e3_added_o_overlap_residue_count"] == 2
    assert values["e3_steric_assessment"] == "steric_overlap_detected"
    assert values["e3_oh_h_best_sampled_protein_overlap_A"] is None
    unassessed = flatten_e3_result(None)
    assert unassessed["e3_status"] == "not_assessed"
    assert all(v is None for k, v in unassessed.items() if k != "e3_status")
    failed = flatten_e3_result({"status": "failed", "summary": saved()["summary"]})
    assert failed["e3_steric_warning"] is None
    assert failed["e3_added_o_closest_distance_A"] is None


def test_optional_sampled_hydrogen_metric_is_separate_from_primary_warning():
    result = saved()
    result["summary"]["e3_oh_h_best_sampled_protein_overlap_A"] = .17
    values = flatten_e3_result(result)
    assert values["e3_oh_h_best_sampled_protein_overlap_A"] == .17
    assert values["e3_steric_warning"] is True
    result["summary"]["e3_oh_h_best_sampled_protein_overlap_A"] = -1
    with pytest.raises(ValueError, match="OH-H overlap"):
        flatten_e3_result(result)


def test_polar_overlap_requires_geometry_and_unknown_warning_stays_null():
    result = saved()
    summary = result["summary"]
    summary.update(e3_added_o_overlap_residue_count=1,
                   e3_added_o_non_hbond_overlap_residue_count=0,
                   e3_added_o_hbond_supported_overlap_residue_count=0,
                   e3_added_o_unassessed_overlap_residue_count=1,
                   e3_steric_assessment="unassessed_polar_geometry",
                   e3_steric_warning=None)
    values = flatten_e3_result(result)
    assert values["e3_steric_warning"] is None
    assert values["e3_steric_assessment"] == "unassessed_polar_geometry"
    summary["e3_steric_warning"] = False
    with pytest.raises(ValueError, match="three-state warning"):
        flatten_e3_result(result)
    summary.update(e3_added_o_unassessed_overlap_residue_count=0,
                   e3_added_o_hbond_supported_overlap_residue_count=1,
                   e3_steric_assessment="hbond_supported_overlaps_only",
                   e3_steric_warning=False)
    assert flatten_e3_result(result)["e3_steric_warning"] is False


def test_unclassified_overlap_cannot_be_silently_reported_as_no_clash():
    result = saved()
    result["summary"].update(e3_added_o_non_hbond_overlap_residue_count=0,
                              e3_added_o_hbond_supported_overlap_residue_count=0,
                              e3_added_o_unassessed_overlap_residue_count=0,
                              e3_steric_assessment="no_large_overlap_in_fixed_pose",
                              e3_steric_warning=False)
    with pytest.raises(ValueError, match="three-state warning"):
        flatten_e3_result(result)


def test_payload_joins_by_identity_preserves_prior_metrics_and_is_idempotent():
    source = payload()
    old = copy.deepcopy(source)
    result = append_e3_to_payload(source, {"mutant.pdb": saved()})
    assert source == old
    assert result["rows"][0]["rank"] == 9
    assert result["rows"][0]["e2_score"] == 4
    assert result["rows"][0]["rosetta_interface_delta_REU"] == -2.3
    assert result["rows"][1]["e3_status"] == "not_assessed"
    assert result["rows"][2]["e3_status"] == "excluded_input_rejection"
    assert result["parent"]["e3_status"] == "not_assessed"
    sections = [c["section"] for c in result["columns"]]
    assert sections[1] == "E2 architecture"
    assert sections[2] == SECTION
    assert sections[-1] == "HCDR3 pLDDT"
    assert append_e3_to_payload(result, {"mutant.pdb": saved()}) == result
    assert not any("energy" in c["key"].lower() or "REU" in c["key"] for c in e3_columns())


def test_unknown_or_mismatched_results_fail_instead_of_wrong_join():
    with pytest.raises(ValueError, match="not present"):
        append_e3_to_payload(payload(), {"wrong.pdb": saved("wrong.pdb")})
    with pytest.raises(ValueError, match="Mismatched"):
        append_e3_to_payload(payload(), {"mutant.pdb": saved("wrong.pdb")})


def test_upgrading_the_e3_section_removes_obsolete_e3_columns_only():
    source = payload()
    old_key = "e3_added_o_polar_close_residue_count"
    source["columns"].append({"key": old_key, "section": SECTION})
    source["rows"][0][old_key] = 5
    upgraded = append_e3_to_payload(source, {"mutant.pdb": saved()})
    assert old_key not in upgraded["rows"][0]
    assert source["rows"][0][old_key] == 5
    assert upgraded["rows"][0]["rosetta_interface_delta_REU"] == -2.3


@pytest.mark.parametrize("key,value", [("e3_added_o_closest_distance_A", float("nan")),
                                     ("e3_added_o_overlap_residue_count", 1.2),
                                     ("e3_steric_warning", "False")])
def test_invalid_measurements_are_not_silently_exported(key, value):
    result = saved()
    result["summary"][key] = value
    with pytest.raises(ValueError):
        flatten_e3_result(result)


def test_loader_is_reporting_only_and_checks_candidate(tmp_path):
    assert load_saved_e3(tmp_path / "absent.json") is None
    path = tmp_path / "e3_discrimination.json"
    path.write_text(json.dumps(saved()))
    assert load_saved_e3(path, candidate="mutant.pdb")["status"] == "complete"
    with pytest.raises(ValueError, match="different candidate"):
        load_saved_e3(path, candidate="other.pdb")


def test_markdown_preserves_e2_and_replaces_only_its_own_block(tmp_path):
    path = tmp_path / "automatic_summary.md"
    original = "# Candidate\n\n<!-- e2-architecture:begin -->\nE2 result\n<!-- e2-architecture:end -->\n"
    path.write_text(original)
    append_e3_report(path, saved())
    append_e3_report(path, saved())
    text = path.read_text()
    assert original.rstrip() in text
    assert text.count("<!-- e3-discrimination:begin -->") == 1
    assert "0.4 A" in text
    assert "not proof" in text
    assert "Not assessed" in render_e3_markdown(None)


def test_geometry_table_exposes_real_passes_failures_and_unassessed_contacts():
    result = saved()
    result["contacts"] = [
        {"protein_atom": "H:112:ASP:OD2", "overlap_at_or_above_warning": True,
         "hbond_geometry": {"status": "supported", "reason": "donor geometry passes",
                            "trials": [{"donor": "X:1:E3L:O3", "hydrogen": "X:1:E3L:H25",
                                        "acceptor": "H:112:ASP:OD2", "donor_acceptor_distance_A": 2.8,
                                        "hydrogen_acceptor_distance_A": 1.8, "angle_DHA_deg": 169.5,
                                        "passes_geometry": True}]}},
        {"protein_atom": "L:35:SER:OG", "overlap_at_or_above_warning": True,
         "hbond_geometry": {"status": "not_supported", "reason": "donor angle fails",
                            "trials": [{"donor": "X:1:E3L:O3", "hydrogen": "X:1:E3L:H25",
                                        "acceptor": "L:35:SER:OG", "donor_acceptor_distance_A": 2.7,
                                        "hydrogen_acceptor_distance_A": 2.4, "angle_DHA_deg": 70,
                                        "passes_geometry": False}]}},
        {"protein_atom": "H:99:HIS:ND1", "overlap_at_or_above_warning": True,
         "hbond_geometry": {"status": "unassessed", "reason": "hydrogen typing missing", "trials": []}},
        {"protein_atom": "H:1:GLY:O", "overlap_at_or_above_warning": False,
         "hbond_geometry": {"status": "unassessed", "reason": "not an overlap", "trials": []}},
    ]
    markdown = render_e3_markdown(result)
    assert "Geometry of overlapping contacts" in markdown
    assert "2.800" in markdown and "1.800" in markdown and "169.5" in markdown
    assert "donor angle fails" in markdown and "70.0" in markdown
    assert "hydrogen typing missing" in markdown
    assert "H:1:GLY:O" not in markdown


def test_standard_generator_requires_explicit_opt_in(monkeypatch, tmp_path):
    from fv_energy_analysis import fingerprint_reports as reports
    batch = tmp_path / "batch"
    case = batch / "mutant"
    case.mkdir(parents=True)
    (case / "automatic_summary.json").write_text("{}")
    (case / "automatic_summary.md").write_text("# Candidate\n")
    (case / "e3_discrimination.json").write_text(json.dumps(saved()))
    with (batch / "batch_summary.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["candidate", "status", "result_folder"])
        writer.writeheader()
        writer.writerow({"candidate": "mutant.pdb", "status": "complete", "result_folder": str(case)})
    monkeypatch.setattr(reports, "candidate_fingerprints", lambda *args: [])
    monkeypatch.setattr(reports, "_raw_count_rows", lambda *args: [])
    monkeypatch.setattr(reports, "_saved_architecture", lambda *args: {"e2_architecture_status": "complete"})
    reports.generate(batch)
    grouped_path = batch / "batch_summary_grouped_fingerprints.json"
    first = json.loads(grouped_path.read_text())
    assert not any(c["section"] == SECTION for c in first["columns"])
    assert "e3_discrimination" not in json.loads((case / "automatic_summary.json").read_text())
    reports.generate(batch, include_e3_discrimination=True)
    second = json.loads(grouped_path.read_text())
    assert any(c["section"] == SECTION for c in second["columns"])
    assert second["rows"][0]["e3_steric_warning"] is True
    assert "<!-- e3-discrimination:begin -->" in (case / "automatic_summary.md").read_text()
    assert second["e3_discrimination_protocol"]["energies_assessed"] is False
