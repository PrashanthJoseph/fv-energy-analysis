"""Selected reports preserve canonical data and never rerun scientific work."""
import json
from contextlib import closing
from pathlib import Path

import pytest

from fv_energy_analysis.fingerprint_reports import _grouped_columns
from fv_energy_analysis.report import write_csv, write_json
from fv_energy_analysis.report_options import export_selected_report, prepare_selected_report, project_payload
from fv_energy_analysis.report_options import validate_report_selection


@pytest.fixture
def canonical(tmp_path):
    batch = tmp_path / "batch"
    batch.mkdir()
    payload = {"title": "Saved analysis", "columns": _grouped_columns(weights="ligand"),
               "scoring_profile": "rosettaligand_legacy", "parent": {"candidate": "Parent reference"},
               "rows": [{"rank": 1, "candidate": "example.pdb", "mutations": "E99F",
                         "rosetta_interface_REU": -12.5, "plip_hbond_total": 0}]}
    write_json(batch / "batch_summary_grouped_fingerprints.json", payload)
    write_csv(batch / "batch_summary_grouped_fingerprints.csv", payload["rows"])
    write_json(batch / "batch_manifest.json", {"candidate_count": 1, "completed_count": 1, "failure_count": 0})
    return batch, payload


def test_section_aliases_and_column_union_keep_identity_and_source_order(canonical):
    _, payload = canonical
    projected = project_payload(payload, ["PLIP"], ["rosetta_interface_REU"])
    keys = [c["key"] for c in projected["columns"]]
    assert keys[:4] == ["rank", "mutations", "candidate", "rosetta_interface_REU"]
    assert "plip_hbond_total" in keys and "arpeggio_hbond" not in keys
    assert projected["rows"][0]["plip_hbond_total"] == 0
    assert projected["rows"][0]["rosetta_interface_REU"] == -12.5
    assert [c["key"] for c in project_payload(payload, ["Rosetta"])["columns"]][3] == "rosetta_interface_REU"
    assert len(payload["columns"]) > len(projected["columns"])


@pytest.mark.parametrize("kwargs", [{"sections": ["typo"]}, {"columns": ["fake_energy"]},
                                   {"sections": ["PLIP", "PLIP"]}])
def test_invalid_selection_is_explicit(canonical, kwargs):
    with pytest.raises(ValueError):
        project_payload(canonical[1], **kwargs)


def test_config_selection_is_checked_before_sources_or_analysis_exist():
    validate_report_selection({"weights": "ligand", "report_sections": ["rosetta", "plip", "plddt"]})
    validate_report_selection({"weights": "custom.wts", "weights_identity": {"legacy_pre_talaris": True},
                               "report_columns": ["rosetta_fa_pair_REU"]})
    with pytest.raises(ValueError, match="Unknown report section"):
        validate_report_selection({"report_sections": ["E3 discrimination"], "e3_discrimination": False})
    with pytest.raises(ValueError, match="Unknown report columns"):
        validate_report_selection({"weights": "ref2015", "report_columns": ["rosetta_fa_pair_REU"]})
    with pytest.raises(ValueError, match="Report formats"):
        validate_report_selection({"report_formats": ["pdf"]})


def test_csv_values_are_checked_even_if_its_manifest_hash_is_changed(canonical):
    from hashlib import sha256
    output = export_selected_report(canonical[0], sections=["PLIP"], formats=["csv"])
    path = output / "batch_summary_grouped_fingerprints.csv"
    path.write_text(path.read_text().replace(",0,", ",9,"))
    manifest_path = output / "selected_report_manifest.json"
    manifest = json.loads(manifest_path.read_text())
    manifest["csv_sha256"] = sha256(path.read_bytes()).hexdigest()
    write_json(manifest_path, manifest)
    with pytest.raises(ValueError, match="CSV values"):
        prepare_selected_report(output)


def test_csv_projection_keeps_canonical_bytes_and_binds_source(canonical):
    batch, _ = canonical
    before = {p.name: p.read_bytes() for p in batch.iterdir() if p.is_file()}
    output = export_selected_report(batch, sections=["Rosetta"], formats=["csv"])
    assert output == batch / "selected_report"
    payload, manifest, failures = prepare_selected_report(output)
    assert manifest["scientific_analysis_rerun"] is False and failures == []
    assert len(payload["rows"]) == 1
    assert {p.name: p.read_bytes() for p in batch.iterdir() if p.is_file()} == before
    assert not (output / "batch_summary_grouped_fingerprints.xlsx").exists()
    source = batch / "batch_summary_grouped_fingerprints.json"
    source.write_text(source.read_text() + "\n")
    with pytest.raises(ValueError, match="canonical source data changed"):
        prepare_selected_report(output)


def test_export_rejects_overwrite_of_canonical_or_different_view(canonical):
    batch, _ = canonical
    with pytest.raises(ValueError, match="separate"):
        export_selected_report(batch, output=batch, formats=["csv"])
    export_selected_report(batch, sections=["PLIP"], formats=["csv"])
    with pytest.raises(ValueError, match="different source or selection"):
        export_selected_report(batch, sections=["Rosetta"], formats=["csv"])


def test_selected_value_tampering_is_detected(canonical):
    output = export_selected_report(canonical[0], sections=["PLIP"], formats=["csv"])
    path = output / "batch_summary_grouped_fingerprints.json"
    saved = json.loads(path.read_text())
    saved["rows"][0]["plip_hbond_total"] = 7
    write_json(path, saved)
    with pytest.raises(ValueError, match="differs from its canonical"):
        prepare_selected_report(output)


def test_identity_only_xlsx_roundtrip(canonical):
    pytest.importorskip("openpyxl")
    from openpyxl import load_workbook
    from fv_energy_analysis.workbook_validation import verify
    output = export_selected_report(canonical[0], sections=[], formats=["csv", "xlsx"])
    result = verify(output, write_report=False)
    assert result["columns"] == 3 and result["filter_range"] == "A5:C6"
    with closing(load_workbook(output / "batch_summary_grouped_fingerprints.xlsx")) as book:
        assert book["Fingerprints"]["C6"].value == "example.pdb"


def test_skipped_fingerprints_are_blank_even_when_stale_data_exists(tmp_path):
    from fv_energy_analysis.fingerprint_counts import candidate_fingerprints, _raw_count_rows
    write_json(tmp_path / "workflow_manifest.json", {"ligand_resname": "E2L", "analysis_status": {
        "arpeggio": "not_run", "plip": "not_run", "discovery_studio": "not_run"}})
    for name in ("arpeggio_parent", "arpeggio_mutant"):
        (tmp_path / name).mkdir()
        (tmp_path / name / ("parent_arpeggio.json" if name.endswith("parent") else "mutant_arpeggio.json")).write_text("invalid stale JSON")
    summary = {"ligand_resname": "E2L", "plip_ligand_interactions": {
        role: {"interaction_counts": {"hydrogen_bonds": 9}, "total_interaction_count": 9}
        for role in ("parent", "mutant")}}
    rows = candidate_fingerprints("example.pdb", tmp_path, tmp_path / "missing_ds.csv", summary)
    assert all(row["mutant_total"] is None for row in rows if row["tool"] != "Rosetta")
    native = _raw_count_rows("example.pdb", tmp_path, summary, tmp_path / "missing_ds.csv")
    assert all(row["mutant_count"] is None for row in native if row["tool"] == "PLIP")
