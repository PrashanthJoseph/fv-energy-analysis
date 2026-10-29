"""Postprocessing must use the exact recorded structure and preserve failures."""
from hashlib import sha256
import csv
import json

from fv_energy_analysis.e3_discrimination import resolve_saved_complexes, write_flat_csv


def test_saved_refresh_uses_recorded_original_pdb_and_keeps_rejections(tmp_path):
    combined, source = tmp_path / "combined", tmp_path / "refresh"
    combined.mkdir()
    case = source / "valid"
    case.mkdir(parents=True)
    original = tmp_path / "original_prepared.pdb"
    original.write_text("prepared coordinate input")
    provenance = {"pdb": str(original), "pdb_sha256": sha256(original.read_bytes()).hexdigest()}
    result = {"e2_architecture_status": "complete", "e2_provenance": provenance}
    (case / "estrogen_binding_architecture.json").write_text(json.dumps({"result": result}))
    case2 = source / "valid2"
    case2.mkdir()
    (case2 / "estrogen_binding_architecture.json").write_text(json.dumps({"mutant": result, "parent": result}))
    payload = {"rows": [{"candidate": "valid.pdb", "status": "complete"},
                        {"candidate": "valid2.pdb", "status": "complete"},
                        {"candidate": "invalid.pdb", "status": "failed", "error": "wrong stereochemistry"}]}
    (combined / "batch_summary_grouped_fingerprints.json").write_text(json.dumps(payload))
    with (combined / "candidate_sources.csv").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=["candidate", "source"])
        writer.writeheader()
        writer.writerows({"candidate": r["candidate"], "source": str(source)} for r in payload["rows"])
    _, tasks, excluded, parent = resolve_saved_complexes(combined)
    assert len(tasks) == 2 and all(t["pdb"] == str(original) for t in tasks)
    assert parent["pdb_sha256"] == provenance["pdb_sha256"]
    assert excluded == [{"candidate": "invalid.pdb", "status": "excluded_input_rejection",
                         "error": "wrong stereochemistry", "analysis_scope": "fixed_pose_steric_only"}]


def test_flat_export_follows_schema_not_dictionary_insertion_order(tmp_path):
    path = tmp_path / "out.csv"
    write_flat_csv(path, {"columns": [{"key": "candidate"}, {"key": "score"}],
                          "rows": [{"score": None, "candidate": "rejected", "extra": 4}]})
    with path.open(encoding="utf-8-sig") as stream:
        reader = csv.DictReader(stream)
        assert reader.fieldnames == ["candidate", "score"]
        assert list(reader) == [{"candidate": "rejected", "score": ""}]
