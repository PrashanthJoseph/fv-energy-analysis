import csv
from pathlib import Path

from fv_energy_analysis.issues import record_issue


def test_issue_log_preserves_actual_ligand_rejection_and_appends(tmp_path):
    folder = tmp_path / "case"
    folder.mkdir()
    (folder / "workflow.log").write_text("Traceback...\nValueError: Ligand 3D stereochemistry differs from reference\n")
    result = {"candidate": "mutant.pdb", "result_folder": str(folder),
              "error": "Command failed: prepare_ligand_complexes.py ..."}
    path = record_issue(tmp_path, result)
    record_issue(tmp_path, result)
    with path.open(newline="") as stream:
        rows = list(csv.DictReader(stream))
    assert len(rows) == 2
    assert rows[0]["category"] == "ligand_validation"
    assert rows[0]["stage"] == "ligand_preparation"
    assert rows[0]["reason"] == "ValueError: Ligand 3D stereochemistry differs from reference"
    assert Path(rows[0]["workflow_log"]) == folder / "workflow.log"


def test_crashed_worker_does_not_inherit_an_old_validation_error(tmp_path):
    (tmp_path / "workflow.log").write_text("ValueError: previous attempt\n")
    path = record_issue(tmp_path, {"candidate": "m", "result_folder": str(tmp_path), "error": "Worker crashed: broken pool"})
    with path.open(newline="") as stream:
        row = next(csv.DictReader(stream))
    assert row["reason"] == "Worker crashed: broken pool" and row["stage"] == "worker"


def test_targeted_diagnostic_is_not_replaced_by_old_generic_error(tmp_path):
    (tmp_path / "workflow.log").write_text("ValueError: old generic error\n")
    path = record_issue(tmp_path, {"candidate": "m", "result_folder": str(tmp_path),
                                  "error": "Ligand identity mismatch: C18 expected S observed R"})
    with path.open(newline="") as stream:
        row = next(csv.DictReader(stream))
    assert row["reason"] == "Ligand identity mismatch: C18 expected S observed R"
