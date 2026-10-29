"""Append-only, human-readable issue records; never bypass failed validation."""
from __future__ import annotations

import csv
from datetime import datetime, timezone
from pathlib import Path
import re

FIELDS = ("recorded_utc", "candidate", "stage", "category", "reason", "workflow_log", "result_folder")


def record_issue(batch: Path, result: dict) -> Path:
    folder = Path(str(result["result_folder"]))
    log_path = folder / "workflow.log"
    error = str(result.get("error") or "Analysis did not complete")
    reason = error
    if log_path.is_file() and not error.startswith(("Worker crashed", "Ligand identity mismatch")):
        with log_path.open("rb") as stream:
            stream.seek(max(0, log_path.stat().st_size - 32768))
            tail = stream.read().decode("utf-8", errors="replace")
        exceptions = re.findall(r"^(?:\w*Error|\w*Exception): .+$", tail, flags=re.MULTILINE)
        if exceptions:
            reason = exceptions[-1]
    stage = next((name for marker, name in (
        ("prepare_ligand_complexes.py", "ligand_preparation"),
        ("pdbe-arpeggio", "arpeggio"), ("plip", "plip"),
        ("Worker crashed", "worker"), ("aggregate", "aggregation"),
    ) if marker in error), "analysis")
    category = "ligand_validation" if "ligand" in reason.lower() and stage == "ligand_preparation" else "analysis_failure"
    path = batch / "batch_issues.csv"
    needs_header = not path.exists() or path.stat().st_size == 0
    with path.open("a", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=FIELDS)
        if needs_header:
            writer.writeheader()
        writer.writerow(dict(zip(FIELDS, (
            datetime.now(timezone.utc).isoformat(), str(result["candidate"]),
            stage, category, reason, str(log_path), str(folder),
        ))))
    return path
