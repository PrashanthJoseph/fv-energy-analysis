"""Replace failed candidates with compatible completed analyses by identity.

Only global ranks change for the existing completed cohort. New candidates
must match the exact failed inventory and the same parent, ligand parameters,
RosettaLigand protocol, detector build and E2/E3 definitions.
"""
from __future__ import annotations

import argparse
from collections import Counter
from copy import deepcopy
import csv
from datetime import datetime, timezone
import math
from pathlib import Path

from .architecture_batch import atom_snapshot
from .combine_batches import (
    _scientific_context, _typed, _validate_table, candidate_name, digest, read_json,
)
from .e3_reports import SECTION as E3_SECTION, append_e3_to_payload, flatten_e3_result
from .paths import native_path
from .report import write_json


TABLE = "batch_summary_grouped_fingerprints"
STATUS_COLUMNS = {"status", "error"}
_E3_IDENTITY = "PROQIPRRNZUXQM-ZXXIGWHRSA-N"


def _schema(columns):
    return [tuple(column[k] for k in ("key", "label", "section", "kind")) for column in columns]


def _write_table(folder, payload):
    write_json(folder / f"{TABLE}.json", payload)
    csv_path = folder / f"{TABLE}.csv"
    with csv_path.open("w", encoding="utf-8-sig", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=[c["key"] for c in payload["columns"]], extrasaction="ignore")
        writer.writeheader()
        writer.writerows(payload["rows"])
    return csv_path


def _csv_path(folder, manifest):
    recorded = native_path(manifest.get("csv", folder / f"{TABLE}.csv"))
    expected = manifest.get("csv_sha256")
    if recorded.is_file():
        if expected and digest(recorded) != expected:
            raise ValueError(f"Source CSV hash differs: {recorded}")
        return recorded
    matches = [p for p in folder.glob("*.csv") if expected and digest(p) == expected]
    if len(matches) != 1:
        raise ValueError(f"Cannot uniquely locate the saved source CSV: {folder}")
    return matches[0]


def _validate_combined(folder, payload, manifest):
    rows, columns = payload["rows"], payload["columns"]
    if manifest.get("status") not in {"complete", "complete_with_failures"}:
        raise ValueError("Base combined analysis has not completed")
    names = [candidate_name(r["candidate"]) for r in rows]
    if len(set(n.casefold() for n in names)) != len(names):
        raise ValueError("Duplicate candidate identities in base dataset")
    completed = sum(row.get("status") == "complete" for row in rows)
    failed = sum(row.get("status") == "failed" for row in rows)
    if (len(rows) != manifest["candidate_count"] or completed != manifest["completed_count"]
            or failed != manifest["failure_count"] or completed + failed != len(rows)):
        raise ValueError("Base candidate statuses/counts disagree with its manifest")
    with _csv_path(folder, manifest).open(encoding="utf-8-sig", newline="") as stream:
        reader = csv.DictReader(stream)
        actual = list(reader)
        if reader.fieldnames != [c["key"] for c in columns] or len(actual) != len(rows):
            raise ValueError("Base CSV/JSON schema or row count differs")
    for index, (row, record) in enumerate(zip(rows, actual), 1):
        if any(_typed(row.get(c["key"]), c) != _typed(record.get(c["key"]), c) for c in columns):
            raise ValueError(f"Base CSV/JSON value differs at row {index}")
    if payload.get("scoring_profile") != "rosettaligand_legacy":
        raise ValueError("Expected the legacy RosettaLigand base dataset")
    if not payload.get("e3_discrimination_protocol", {}).get("hydrogen_bond_geometry_assessed"):
        raise ValueError("Base must use the geometry-checked E3 protocol")
    for row in [payload["parent"], *(r for r in rows if r["status"] == "complete")]:
        if row.get("e3_status") != "complete":
            raise ValueError(f"Base completed row lacks its E3 assessment: {row.get('candidate')}")


def _compare_parent(base, replacement, columns):
    """Scientific values must match; absent imported annotations reuse baseline."""
    imported_sections = {"HCDR3 pLDDT", "Boltz-2", "GNINA"}
    for column in columns:
        key = column["key"]
        if key in {"rank", "candidate", "status", "error"}:
            continue
        old, new = _typed(base.get(key), column), _typed(replacement.get(key), column)
        if column["section"] in imported_sections and new is None:
            continue
        if old != new:
            raise ValueError(f"Replacement parent baseline differs: {key}: {old!r} != {new!r}")


def _geometry_context(result):
    protocol = result["protocol"]
    return {key: protocol.get(key) for key in (
        "method", "scope", "warning_overlap_A", "radii_A", "hbond_geometry_criteria",
        "geometry_combination_policy",
    )}


def _validate_e3(result, candidate, prepared, expected_context):
    if result.get("candidate") != candidate or result.get("status") != "complete":
        raise ValueError(f"E3 result candidate/status mismatch: {candidate}")
    values = flatten_e3_result(result)
    if result.get("schema_version") != 2 or _geometry_context(result) != expected_context:
        raise ValueError(f"E3 geometry protocol differs: {candidate}")
    chemistry = result["chemistry"]
    if (chemistry.get("inchikey") != _E3_IDENTITY or chemistry.get("formula") != "C18H24O3"
            or chemistry.get("formal_charge") != 0 or chemistry.get("ligand_hydrogen_count") != 24):
        raise ValueError(f"E3 identity/protonation differs: {candidate}")
    probe, source = native_path(result["probe_pdb"]), native_path(result["input_pdb"])
    if digest(probe) != result["probe_pdb_sha256"] or digest(source) != result["input_pdb_sha256"]:
        raise ValueError(f"E3 probe/input hash differs: {candidate}")
    if atom_snapshot(source) != atom_snapshot(prepared):
        raise ValueError(f"E3 input differs from the scored prepared complex: {candidate}")
    ligand = [line for line in probe.read_text(encoding="utf-8").splitlines()
              if line.startswith("HETATM") and line[17:20] == "E3L" and line[21:22] == "X"]
    elements = Counter(line[76:78].strip() for line in ligand)
    if elements != {"C": 18, "O": 3, "H": 24}:
        raise ValueError(f"E3 probe must contain 21 heavy atoms and 24 H: {candidate}")
    return values


def merge_rows(base_rows, replacement_rows, replacement_e3):
    """Pure identity merge, refusing any change to an already complete case."""
    old = {candidate_name(r["candidate"]): r for r in base_rows}
    if len(old) != len(base_rows) or len({n.casefold() for n in old}) != len(old):
        raise ValueError("Duplicate candidate identities in base dataset")
    new = {candidate_name(r["candidate"]): r for r in replacement_rows}
    if len(new) != len(replacement_rows) or len({n.casefold() for n in new}) != len(new):
        raise ValueError("Duplicate candidate identities in replacements")
    expected = {name for name, row in old.items() if row.get("status") == "failed"}
    if set(new) != expected:
        raise ValueError(f"Replacement inventory differs from failed candidates; missing={sorted(expected-set(new))}, unexpected={sorted(set(new)-expected)}")
    if set(replacement_e3) != expected:
        raise ValueError("Replacement E3 inventory differs from the exact failed candidate set")
    rows = []
    for name, previous in old.items():
        if name not in expected:
            rows.append(deepcopy(previous))
            continue
        row = deepcopy(new[name])
        if row.get("status") not in (None, "complete"):
            raise ValueError(f"Replacement analysis is not complete: {name}")
        row.update(candidate=name, status="complete", error=None)
        row.update(flatten_e3_result(replacement_e3[name]))
        if row["e3_status"] != "complete":
            raise ValueError(f"Replacement E3 analysis is not complete: {name}")
        rows.append(row)
    for row in rows:
        for key in ("rosetta_interface_REU", "rosetta_interface_delta_REU"):
            if isinstance(row.get(key), bool) or not math.isfinite(float(row.get(key))):
                raise ValueError(f"Missing/nonfinite interface energy: {row['candidate']}: {key}")
    rows.sort(key=lambda row: (float(row["rosetta_interface_delta_REU"]), row["candidate"]))
    for rank, row in enumerate(rows, 1):
        row["rank"] = rank
    current_by_name = {row["candidate"]: row for row in rows}
    for name, previous in old.items():
        if name in expected:
            continue
        current = current_by_name[name]
        if any(previous.get(key) != current.get(key) for key in previous if key != "rank"):
            raise ValueError(f"Existing completed metrics changed: {name}")
    return rows, expected


def replace_failed_candidates(base, replacement, output, replacement_e3, expected_count=None):
    base, replacement, output, replacement_e3 = [native_path(p).resolve() for p in
                                                (base, replacement, output, replacement_e3)]
    if any(output == p or p in output.parents for p in (base, replacement, replacement_e3)):
        raise ValueError("Use a fresh output folder outside all source folders")
    if output.exists() and any(output.iterdir()):
        raise FileExistsError(f"Output must be fresh: {output}")
    payload_path, manifest_path = base / f"{TABLE}.json", base / "combined_manifest.json"
    base_hash = digest(payload_path)
    replacement_hash = digest(replacement / f"{TABLE}.json")
    payload, manifest = read_json(payload_path), read_json(manifest_path)
    _validate_combined(base, payload, manifest)
    replacement_payload = read_json(replacement / f"{TABLE}.json")
    replacement_manifest = read_json(replacement / "batch_manifest.json")
    _validate_table(replacement, replacement_payload, replacement_manifest)
    scientific_columns = [c for c in payload["columns"] if c["section"] != E3_SECTION and c["key"] not in STATUS_COLUMNS]
    replacement_columns = [c for c in replacement_payload["columns"] if c["section"] != E3_SECTION and c["key"] not in STATUS_COLUMNS]
    if _schema(scientific_columns) != _schema(replacement_columns):
        raise ValueError("Replacement scientific column schema differs from the final dataset")
    context = _scientific_context(replacement, replacement_manifest)
    if context != manifest["scientific_protocol"]:
        raise ValueError("Replacement scientific protocol, parent or parameter hashes differ")
    _compare_parent(payload["parent"], replacement_payload["parent"], scientific_columns)
    e3_manifest = read_json(replacement_e3 / "e3_batch_manifest.json")
    if e3_manifest.get("status") != "complete" or e3_manifest.get("status_counts", {}).get("failed", 0):
        raise ValueError("Replacement E3 batch has not completed successfully")
    base_e3_manifest = read_json(native_path(manifest["e3_discrimination"]))
    expected_geometry = _geometry_context(base_e3_manifest["parent"])
    folders = {candidate_name(case["candidate"]): native_path(case["result_folder"])
               for case in replacement_manifest["candidates"] if case["status"] == "complete"}
    e3_results = {}
    for path in replacement_e3.glob("cases/*/e3_discrimination.json"):
        result = read_json(path)
        name = candidate_name(result["candidate"])
        if name in e3_results:
            raise ValueError(f"Duplicate replacement E3 result: {name}")
        if name not in folders:
            raise ValueError(f"Unexpected replacement E3 result: {name}")
        summary = read_json(folders[name] / "final/summary.json")
        prepared = native_path(summary["mutant"])
        if not prepared.is_file():
            prepared = folders[name] / "prepared" / prepared.name
        _validate_e3(result, name, prepared, expected_geometry)
        e3_results[name] = result
    parent_e3_path = replacement_e3 / "parent/e3_discrimination.json"
    if parent_e3_path.is_file():
        parent_e3 = read_json(parent_e3_path)
        first_folder = next(iter(folders.values()))
        parent_prepared = native_path(read_json(first_folder / "final/summary.json")["parent"])
        _validate_e3(parent_e3, "parent_m7", parent_prepared, expected_geometry)
        if flatten_e3_result(parent_e3) != {key: payload["parent"].get(key) for key in flatten_e3_result(parent_e3)}:
            raise ValueError("Replacement E3 parent baseline differs")
    rows, replaced = merge_rows(payload["rows"], replacement_payload["rows"], e3_results)
    if expected_count is not None and len(rows) != expected_count:
        raise ValueError(f"Expected {expected_count} candidates, found {len(rows)}")
    with (base / "candidate_sources.csv").open(encoding="utf-8-sig", newline="") as stream:
        origins = list(csv.DictReader(stream))
    by_name = {r["candidate"]: r for r in origins}
    if len(by_name) != len(origins) or set(by_name) != {r["candidate"] for r in rows}:
        raise ValueError("Base candidate provenance does not cover the exact inventory")
    new_rank = {r["candidate"]: r.get("rank") for r in replacement_payload["rows"]}
    for name in replaced:
        by_name[name] = {"candidate": name, "source": str(replacement), "source_rank": new_rank[name]}
    grouped_origins = Counter(by_name[r["candidate"]]["source"] for r in rows)
    previous_sources = {record["source"]: record for record in manifest.get("sources", [])}
    source_records = []
    for source_name, count in grouped_origins.items():
        record = deepcopy(previous_sources.get(source_name, {}))
        if record:
            record["original_source_candidate_count"] = record["candidate_count"]
            record["selection"] = "retained completed rows; recovered failed candidates use their replacement cohort"
        else:
            record.update(manifest_sha256=digest(replacement / "batch_manifest.json"),
                          payload_sha256=replacement_hash,
                          csv_sha256=digest(replacement / f"{TABLE}.csv"))
        record.update(source=source_name, candidate_count=count, completed_count=count, failure_count=0)
        source_records.append(record)
    if digest(payload_path) != base_hash or digest(replacement / f"{TABLE}.json") != replacement_hash:
        raise ValueError("A source dataset changed during replacement validation")
    output.mkdir(parents=True, exist_ok=True)
    # Bind the final E3 export to a complete pre-E3 snapshot. The existing
    # verifier can then check all non-E3 values without exceptions for recovery.
    snapshot_folder = output / "pre_e3_combined"
    snapshot_folder.mkdir()
    snapshot = deepcopy(payload)
    snapshot["columns"] = [c for c in payload["columns"] if c["section"] != E3_SECTION]
    e3_keys = {c["key"] for c in payload["columns"] if c["section"] == E3_SECTION}
    snapshot["rows"] = [{k: v for k, v in r.items() if k not in e3_keys} for r in rows]
    snapshot["parent"] = {k: v for k, v in payload["parent"].items() if k not in e3_keys}
    snapshot.pop("e3_discrimination_protocol", None)
    snapshot["combined_sources"] = source_records
    snapshot.update(title="Complete m7 RosettaLigand analysis", summary_subtitle=f"All {len(rows)} candidates complete; global rank by mutant-minus-parent interface energy.")
    snapshot_csv = _write_table(snapshot_folder, snapshot)
    combined = deepcopy(payload)
    combined.update(rows=rows, title="Complete m7 RosettaLigand analysis",
                    summary_subtitle=f"All {len(rows)} candidates complete; recovered {len(replaced)} prior exclusions. E3 overlaps use explicit H-bond geometry.",
                    note="All existing completed metrics are preserved. Previously rejected candidates were replaced by corrected analyses; global ranks were recalculated.")
    combined["combined_sources"] = source_records
    combined["e3_discrimination_protocol"].update({
        "source_combined_batch": str(snapshot_folder), "source_payload_sha256": digest(snapshot_folder / f"{TABLE}.json"),
        "replacement_e3_batch": str(replacement_e3), "previous_e3_batch": str(base),
        "replacement_count": len(replaced)})
    final_csv = _write_table(output, combined)
    for folder in (output, snapshot_folder):
        with (folder / "candidate_sources.csv").open("w", encoding="utf-8-sig", newline="") as stream:
            writer = csv.DictWriter(stream, fieldnames=["candidate", "source", "source_rank"])
            writer.writeheader()
            writer.writerows(by_name[r["candidate"]] for r in rows)
    recovery = {"base": str(base), "base_payload_sha256": base_hash,
                "replacement": str(replacement), "replacement_payload_sha256": replacement_hash,
                "replacement_e3": str(replacement_e3), "replacement_e3_manifest_sha256": digest(replacement_e3 / "e3_batch_manifest.json"),
                "replaced_candidates": sorted(replaced), "replaced_count": len(replaced),
                "existing_completed_metrics_preserved": True, "only_existing_metric_changed": "rank"}
    common_manifest = deepcopy(manifest)
    common_manifest.update(status="complete", completed_utc=datetime.now(timezone.utc).isoformat(),
                           candidate_count=len(rows), completed_count=len(rows), failure_count=0,
                           expected_count=len(rows), metrics_changed=["rank", "recovered failed candidate analyses"],
                           replacement=recovery, scoring_rerun=False, detectors_rerun=[], sources=source_records)
    for folder, table, csv_path in ((snapshot_folder, snapshot, snapshot_csv), (output, combined, final_csv)):
        saved = deepcopy(common_manifest)
        saved.update(columns=len(table["columns"]), csv=str(csv_path), csv_sha256=digest(csv_path),
                     interaction_fingerprints={"grouped_summary_xlsx": str(folder / f"{TABLE}.xlsx")})
        if folder == snapshot_folder:
            saved.pop("e3_discrimination", None)
        else:
            saved["e3_discrimination"] = str(output / "e3_final_manifest.json")
        write_json(folder / "combined_manifest.json", saved)
    write_json(output / "e3_final_manifest.json", {
        "status": "complete", "candidate_count": len(rows), "completed_count": len(rows),
        "parent": base_e3_manifest["parent"], "existing_e3_batch": str(base),
        "replacement_e3_batch": str(replacement_e3), "replacement_count": len(replaced),
        "protocol": expected_geometry, "sources_recomputed": False})
    write_json(output / "replacement_verification.json", recovery)
    (output / "final_summary.md").write_text(
        f"# Complete m7 analysis\n\nAll {len(rows):,} candidates are complete, including {len(replaced)} corrected prior exclusions. "
        "The previously completed analyses, imported model/confidence scores and E3 results were reused unchanged. "
        "Only global ranks were recalculated. New candidates use the same parent, one OptH pass, RosettaLigand weights, "
        "native DS E2 recognition and geometry-checked E3 steric protocol.\n", encoding="utf-8")
    return final_csv


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--base", type=Path, required=True)
    parser.add_argument("--replacement", type=Path, required=True)
    parser.add_argument("--replacement-e3", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--expected-count", type=int)
    args = parser.parse_args()
    print(replace_failed_candidates(args.base, args.replacement, args.output, args.replacement_e3, args.expected_count), flush=True)


if __name__ == "__main__":
    main()
