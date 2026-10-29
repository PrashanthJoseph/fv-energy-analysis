"""Merge saved geometry/model metrics with a completed RosettaLigand rescore.

This is reporting only: it never scores structures or reruns fingerprint tools.
The XLSX renderer consumes the resulting grouped JSON payload.
"""
from __future__ import annotations

import argparse
import copy
import csv
import json
import math
from pathlib import Path
import shutil

PROFILE = "rosettaligand_legacy"
NETWORK_KEYS = (
    "network_added_pair_type_events", "network_removed_pair_type_events",
    "network_changed_unique_pairs_energy_delta_REU", "comparison_quality_warning",
)
SCORE_FIELDS = {
    "rosetta_interface_REU": "interface_REU",
    "rosetta_hbond_total": "ligand_hbonds",
    "rosetta_fa_atr_REU": "fa_atr_REU",
    "rosetta_fa_rep_REU": "fa_rep_REU",
    "rosetta_fa_elec_REU": "fa_elec_REU",
    "rosetta_solvation_REU": "fa_sol_REU",
    "rosetta_fa_pair_REU": "fa_pair_REU",
    "rosetta_hbond_bb_sc_REU": "hbond_bb_sc_REU",
    "rosetta_hbond_sc_REU": "hbond_sc_REU",
}


def report_columns(columns: list[dict]) -> list[dict]:
    """Preserve the layout and expose both additional active legacy terms."""
    result = []
    for original in columns:
        column = copy.deepcopy(original)
        if column["section"] == "Rosetta":
            column["section"] = "RosettaLigand"
        if column["key"] == "rosetta_hbond_sc_REU":
            for key, label in (
                ("rosetta_fa_pair_REU", "fa_pair (REU)"),
                ("rosetta_hbond_bb_sc_REU", "hbond_bb_sc (REU)"),
            ):
                extra = copy.deepcopy(column)
                extra.update(key=key, label=label)
                result.append(extra)
        result.append(column)
    return result


def _records(path: Path) -> list[dict]:
    with path.open(encoding="utf-8-sig", newline="") as stream:
        return list(csv.DictReader(stream))


def _number(value: str, key: str) -> float | int:
    number = float(value)
    if not math.isfinite(number):
        raise ValueError(f"Non-finite rescore value: {key}")
    if key == "rosetta_hbond_total":
        if not number.is_integer():
            raise ValueError("H-bond count must be an integer")
        return int(number)
    return number


def _score_values(record: dict, state: str) -> dict:
    values = {key: _number(record[f"{state}_{field}"], key)
              for key, field in SCORE_FIELDS.items()}
    values["rosetta_interface_delta_REU"] = (
        0.0 if state == "parent" else _number(record["interface_delta_REU"], "delta")
    )
    # The displayed seven interface terms cover all non-zero legacy terms here.
    total = sum(values[key] for key in SCORE_FIELDS
                if key not in ("rosetta_interface_REU", "rosetta_hbond_total"))
    if not math.isclose(total, values["rosetta_interface_REU"], abs_tol=1e-5):
        raise ValueError(f"Displayed interface decomposition does not sum: {total}")
    return values


def prepare_report(source: Path, rescore: Path) -> dict:
    """Create a new grouped CSV/JSON and rejection metadata, without XLSX I/O."""
    source, rescore = source.resolve(), rescore.resolve()
    if source == rescore:
        raise ValueError("The rescore report must not overwrite the source batch")
    source_payload = json.loads((source / "batch_summary_grouped_fingerprints.json").read_text(encoding="utf-8-sig"))
    score_manifest = json.loads((rescore / "rescore_manifest.json").read_text())
    if score_manifest["status"] != "complete" or score_manifest["failure_count"]:
        raise ValueError("A complete, failure-free rescore is required")
    scores = _records(rescore / "rosettaligand_rescore.csv")
    old_rows = {row["candidate"]: row for row in source_payload["rows"]}
    if len(old_rows) != len(source_payload["rows"]):
        raise ValueError("Duplicate source candidate")
    names = [row["candidate"] for row in scores]
    if len(names) != len(set(names)) or set(names) != set(old_rows):
        raise ValueError("Rescored candidates must match the source workbook exactly")
    networks = json.loads((rescore / "rescore_network_summary.json").read_text())
    if set(networks) != set(names):
        raise ValueError("Network summaries must match the rescored candidates exactly")
    parent = _score_values(scores[0], "parent")
    new_rows = []
    changed_keys = {"rank", *SCORE_FIELDS, "rosetta_interface_delta_REU", *NETWORK_KEYS}
    for rank, score in enumerate(scores, 1):
        if score["weights"] != "ligand" or int(score["rank_by_interface_delta"]) != rank:
            raise ValueError("Unexpected scoring weights or rank order")
        if _score_values(score, "parent") != parent:
            raise ValueError("Parent score differs between candidates")
        if networks[score["candidate"]].get("scorefunction") != "ligand":
            raise ValueError("Network energies do not use the ligand score function")
        row = copy.deepcopy(old_rows[score["candidate"]])
        if row["mutations"] != score["substitutions"]:
            raise ValueError(f"Substitution mismatch: {score['candidate']}")
        row.update(_score_values(score, "mutant"))
        row["rank"] = rank
        row.update({key: networks[score["candidate"]][key] for key in NETWORK_KEYS})
        if not isinstance(row["comparison_quality_warning"], bool):
            raise ValueError("Quality warning must be a boolean")
        if any(row[key] != value for key, value in old_rows[score["candidate"]].items()
               if key not in changed_keys):
            raise AssertionError("A geometry or model-confidence metric changed")
        expected_delta = row["rosetta_interface_REU"] - parent["rosetta_interface_REU"]
        if not math.isclose(expected_delta, row["rosetta_interface_delta_REU"], abs_tol=1e-10):
            raise ValueError("Invalid mutant-minus-parent interface delta")
        new_rows.append(row)
    if any(a["rosetta_interface_delta_REU"] > b["rosetta_interface_delta_REU"]
           for a, b in zip(new_rows, new_rows[1:])):
        raise ValueError("Candidates are not ordered by interface delta")

    payload = copy.deepcopy(source_payload)
    payload.update(
        title="m7 mutants: RosettaLigand rescore",
        scoring_profile=PROFILE,
        columns=report_columns(source_payload["columns"]),
        rows=new_rows,
        source_note=("Scoring-only recalculation with historical RosettaLigand ligand.wts and pre-talaris energy options. "
                     "Coordinates, ligand parameters and charges are unchanged. Rosetta H-bonds and mutation-network energies "
                     "were rebuilt with this score function. Arpeggio, PLIP, Discovery Studio, structural and model-confidence "
                     "metrics are reused unchanged from the corrected batch."),
        note="Ranked by RosettaLigand mutant-minus-parent interface energy. Lower is more favorable. REU are not experimental binding free energies.",
        summary_subtitle=(f"{len(scores):,} rescored. {score_manifest['excluded_source_failures']:,} excluded by source ligand validation. "
                          "Mutant values are shown; the pinned row is parent m7. Rank uses RosettaLigand interface delta."),
    )
    payload["parent"].update(parent)
    payload["scoring"] = {"weights": "ligand", "legacy_pre_talaris": True,
                          "source_batch": str(source), "manifest": str(rescore / "rescore_manifest.json"),
                          "extra_interface_terms": ["fa_pair", "hbond_bb_sc"]}
    (rescore / "batch_summary_grouped_fingerprints.json").write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    with (rescore / "batch_summary_grouped_fingerprints.csv").open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=[column["key"] for column in payload["columns"]], extrasaction="ignore")
        writer.writeheader()
        writer.writerows(new_rows)

    # Keep rejected source structures visible without presenting their old scores.
    status_rows = _records(source / "batch_summary.csv")
    excluded = [row for row in status_rows if row["status"] != "complete"]
    if len(status_rows) != len(new_rows) + len(excluded):
        raise ValueError("Source batch contains unaccounted candidates")
    if len(excluded) != score_manifest["excluded_source_failures"]:
        raise ValueError("Source rejection count differs from the rescore manifest")
    with (rescore / "batch_summary.csv").open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=["candidate", "candidate_path", "status", "error"])
        writer.writeheader()
        for row in status_rows:
            writer.writerow({"candidate": row["candidate"], "candidate_path": row["candidate_path"],
                             "status": "complete" if row["status"] == "complete" else "excluded_source_validation",
                             "error": row.get("error", "")})
    shutil.copyfile(source / "ligand_validation_diagnostics.csv", rescore / "ligand_validation_diagnostics.csv")
    manifest = {"status": "complete_with_source_exclusions", "candidate_count": len(status_rows),
                "completed_count": len(new_rows), "failure_count": len(excluded),
                "source_excluded_count": len(excluded), "rescore_failure_count": 0,
                "scoring_profile": PROFILE, "source_batch": str(source), "output": str(rescore)}
    (rescore / "batch_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    return {"rows": len(new_rows), "columns": len(payload["columns"]),
            "source_exclusions": len(excluded), "profile": PROFILE}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("source", type=Path)
    parser.add_argument("rescore", type=Path)
    args = parser.parse_args()
    print(json.dumps(prepare_report(args.source, args.rescore)))


if __name__ == "__main__":
    main()
