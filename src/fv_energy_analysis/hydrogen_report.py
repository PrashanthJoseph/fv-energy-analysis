"""Merge a completed H-only rescore with the previous grouped report.

This module performs reporting only. It never initializes Rosetta or reruns
fingerprint tools. Reused hydrogen-sensitive metrics are explicitly labelled
pre-H-opt rather than presented as results for the optimized structures.
"""
from __future__ import annotations

import argparse
import copy
import csv
import json
import math
from pathlib import Path

from .paths import native_path


SCORE_FIELDS = {
    "rosetta_fa_atr_REU": "fa_atr",
    "rosetta_fa_rep_REU": "fa_rep",
    "rosetta_fa_elec_REU": "fa_elec",
    "rosetta_solvation_REU": "fa_sol",
    "rosetta_fa_pair_REU": "fa_pair",
    "rosetta_hbond_bb_sc_REU": "hbond_bb_sc",
    "rosetta_hbond_sc_REU": "hbond_sc",
}
NETWORK_KEYS = {
    "network_added_pair_type_events", "network_removed_pair_type_events",
    "network_changed_unique_pairs_energy_delta_REU",
}
H_COLUMNS = [
    {"key": "before_interface_REU", "label": "Before interface (REU)", "kind": "energy"},
    {"key": "before_interface_delta_REU", "label": "Before interface delta (REU)", "kind": "delta_energy"},
    {"key": "hydrogen_interface_change_REU", "label": "H-opt interface change (REU)", "kind": "delta_energy"},
    {"key": "bound_before_REU", "label": "Bound before (REU)", "kind": "energy"},
    {"key": "bound_after_REU", "label": "Bound after (REU)", "kind": "energy"},
    {"key": "hydrogens_moved", "label": "Hydrogens moved", "kind": "integer"},
    {"key": "max_heavy_displacement_A", "label": "Max heavy displacement (A)", "kind": "number"},
]
STALE_SECTIONS = {"Arpeggio", "PLIP", "Discovery Studio", "E2 architecture"}
GROUPED_BASENAME = "batch_summary_grouped_fingerprints"


def _json(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8-sig"))


def _number(value: object, key: str, integer: bool = False) -> float | int:
    number = float(value)
    if not math.isfinite(number):
        raise ValueError(f"Non-finite score: {key}")
    if integer:
        if number < 0 or not number.is_integer():
            raise ValueError(f"Expected a nonnegative integer: {key}")
        return int(number)
    return number


def _close(actual: float, expected: float, label: str, tolerance: float = 1e-5) -> None:
    if not math.isclose(actual, expected, rel_tol=0, abs_tol=tolerance):
        raise ValueError(f"{label}: {actual} differs from {expected}")


def _score_values(record: dict, state: str) -> dict:
    values = {key: _number(record[f"{state}_{term}_REU"], key)
              for key, term in SCORE_FIELDS.items()}
    values["rosetta_interface_REU"] = _number(record[f"interface_{state}_REU"], "interface")
    values["rosetta_hbond_total"] = _number(record[f"{state}_ligand_hbonds"], "ligand_hbonds", True)
    values["rosetta_interface_delta_REU"] = (0.0 if state == "parent"
                                            else _number(record["interface_delta_REU"], "interface_delta"))
    _close(sum(values[key] for key in SCORE_FIELDS), values["rosetta_interface_REU"],
           "Displayed interface decomposition does not sum")
    return values


def report_columns(columns: list[dict]) -> list[dict]:
    """Preserve grouped layout while making stale-coordinate sections explicit."""
    result = []
    for original in columns:
        column = copy.deepcopy(original)
        if column["section"] in STALE_SECTIONS:
            column["section"] += " (pre-H-opt)"
        if column["key"] in NETWORK_KEYS:
            column["section"] = "Network (pre-H-opt)"
            column["label"] = "Prior " + column["label"]
        elif column["section"] == "Structure & network":
            column["section"] = "Structure"
        result.append(column)
    # Keep the current quality warning with current structural comparisons,
    # before the three explicitly historical network columns.
    quality = [column for column in result if column["key"] == "comparison_quality_warning"]
    result = [column for column in result if column["key"] != "comparison_quality_warning"]
    if quality:
        insert = next((i for i, c in enumerate(result) if c["key"] in NETWORK_KEYS), len(result))
        result[insert:insert] = quality
    insert = next((i for i, c in enumerate(result)
                   if c["section"].startswith("Arpeggio")), len(result))
    result[insert:insert] = [dict(column, section="Hydrogen optimization") for column in H_COLUMNS]
    keys = [column["key"] for column in result]
    if len(keys) != len(set(keys)):
        raise ValueError("Duplicate grouped-report column")
    return result


def _hydrogen_values(record: dict, state: str) -> dict:
    return {
        "before_interface_REU": _number(record[f"before_interface_{state}_REU"], "before_interface"),
        "before_interface_delta_REU": (0.0 if state == "parent" else
                                        _number(record["before_interface_delta_REU"], "before_delta")),
        "hydrogen_interface_change_REU": (
            _number(record[f"interface_{state}_REU"], "interface")
            - _number(record[f"before_interface_{state}_REU"], "before_interface")),
        "bound_before_REU": _number(record[f"before_{state}_bound_total_REU"], "bound_before"),
        "bound_after_REU": _number(record[f"{state}_bound_total_REU"], "bound_after"),
    }


def _parent_values(saved: dict, record: dict, protocol: dict) -> dict:
    if saved.get("status") != "complete" or saved["identity"]["protocol"] != protocol:
        raise ValueError("Parent optimization protocol differs from the batch")
    values = _score_values(record, "parent")
    values.update(_hydrogen_values(record, "parent"))
    values["hydrogens_moved"] = _number(saved["optimization"]["moved_hydrogen_count"], "parent_hydrogens_moved", True)
    values["max_heavy_displacement_A"] = _number(saved["optimization"]["max_heavy_displacement_A"], "parent_heavy_displacement")
    metrics, before = saved["metrics"], saved["before_metrics"]
    _close(values["rosetta_interface_REU"], metrics["interface_delta_X"], "Parent interface")
    for key, term in SCORE_FIELDS.items():
        _close(values[key], metrics[f"if_X_{term}"], f"Parent {term}")
    _close(values["before_interface_REU"], before["interface_delta_X"], "Parent before interface")
    _close(values["bound_before_REU"], before["bound_total"], "Parent bound before")
    _close(values["bound_after_REU"], metrics["bound_total"], "Parent bound after")
    if values["rosetta_hbond_total"] != saved["ligand_hbond_count"]:
        raise ValueError("Parent H-bond count mismatch")
    return values


def prepare_report(source: Path, output: Path) -> dict:
    """Write grouped JSON/CSV only; score calculations and source files untouched."""
    source, output = Path(source).resolve(), Path(output).resolve()
    if source == output or source in output.parents:
        raise ValueError("Hydrogen report must be outside the source report directory")
    original = _json(source / f"{GROUPED_BASENAME}.json")
    source_manifest = _json(source / "rescore_manifest.json")
    manifest = _json(output / "hydrogen_rescore_manifest.json")
    if (manifest.get("status") != "complete" or manifest.get("failure_count") != 0
            or source_manifest.get("status") != "complete" or source_manifest.get("failure_count") != 0):
        raise ValueError("Complete, failure-free source and hydrogen rescores are required")
    protocol = manifest["protocol"]
    if (protocol["source_protocol"] != source_manifest["protocol"]
            or protocol["source_protocol"].get("weights") != "ligand"
            or protocol["source_protocol"].get("legacy_pre_talaris") is not True
            or protocol["source_protocol"].get("ligand_chain") != "X"
            or protocol.get("passes_per_complex") != 1
            or manifest.get("fingerprints_regenerated") is not False):
        raise ValueError("Unexpected scoring/optimization/fingerprint protocol")
    with (output / "rosettaligand_hydrogen_optimized_rescore.csv").open(encoding="utf-8-sig", newline="") as stream:
        scores = list(csv.DictReader(stream))
    old_rows = {row["candidate"]: row for row in original["rows"]}
    names = [row["candidate"] for row in scores]
    if (not scores or len(old_rows) != len(original["rows"])
            or len(names) != len(set(names)) or set(names) != set(old_rows)):
        raise ValueError("Hydrogen-rescored candidates must match the source workbook exactly")
    if any(count != len(scores) for count in (manifest["candidate_count"], manifest["completed_count"], source_manifest["selected_count"])):
        raise ValueError("Manifest candidate count mismatch")
    parent_saved = _json(output / "parent/result.json")
    parent = _parent_values(parent_saved, scores[0], protocol)
    _close(original["parent"]["rosetta_interface_REU"], parent["before_interface_REU"], "Source parent before interface")
    if original["parent"]["rosetta_hbond_total"] != parent_saved["before_ligand_hbond_count"]:
        raise ValueError("Source parent before H-bond count mismatch")
    changed = {*SCORE_FIELDS, "rank", "rosetta_interface_REU", "rosetta_interface_delta_REU", "rosetta_hbond_total", "comparison_quality_warning"}
    new_rows = []
    for score in scores:
        name = score["candidate"]
        if score["weights"] != "ligand" or score["status"] != "complete":
            raise ValueError(f"Unexpected score weights or status: {name}")
        if _score_values(score, "parent") != {key: parent[key] for key in _score_values(score, "parent")}:
            raise ValueError(f"Parent score differs between candidates: {name}")
        if _hydrogen_values(score, "parent") != {key: parent[key] for key in _hydrogen_values(score, "parent")}:
            raise ValueError(f"Parent before/bound scores differ between candidates: {name}")
        row = copy.deepcopy(old_rows[name])
        if row["mutations"] != score["substitutions"]:
            raise ValueError(f"Substitution mismatch: {name}")
        _close(row["rosetta_interface_REU"], _number(score["before_interface_mutant_REU"], "before_interface"), "Source mutant before interface")
        if row["rosetta_hbond_total"] != _number(score["before_mutant_ligand_hbonds"], "before_ligand_hbonds", True):
            raise ValueError(f"Source mutant before H-bond count mismatch: {name}")
        row.update(_score_values(score, "mutant"))
        row.update(_hydrogen_values(score, "mutant"))
        row["hydrogens_moved"] = _number(score["moved_hydrogen_count"], "hydrogens_moved", True)
        row["max_heavy_displacement_A"] = _number(score["max_heavy_displacement_A"], "heavy_displacement")
        if not 0 <= row["max_heavy_displacement_A"] <= 1e-6:
            raise ValueError(f"Heavy atoms moved: {name}")
        # Same documented heuristic as rescore_network.comparison_warning;
        # implemented without importing structure/scoring dependencies.
        row["comparison_quality_warning"] = (
            abs(row["bound_after_REU"] - parent["bound_after_REU"]) > 100
            or row["rosetta_fa_rep_REU"] - parent["rosetta_fa_rep_REU"] > 10)
        _close(row["rosetta_interface_REU"] - parent["rosetta_interface_REU"], row["rosetta_interface_delta_REU"], "Mutant-minus-parent interface delta")
        _close(row["before_interface_REU"] - parent["before_interface_REU"], row["before_interface_delta_REU"], "Before mutant-minus-parent interface delta")
        _close(row["hydrogen_interface_change_REU"], _number(score["hydrogen_optimization_interface_change_REU"], "interface_change"), "Hydrogen interface change")
        if any(row[key] != value for key, value in old_rows[name].items() if key not in changed):
            raise AssertionError("A reused geometry/network/confidence metric changed")
        new_rows.append(row)
    new_rows.sort(key=lambda row: (row["rosetta_interface_delta_REU"], row["candidate"]))
    for rank, row in enumerate(new_rows, 1):
        row["rank"] = rank
    columns = report_columns(original["columns"])
    payload = copy.deepcopy(original)
    payload.update(
        title="m7 mutants: hydrogen-optimized RosettaLigand rescore",
        note="Ranked by H-optimized RosettaLigand mutant-minus-H-optimized-parent interface energy. Lower is more favorable; REU are not experimental binding free energies.",
        columns=columns, rows=new_rows,
        source_note=("One fixed-seed supported-proton-chi optimization pass per complex, including parent m7 once; heavy atoms and chemical states fixed. "
                     "Current RosettaLigand ligand.wts energies, Rosetta H-bonds and quality warnings were recalculated from optimized structures. "
                     "Arpeggio, PLIP, Discovery Studio, E2 architecture and mutation-network metrics are copied from the pre-H-optimization report, not rerun; "
                     "their hydrogen-sensitive counts and network energies must not be interpreted as optimized-structure results. Heavy-atom RMSDs, pLDDT, Boltz-2 and GNINA values are unchanged."),
        summary_subtitle=(f"{len(scores):,} hydrogen-optimized mutants; parent m7 optimized once. Current Rosetta scores; explicitly marked pre-H-opt fingerprints/network. "
                          f"{source_manifest['excluded_source_failures']:,} source-validation exclusions retained separately."),
        source_exclusions=source_manifest["excluded_source_failures"],
        prior_report_source_note=original.get("source_note", ""),
    )
    payload["parent"].update(parent)
    payload["scoring"] = {
        "weights": "ligand", "legacy_pre_talaris": True,
        "source_report": str(source), "manifest": str(output / "hydrogen_rescore_manifest.json"),
        "optimization": protocol, "fingerprints_regenerated": False,
        "source_batch": str(native_path(manifest["source_batch"])),
        "rejections_csv": str(source / "batch_summary.csv"),
        "ligand_validation_csv": str(source / "ligand_validation_diagnostics.csv"),
    }
    payload["column_provenance"] = {
        column["key"]: ("pre-H-optimization report; not recalculated" if "pre-H-opt" in column["section"]
                        else "optimized structures; current Rosetta scoring" if column["section"] in {"RosettaLigand", "Hydrogen optimization"}
                        or column["key"] == "comparison_quality_warning"
                        else "unchanged source heavy-atom/model-confidence metric")
        for column in columns
    }
    payload["column_definitions"] = {
        "before_interface_REU": "Interface energy on the original prepared complex, using the same RosettaLigand weights.",
        "before_interface_delta_REU": "Original mutant interface energy minus original parent interface energy.",
        "hydrogen_interface_change_REU": "Optimized minus original interface energy for this same complex. Negative favors its interface; optimization selected by bound-total energy, not this metric.",
        "bound_before_REU": "Whole bound-complex score before hydrogen optimization; not binding affinity.",
        "bound_after_REU": "Whole bound-complex score from the saved optimized PDB; not binding affinity.",
        "hydrogens_moved": "Number of explicit hydrogen atoms whose coordinates moved by more than 1e-6 A. Standard supported proton chis were sampled; other H coordinates were not freely minimized.",
        "max_heavy_displacement_A": "Maximum heavy-atom displacement between input and saved optimized PDB; required <=1e-6 A.",
        "comparison_quality_warning": "Recomputed current-score heuristic: absolute mutant-minus-parent bound-total score >100 REU, or mutant-minus-parent interface fa_rep >10 REU. Not a statistical significance test.",
        "network_added_pair_type_events": "Prior pre-H-opt mutation-network pair/type events added versus the original parent; not recalculated after hydrogen optimization.",
        "network_removed_pair_type_events": "Prior pre-H-opt mutation-network pair/type events removed versus the original parent; not recalculated after hydrogen optimization.",
        "network_changed_unique_pairs_energy_delta_REU": "Prior RosettaLigand mutation-network energy delta on pre-H-opt coordinates; not the current optimized network energy.",
    }
    (output / f"{GROUPED_BASENAME}.json").write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    with (output / f"{GROUPED_BASENAME}.csv").open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=[column["key"] for column in columns], extrasaction="ignore")
        writer.writeheader()
        writer.writerows(new_rows)
    return {"rows": len(new_rows), "columns": len(columns), "source_exclusions": payload["source_exclusions"], "output": str(output)}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("source", type=Path, help="Previous completed RosettaLigand grouped-report directory")
    parser.add_argument("output", type=Path, help="Completed hydrogen-rescore directory; grouped files written here")
    args = parser.parse_args()
    print(json.dumps(prepare_report(args.source, args.output)))


if __name__ == "__main__":
    main()
