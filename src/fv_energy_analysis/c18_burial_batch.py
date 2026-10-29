"""Compare saved DS C18 contacts with C18 solvent burial without rescoring.

Inputs are a combined batch and its candidate_sources.csv. Every measured PDB
is resolved from the saved architecture provenance and checked against its hash.
The production batch is read only; this module writes a separate audit folder.
"""
from __future__ import annotations

import argparse
from collections import Counter
from concurrent.futures import ProcessPoolExecutor
from datetime import datetime, timezone
from hashlib import sha256
import csv
import json
from pathlib import Path
import time

import numpy as np

from .paths import native_path


def read_json(path):
    return json.loads(Path(path).read_text(encoding="utf-8-sig"))


def digest(path):
    return sha256(Path(path).read_bytes()).hexdigest()


def write_json(path, data):
    Path(path).write_text(json.dumps(data, indent=2, allow_nan=False) + "\n", encoding="utf-8")


def average_ranks(values):
    """One-based average ranks, including ties, without a scipy dependency."""
    values = np.asarray(values, dtype=float)
    order = np.argsort(values, kind="stable")
    ranked = np.empty(len(values), dtype=float)
    start = 0
    while start < len(values):
        end = start + 1
        while end < len(values) and values[order[end]] == values[order[start]]:
            end += 1
        ranked[order[start:end]] = (start + 1 + end) / 2
        start = end
    return ranked


def correlation(x, y):
    x, y = np.asarray(x, dtype=float), np.asarray(y, dtype=float)
    if len(x) < 2 or np.ptp(x) == 0 or np.ptp(y) == 0:
        return {"coefficient": None, "reason": "undefined: at least one variable is constant"}
    return {"coefficient": float(np.corrcoef(x, y)[0, 1]), "reason": None}


def describe(rows):
    values = np.array([r["percent_buried"] for r in rows], dtype=float)
    if not len(values):
        return {"n": 0}
    return {
        "n": len(values), "mean_percent": float(values.mean()),
        "sd_percent": float(values.std(ddof=1)) if len(values) > 1 else None,
        "min_percent": float(values.min()), "p05_percent": float(np.quantile(values, .05)),
        "q1_percent": float(np.quantile(values, .25)), "median_percent": float(np.median(values)),
        "q3_percent": float(np.quantile(values, .75)), "max_percent": float(values.max()),
        "fully_buried_n": sum(r["complex_SASA_A2"] == 0 for r in rows),
    }


def summarize(rows):
    valid = [r for r in rows if r["status"] == "complete"]
    burial = [r["percent_buried"] for r in valid]
    flags = [int(r["ds_c18_packing"]) for r in valid]
    counts = [r["ds_c18_partner_count"] for r in valid]
    groups = {str(flag).lower(): describe([r for r in valid if r["ds_c18_packing"] == flag])
              for flag in (False, True)}
    return {
        "candidate_count": len(rows), "measured_count": len(valid),
        "excluded_input_rejection_count": sum(r["status"] == "excluded_input_rejection" for r in rows),
        "measurement_failure_count": sum(r["status"] == "measurement_failed" for r in rows),
        "ds_flag_counts": dict(Counter(str(bool(v)).lower() for v in flags)),
        "ds_partner_count_distribution": dict(sorted(Counter(counts).items())),
        "burial_resolution_note": "Zero sampled complex SASA is a finite-resolution ceiling, not an infinite-resolution proof.",
        "overall_burial": describe(valid), "by_ds_flag": groups,
        "by_ds_partner_count": {str(k): describe([r for r in valid if r["ds_c18_partner_count"] == k])
                                for k in sorted(set(counts))},
        "point_biserial_flag_vs_burial": correlation(flags, burial),
        "spearman_partner_count_vs_burial": correlation(average_ranks(counts), average_ranks(burial)),
        "pearson_partner_count_vs_burial": correlation(counts, burial),
        "mean_burial_difference_true_minus_false_percentage_points": (
            groups["true"].get("mean_percent", 0) - groups["false"].get("mean_percent", 0)
            if groups["true"]["n"] and groups["false"]["n"] else None),
        "fully_buried_but_ds_false_count": sum(r["complex_SASA_A2"] == 0 and not r["ds_c18_packing"]
                                                   for r in valid),
        "not_fully_buried_but_ds_true_count": sum(r["complex_SASA_A2"] > 0 and r["ds_c18_packing"]
                                                  for r in valid),
        "inferential_p_values": None,
        "inferential_note": "Descriptive analysis of the full related-mutant panel; structures are not independent experimental replicates.",
    }


def resolve_tasks(source):
    payload_path = source / "batch_summary_grouped_fingerprints.json"
    payload = read_json(payload_path)
    with (source / "candidate_sources.csv").open(encoding="utf-8-sig", newline="") as handle:
        associations = list(csv.DictReader(handle))
    lookup = {r["candidate"]: native_path(r["source"]) for r in associations}
    if len(lookup) != len(associations) or set(lookup) != {r["candidate"] for r in payload["rows"]}:
        raise ValueError("Candidate/source associations must be unique and cover the full dataset")
    tasks, excluded, parent = [], [], None
    for row in payload["rows"]:
        base = {"candidate": row["candidate"], "mutations": row.get("mutations"),
                "source_cohort": str(lookup[row["candidate"]]),
                "ds_c18_packing": row.get("e2_c18_hydrophobic_contact"),
                "ds_c18_partner_count": row.get("e2_c18_hydrophobic_partner_count")}
        if row["status"] != "complete":
            excluded.append({**base, "status": "excluded_input_rejection", "error": row.get("error")})
            continue
        folder = lookup[row["candidate"]] / Path(row["candidate"]).stem
        architecture = read_json(folder / "estrogen_binding_architecture.json")
        result = architecture.get("result", architecture.get("mutant"))
        if result is None or result.get("e2_architecture_status") != "complete":
            raise ValueError(f"Missing completed architecture for {row['candidate']}")
        for key in ("e2_c18_hydrophobic_contact", "e2_c18_hydrophobic_partner_count"):
            if result[key] != row[key]:
                raise ValueError(f"Combined/saved architecture mismatch: {row['candidate']}: {key}")
        if base["ds_c18_packing"] != (base["ds_c18_partner_count"] >= 1):
            raise ValueError(f"Inconsistent DS flag/count: {row['candidate']}")
        provenance = result["e2_provenance"]
        path = native_path(provenance["pdb"])
        if not path.is_file():
            raise FileNotFoundError(path)
        task = {**base, "pdb": str(path), "pdb_sha256": provenance["pdb_sha256"],
                "ligand_atom": result["e2_atom_roles"]["c18_methyl_atom"]}
        if task["ligand_atom"] != "C10":
            raise ValueError("Expected the validated E2L steroid-C18 mapping to file C10")
        tasks.append(task)
        if parent is None and "parent" in architecture:
            pr = architecture["parent"]
            pp = pr["e2_provenance"]
            parent = {"candidate": "parent_m7", "mutations": "Parent", "source_cohort": str(folder),
                      "pdb": str(native_path(pp["pdb"])), "pdb_sha256": pp["pdb_sha256"],
                      "ligand_atom": pr["e2_atom_roles"]["c18_methyl_atom"],
                      "ds_c18_packing": pr["e2_c18_hydrophobic_contact"],
                      "ds_c18_partner_count": pr["e2_c18_hydrophobic_partner_count"]}
    if parent is None:
        raise ValueError("No saved parent architecture found")
    return payload_path, tasks, excluded, parent


def measure(task):
    from .c18_burial import analyze_c18_burial
    try:
        before = digest(task["pdb"])
        if before != task["pdb_sha256"]:
            raise ValueError("Prepared PDB no longer matches saved DS-architecture input hash")
        value = analyze_c18_burial(task["pdb"], methyl_atom=task["ligand_atom"], n_points=task["n_points"])
        if digest(task["pdb"]) != before:
            raise ValueError("Source PDB changed during measurement")
        return {**task, **value, "status": "complete", "error": None}
    except Exception as exc:
        return {**task, "status": "measurement_failed", "error": f"{type(exc).__name__}: {exc}"}


def render_report(result):
    summary = result["summary"]
    def fmt(v):
        return "undefined (constant variable)" if v is None else f"{v:.4f}"
    lines = ["# C18 burial versus Discovery Studio packing", "",
             f"Measured {summary['measured_count']:,} mutants; kept {summary['excluded_input_rejection_count']} prior chemistry rejections separate. "
             f"New measurement failures: {summary['measurement_failure_count']}.", "",
             "The existing DS flag means at least one unique residue with a C18-specific alkyl/pi-alkyl annotation. "
             "It uses the saved native DS surface-area scale 1.00. The diagnostic 0.65 setting was not applied here.", "",
             "Burial (%) = 100 × (SASA of C18 in the complete isolated ligand − SASA of C18 in the complex) / "
             "SASA of C18 in the complete isolated ligand. Same coordinates in both states; heavy atoms; "
             f"1.4 Å solvent probe; Biopython Shrake–Rupley, {result['method']['n_points']:,} points per atom; C/N/O/S radii 1.70/1.55/1.52/1.80 Å. "
             "Protein backbone is included; hydrogens/water are excluded. Only mathematically possible protein occluders "
             "are retained during the calculation, with the full ligand in both states.", "",
             "| DS C18 flag | N | Mean burial % | Median % | Range % | Fully buried N |",
             "|---|---:|---:|---:|---:|---:|"]
    for flag in ("false", "true"):
        group = summary["by_ds_flag"][flag]
        lines.append(f"| {flag.upper()} | {group['n']} | {group['mean_percent']:.3f} | {group['median_percent']:.3f} | "
                     f"{group['min_percent']:.3f}–{group['max_percent']:.3f} | {group['fully_buried_n']} |")
    lines += ["", f"Point-biserial correlation, flag versus burial: **{fmt(summary['point_biserial_flag_vs_burial']['coefficient'])}**. "
              f"Spearman correlation, unique partner count versus burial: **{fmt(summary['spearman_partner_count_vs_burial']['coefficient'])}**.", "",
              f"Fully buried C18 but DS FALSE: **{summary['fully_buried_but_ds_false_count']}**. "
              f"DS TRUE with some remaining C18 solvent exposure: **{summary['not_fully_buried_but_ds_true_count']}**.", "",
              "| Unique DS C18 partners | N | Mean burial % | Median % | Range % | Fully buried N |",
              "|---:|---:|---:|---:|---:|---:|"]
    for count, group in summary["by_ds_partner_count"].items():
        lines.append(f"| {count} | {group['n']} | {group['mean_percent']:.3f} | {group['median_percent']:.3f} | "
                     f"{group['min_percent']:.3f}–{group['max_percent']:.3f} | {group['fully_buried_n']} |")
    parent = result["parent"]
    ceiling = summary["overall_burial"]["fully_buried_n"]
    lines += ["", f"Parent m7: {parent['percent_buried']:.3f}% buried, DS flag {parent['ds_c18_packing']}, "
              f"{parent['ds_c18_partner_count']} unique partners. Parent is excluded from mutant correlations.", "",
              "## Interpretation", "",
              f"{ceiling:,}/{summary['measured_count']:,} mutants "
              f"({100*ceiling/summary['measured_count']:.2f}%) have 100% measured C18 burial. "
              "Consequently burial saturates for much of this panel. The positive correlations show that DS "
              "contacts tend to accompany greater burial, but the fully buried DS-negative cases demonstrate "
              "that the DS flag is not a burial requirement or a replacement for the continuous measurement.", "",
              "Burial measures solvent shielding by all protein atoms, including polar atoms. DS annotations measure "
              "a particular hydrophobic interaction class; the two quantities are not interchangeable. Near-complete "
              "burial creates a ceiling that limits its ability to discriminate these poses. Neither value alone "
              "establishes favorable packing energy, stronger affinity or E2 specificity.", "",
              "This is descriptive analysis of the full panel of related mutants; no independent-replicate p-values "
              "or new screening cutoff were assigned. Fully buried means zero accessible sampled points at this "
              "resolution, not a claim of exact zero at infinite resolution. Existing structures, DS settings, "
              "energies, production CSV and workbook were not modified.", "",
              "Per-candidate values and input PDB hashes: [c18_burial_comparison.csv](c18_burial_comparison.csv). "
              "Statistics, method and provenance: [results.json](results.json).", ""]
    return "\n".join(lines)


def run(args):
    started = time.monotonic()
    source, output = args.source.resolve(), args.output.resolve()
    if output == source or source in output.parents:
        raise ValueError("Use a separate output folder outside the completed combined batch")
    output.mkdir(parents=True, exist_ok=True)
    if (output / "results.json").exists() and not args.resume:
        raise FileExistsError("Use --resume or a fresh audit folder")
    payload, tasks, excluded, parent_task = resolve_tasks(source)
    input_hash = digest(payload)
    method = {"n_points": args.n_points, "full_ligand": True, "probe_A": 1.4,
              "scope": "C18 heavy carbon only", "ds_surface_scale": 1.0}
    identity = {"payload_sha256": input_hash, "method": method,
                "measurement_code_sha256": digest(Path(__file__).with_name("c18_burial.py"))}
    checkpoint = output / "checkpoint.jsonl"
    identity_path = output / "identity.json"
    cached = {}
    if checkpoint.exists():
        if not args.resume or read_json(identity_path) != identity:
            raise ValueError("Existing checkpoint requires --resume with identical inputs/method")
        for line in checkpoint.read_text(encoding="utf-8").splitlines():
            row = json.loads(line)
            if row["status"] == "complete":
                cached[row["candidate"]] = row
    else:
        write_json(identity_path, identity)
    for task in tasks + [parent_task]:
        task["n_points"] = args.n_points
    parent = measure(parent_task)
    if parent["status"] != "complete":
        raise ValueError(f"Parent measurement failed: {parent['error']}")
    pending, measured = [], []
    for task in tasks:
        old = cached.get(task["candidate"])
        if old and old["pdb_sha256"] == task["pdb_sha256"] and digest(task["pdb"]) == task["pdb_sha256"]:
            measured.append(old)
        else:
            pending.append(task)
    print(json.dumps({"stage": "measure", "eligible": len(tasks), "reused": len(measured),
                      "pending": len(pending), "workers": args.jobs}), flush=True)
    with checkpoint.open("a", encoding="utf-8") as handle:
        with ProcessPoolExecutor(max_workers=args.jobs) as pool:
            for row in pool.map(measure, pending, chunksize=8):
                measured.append(row)
                handle.write(json.dumps(row, allow_nan=False) + "\n")
                if len(measured) % 250 == 0:
                    handle.flush()
                    print(json.dumps({"completed": len(measured), "eligible": len(tasks),
                                      "elapsed_seconds": round(time.monotonic()-started, 1)}), flush=True)
    if digest(payload) != input_hash:
        raise ValueError("Combined input dataset changed during measurement")
    rows = sorted(measured + excluded, key=lambda r: r["candidate"])
    columns = ["candidate", "mutations", "status", "ds_c18_packing", "ds_c18_partner_count",
               "percent_buried", "fraction_buried", "isolated_full_ligand_SASA_A2", "complex_SASA_A2",
               "buried_SASA_A2", "ligand_atom", "protein_occluder_atom_count", "n_points",
               "pdb", "pdb_sha256", "source_cohort", "error"]
    with (output / "c18_burial_comparison.csv").open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)
    result = {"completed_utc": datetime.now(timezone.utc).isoformat(),
              "input": str(payload), "input_sha256": input_hash, "method": method,
              "source_coordinates_changed": False, "production_outputs_changed": False,
              "elapsed_seconds": round(time.monotonic()-started, 2), "parent": parent,
              "summary": summarize(rows)}
    write_json(output / "results.json", result)
    (output / "report.md").write_text(render_report(result), encoding="utf-8")
    print(json.dumps(result["summary"], indent=2, allow_nan=False), flush=True)
    if result["summary"]["measurement_failure_count"]:
        raise RuntimeError("Burial measurements failed; see explicit per-candidate errors")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--jobs", type=int, default=3)
    parser.add_argument("--n-points", type=int, default=1920)
    parser.add_argument("--resume", action="store_true")
    args = parser.parse_args()
    if args.jobs < 1 or args.n_points < 1:
        parser.error("--jobs and --n-points must be positive")
    run(args)


if __name__ == "__main__":
    main()
