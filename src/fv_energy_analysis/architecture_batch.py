"""Refresh DS-only E2 architecture from completed single-pose batch artifacts.

No scoring, hydrogen optimization or contact detector is executed. Unchanged
scores/fingerprints/confidence/model data are copied by candidate identity.
"""
from __future__ import annotations

import argparse
from concurrent.futures import ProcessPoolExecutor
from datetime import datetime, timezone
from hashlib import sha256
import json
from pathlib import Path
import time
import traceback

from .atom_roles import load_atom_roles
from .paths import native_path
from .report import write_csv, write_json


def digest(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def read_json(path: Path):
    return json.loads(path.read_text(encoding="utf-8-sig"))


def atom_snapshot(path: Path) -> dict:
    """Compare all named atoms including H, not PDB formatting or serials."""
    result = {}
    for line in path.read_text(encoding="utf-8").splitlines():
        if line.startswith("ENDMDL"):
            break
        if not line.startswith(("ATOM  ", "HETATM")):
            continue
        key = (line[21:22], line[22:27], line[17:20], line[12:16].strip())
        if key in result:
            raise ValueError(f"Duplicate named atom {key!r} in {path}")
        result[key] = (line[76:78].strip(), *[float(line[i:i+8]) for i in (30, 38, 46)])
    if not result:
        raise ValueError(f"No atoms in {path}")
    return result


def verify_opth(prepared: Path, record_path: Path, params_hash: str) -> dict:
    record = read_json(record_path)
    optimized = native_path(record["optimized_pdb"])
    expected_hash = record["artifact_sha256"][record["optimized_pdb"]]
    if record.get("status") != "complete" or digest(optimized) != expected_hash:
        raise ValueError(f"OptH output failed saved artifact hash/status validation: {optimized}")
    protocol = record["identity"]["protocol"]
    source = protocol["source_protocol"]
    if (protocol["passes_per_complex"] != 1 or source["weights"] != "ligand"
            or not source["legacy_pre_talaris"] or source["params_sha256"] != params_hash):
        raise ValueError(f"Incompatible single-OptH RosettaLigand protocol: {record_path}")
    if atom_snapshot(prepared) != atom_snapshot(optimized):
        raise ValueError(f"DS prepared structure differs from saved OptH atoms/coordinates: {prepared}")
    return {"optimized_pdb": str(optimized), "optimized_pdb_sha256": expected_hash,
            "prepared_pdb": str(prepared), "prepared_pdb_sha256": digest(prepared),
            "all_named_atoms_and_coordinates_identical": True, "seed": protocol["seed"],
            "replicates": 1, "opth_record": str(record_path)}


def render_report(candidate: str, result: dict, parent: dict, original: Path) -> str:
    from .fingerprint_reports import ARCHITECTURE_COLUMNS
    lines = [f"# E2 recognition: {candidate}", "",
             "One existing OptH pose. Interaction classifications: Discovery Studio only. "
             "Steric sandwich: project-calibrated coordinate geometry, not a native DS class.", "",
             f"Original complete analysis: [report]({(original / 'automatic_summary.md').as_posix()}). "
             "Use this report for updated E2 recognition flags; the original architecture section is historical.", "",
             "| Metric | Parent m7 | Mutant |", "|---|---|---|"]
    for key, label, _ in ARCHITECTURE_COLUMNS:
        lines.append(f"| {label} | {parent.get(key, 'unavailable')} | {result.get(key, 'unavailable')} |")
    lines.extend(["", "## Definition", "",
                  "O3/O17: DS conventional H-bonds involving the mapped oxygen or its attached hydroxyl hydrogen. "
                  "Carbon H-bonds do not satisfy those flags. Pi-mediated contacts combine DS stacked/T-shaped "
                  "pi-pi, pi-cation and pi-anion partners, counted once per protein residue. They are a descriptor, "
                  "not an additional prerequisite for the four-feature screen.", "",
                  "Steric sandwich (two-sided A-ring enclosure): distinct side-chain residues provide "
                  "opposite-face heavy atoms with nearest A-ring-carbon surface gaps 0–1.0 A. Each atom "
                  "lies at least its own van der Waals radius from the ring plane, and their joining segment "
                  "crosses the projected union of the six carbon van der Waals disks (radius 1.70 A). "
                  "All qualifying atoms are searched; residue pairs are counted once. Overlaps are "
                  "reported separately and do not qualify as close contacts. This is not SASA or a C19-exclusion test.", "",
                  "C18 definition and calibration, matched DS records and wall atom identities are retained in "
                  "the accompanying JSON. Passing these geometric features does not establish improved affinity.", ""])
    return "\n".join(lines)


def _one(task: dict) -> dict:
    from .ds_architecture import VERSION, analyze_ds_architecture
    from .steric_sandwich import METHOD_VERSION
    try:
        folder = Path(task["source_result"])
        summary = read_json(folder / "final" / "summary.json")
        role = task["role"]
        workflow_path = folder / "workflow_manifest.json"
        if role == "mutant":
            score_protocol = read_json(workflow_path).get("scoring_protocol", {})
            if score_protocol.get("weights") != "ligand" or score_protocol.get("params_sha256") != task["params_hash"]:
                raise ValueError(f"Copied full-analysis scores have incompatible weights/params: {workflow_path}")
        pdb = native_path(summary[role])
        if not pdb.is_file():
            pdb = folder / "prepared" / pdb.name
        ds = Path(task["ds_csv"])
        cache_path = Path(task["cache_path"])
        if task.get("resume") and cache_path.is_file():
            cached = read_json(cache_path)
            previous = cached.get("result", {})
            provenance = previous.get("provenance", {})
            if (cached.get("status") == "complete" and cached.get("candidate") == task["candidate"]
                    and previous.get("e2_architecture_definition_version") == VERSION
                    and previous.get("e2_steric_sandwich", {}).get("method") == METHOD_VERSION
                    and previous.get("e2_c18_calibration") == task["c18_calibration"]
                    and provenance.get("prepared_pdb_sha256") == digest(pdb)
                    and provenance.get("ds_csv_sha256") == digest(ds)
                    and provenance.get("optimized_pdb_sha256") == digest(native_path(provenance["optimized_pdb"]))):
                cached["reused"] = True
                return cached
        provenance = verify_opth(pdb, Path(task["opth_record"]), task["params_hash"])
        # DS did not save an input digest in the historical run. Check its
        # retained role/candidate/status associations, without overstating that
        # chain of provenance as a cryptographic DS-input guarantee.
        import csv
        with ds.open(encoding="utf-8-sig", newline="") as handle:
            rows = list(csv.DictReader(handle))
        if len(rows) != int(task["ds_status"]["interaction_count"]):
            raise ValueError(f"DS row count differs from completed detector status: {ds}")
        for row in rows:
            if (row["role"] != role or row["candidate"] != task["candidate"]
                    or row["structure_id"] != task["ds_status"]["structure_id"]
                    or native_path(row["result_folder"]).resolve() != folder.resolve()):
                raise ValueError(f"Mismatched DS structure/candidate/result association in {ds}")
        result = analyze_ds_architecture(
            pdb, ds, roles=task["roles"], ligand_resname="E2L", ligand_chain="X",
            fv_chains=set(summary["fv_chains"]), c18_calibration=task["c18_calibration"],
        )
        if result.get("e2_architecture_status") != "complete":
            raise ValueError(f"Architecture unavailable: {result.get('e2_architecture_status')}")
        result["provenance"] = {
            **provenance, "source_result": str(folder), "ds_csv": str(ds),
            "ds_csv_sha256": digest(ds), "ds_status": task["ds_status"],
            "ds_add_hydrogens": False, "ds_input_digest_recorded_by_original_run": False,
            "ds_reuse_validation": "saved detector associations and OptH/prepared atom-coordinate identity",
            "arpeggio_used_for_architecture": False,
        }
        return {"candidate": task["candidate"], "status": "complete", "result": result, "reused": False,
                "source_result": str(folder)}
    except Exception as error:
        return {"candidate": task["candidate"], "status": "failed",
                "error_type": type(error).__name__, "error": str(error),
                "source_result": task["source_result"], "traceback": traceback.format_exc()}


def run(args) -> Path:
    from .ds_architecture import VERSION
    from .fingerprint_reports import ARCHITECTURE_COLUMNS
    from .fingerprint_counts import _read_csv
    from .steric_sandwich import METHOD_VERSION
    started = time.monotonic()
    source, output = args.source.resolve(), args.output.resolve()
    if source == output or source in output.parents:
        raise ValueError("Use a separate output folder, outside the completed source batch")
    output.mkdir(parents=True, exist_ok=True)
    previous_manifest = read_json(output / "batch_manifest.json") if (output / "batch_manifest.json").is_file() else None
    if previous_manifest is not None and not args.resume:
        raise FileExistsError(f"Output already exists; choose a new folder: {output}")
    manifest = read_json(source / "batch_manifest.json")
    payload = read_json(source / "batch_summary_grouped_fingerprints.json")
    if (payload.get("scoring_profile") != "rosettaligand_legacy"
            or payload.get("scoring", {}).get("weights") != "ligand"
            or payload.get("scoring", {}).get("legacy_pre_talaris") is not True):
        raise ValueError("The copied full-metrics table must use legacy RosettaLigand scoring")
    if (manifest["status"] != "complete" or manifest.get("failure_count")
            or len(payload["rows"]) != manifest["completed_count"]):
        raise ValueError("Source must be a complete, internally consistent scientific batch")
    if args.expected_count and len(payload["rows"]) != args.expected_count:
        raise ValueError("Candidate count differs from requested batch")
    roles = load_atom_roles(args.params, args.template)
    reference_calibration = read_json(args.c18_calibration)
    # This file is produced by the isolated reference study, not inferred from
    # the candidate cohort or a permissive geometry-only fallback.
    if not reference_calibration.get("production_rule"):
        raise ValueError("C18 calibration lacks an explicit production_rule")
    write_json(output / "c18_reference_calibration.json", reference_calibration)
    calibration = {key: reference_calibration[key] for key in (
        "protocol", "production_rule", "reference_paper", "published_concept",
        "reference_copies", "families", "singleton_DS_presence_pass", "singleton_count_min",
        "singleton_count_max", "controls_remaining_pass", "numerical_rule_note")}
    calibration.update(source_calibration_sha256=digest(args.c18_calibration),
                       source_calibration_json=str(args.c18_calibration.resolve()))
    ds_info = manifest.get("discovery_studio", {})
    if ds_info.get("status") != "complete" or ds_info.get("add_hydrogens") is not False:
        raise ValueError("Completed Discovery Studio with AddHydrogens=False is required")
    statuses = _read_csv(source / "discovery_studio_status.csv")
    status_by = {row["candidate"]: row for row in statuses}
    if len(status_by) != len(statuses) or any(row["status"] != "complete" for row in statuses):
        raise ValueError("Incomplete or duplicate Discovery Studio detector status")
    items = {native_path(item["candidate"]).name: item for item in manifest["candidates"]}
    names = [row["candidate"] for row in payload["rows"]]
    if len(set(names)) != len(names) or set(names) != set(items) or set(status_by) != {*names, "parent"}:
        raise ValueError("Report, manifest and Discovery Studio candidate identities differ")
    params_hash = digest(args.params)
    if previous_manifest is not None:
        previous_payload = read_json(output / "batch_summary_grouped_fingerprints.json")
        if (previous_manifest.get("source_manifest_sha256") != digest(source / "batch_manifest.json")
                or previous_manifest.get("source_payload_sha256") != digest(source / "batch_summary_grouped_fingerprints.json")
                or previous_manifest.get("params_sha256") != params_hash
                or previous_manifest.get("template_sha256") != digest(args.template)
                or previous_payload.get("architecture_protocol", {}).get("c18_calibration_sha256") != digest(args.c18_calibration)):
            raise ValueError("Resume source/parameters/calibration differs from the saved batch")
    common = {"params_hash": params_hash, "roles": roles, "c18_calibration": calibration}
    parent_cache = native_path(manifest["optimization"]["parent_cache"])
    parent_task = {**common, "role": "parent", "candidate": "parent", "source_result": str(parent_cache),
                   "ds_csv": str(source / "discovery_studio_parent_interactions.csv"),
                   "ds_status": status_by["parent"], "opth_record": str(args.opth_results / "parent/result.json"),
                   "resume": False, "cache_path": str(output / "parent_architecture.json")}
    parent = _one(parent_task)
    if parent["status"] != "complete":
        write_json(output / "parent_failure.json", parent)
        raise RuntimeError(f"Parent preflight failed: {parent['error']}")
    parent_result = parent["result"]
    write_json(output / "parent_architecture.json", parent_result)
    tasks = [{**common, "role": "mutant", "candidate": name,
              "source_result": str(native_path(items[name]["result_folder"])),
              "ds_csv": str(native_path(items[name]["result_folder"]) / "discovery_studio_interactions.csv"),
              "ds_status": status_by[name],
              "resume": args.resume, "cache_path": str(output / Path(name).stem / "estrogen_binding_architecture.json"),
              "opth_record": str(args.opth_results / "cases" / Path(name).stem / "result.json")}
             for name in names]
    state_path = output / "architecture_status.json"
    write_json(state_path, {"status": "running", "candidate_count": len(names), "jobs": args.jobs,
                            "replicates": 1, "source": str(source), "arpeggio_used": False})
    results = []
    with ProcessPoolExecutor(max_workers=args.jobs) as pool:
        for record in pool.map(_one, tasks, chunksize=8):
            results.append(record)
            case = output / Path(record["candidate"]).stem
            case.mkdir(exist_ok=True)
            if not record.get("reused"):
                write_json(case / "estrogen_binding_architecture.json", record)
            if record["status"] == "complete" and not record.get("reused"):
                (case / "architecture_summary.md").write_text(render_report(
                    record["candidate"], record["result"], parent_result, Path(record["source_result"])), encoding="utf-8")
    by_name = {record["candidate"]: record for record in results}
    failures = [record for record in results if record["status"] != "complete"]
    old_columns = payload["columns"]
    new_arch_columns = [{"key": key, "label": label, "section": "E2 architecture", "kind": kind,
                         "width": 35 if kind == "text" else 23}
                        for key, label, kind in ARCHITECTURE_COLUMNS]
    insertion = next(i for i, column in enumerate(old_columns) if column["section"] == "E2 architecture")
    columns = [column for column in old_columns if column["section"] != "E2 architecture"]
    columns[insertion:insertion] = new_arch_columns
    architecture_keys = [key for key, _, _ in ARCHITECTURE_COLUMNS]
    old_arch_keys = {column["key"] for column in old_columns if column["section"] == "E2 architecture"}
    for row in payload["rows"]:
        record = by_name[row["candidate"]]
        for key in old_arch_keys:
            row.pop(key, None)
        result = record["result"] if record["status"] == "complete" else {
            "e2_architecture_status": "unavailable: " + record["error"]}
        row.update({key: result.get(key) for key in architecture_keys})
    for key in old_arch_keys:
        payload["parent"].pop(key, None)
    payload["parent"].update({key: parent_result.get(key) for key in architecture_keys})
    payload["columns"] = columns
    payload["architecture_protocol"] = {"version": VERSION, "steric_method": METHOD_VERSION, "replicates": 1,
        "source_batch": str(source), "interaction_source": "Discovery Studio",
        "steric_source": "paper-reference-calibrated project coordinate proxy",
        "c18_calibration_sha256": digest(args.c18_calibration), "arpeggio_used": False}
    payload["summary_subtitle"] = f"{len(names)} single-OptH candidates. DS-only recognition metrics and calibrated steric geometry. Other results reused unchanged."
    payload["source_note"] = ("E2: DS conventional O3/O17 H-bonds, reference-calibrated DS C18 packing, and two-sided A-ring enclosure geometry. "
        "Pi-mediated partners are distinct residues (stacked/T-shaped pi-pi, pi-cation, pi-anion), not an extra pass requirement. "
        "Nahoum et al. (2003), doi:10.1096/fj.02-0524fje. Project numerical cutoffs are calibrated proxies, not published affinity thresholds.")
    write_json(output / "batch_summary_grouped_fingerprints.json", payload)
    write_csv(output / "batch_summary_grouped_fingerprints.csv", [{column["key"]: row.get(column["key"])
        for column in columns} for row in payload["rows"]])
    write_csv(output / "e2_architecture_screen.csv", [{"candidate": record["candidate"],
        **{key: record.get("result", {}).get(key) for key in architecture_keys}} for record in results])
    counts = {key: sum(record.get("result", {}).get(key) is True for record in results)
              for key in ("e2_o3_hbond", "e2_o17_hbond", "e2_c18_hydrophobic_contact",
                          "e2_geometric_sandwich", "e2_strict_screen_pass")}
    output_manifest = {"status": "complete" if not failures else "architecture_failed",
        "candidate_count": len(names), "completed_count": len(names), "failure_count": 0,
        "architecture_failure_count": len(failures), "source_batch": str(source),
        "source_manifest_sha256": digest(source / "batch_manifest.json"),
        "source_payload_sha256": digest(source / "batch_summary_grouped_fingerprints.json"),
        "params_sha256": params_hash, "template_sha256": digest(args.template),
        "c18_calibration": str(args.c18_calibration.resolve()), "feature_counts": counts,
        "architecture_definition_version": VERSION, "steric_method": METHOD_VERSION,
        "replicates": 1, "jobs": args.jobs, "detectors_rerun": [], "scoring_rerun": False,
        "cached_architectures_reused": sum(record.get("reused", False) for record in results),
        "architectures_calculated_this_run": sum(not record.get("reused", False) for record in results),
        "completed_utc": datetime.now(timezone.utc).isoformat(),
        "elapsed_seconds": round(time.monotonic() - started, 3),
        "candidates": [{"candidate": record["candidate"], "result_folder": str(output / Path(record["candidate"]).stem),
                         "source_result": record["source_result"], "status": record["status"]} for record in results]}
    write_json(output / "batch_manifest.json", output_manifest)
    write_json(output / "architecture_failures.json", failures)
    write_json(state_path, {key: value for key, value in output_manifest.items() if key != "candidates"})
    (output / "batch_summary.md").write_text("# DS-only E2 recognition update\n\n"
        f"{len(names)} existing single-OptH mutant poses, with parent m7 evaluated once. "
        f"Architecture failures: {len(failures)}. No hydrogen optimization, scoring or detector was rerun.\n\n"
        + "\n".join(f"- {key}: {value}/{len(names)}" for key, value in counts.items())
        + "\n\nThe workbook/CSV retain all original non-architecture metrics unchanged. "
        "Pi-mediated contacts are a descriptor; all-four recognition requires O3, O17, C18 packing and the steric proxy "
        "in this same pose. This is not an affinity or E2-versus-E3 specificity proof.\n", encoding="utf-8")
    from .resources import resource_path
    (output / "E2_DEFINITIONS.md").write_text(
        resource_path("docs/e2_architecture.md").read_text(encoding="utf-8"), encoding="utf-8")
    if failures:
        raise RuntimeError(f"{len(failures)} architecture failures; detailed errors saved without silent false flags")
    return output


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("source", "output", "params", "template", "opth-results", "c18-calibration"):
        parser.add_argument("--" + name, type=Path, required=True)
    parser.add_argument("--jobs", type=int, choices=range(1, 9), default=3)
    parser.add_argument("--expected-count", type=int)
    parser.add_argument("--resume", action="store_true", help="Reuse verified completed architecture cases; retry failures only")
    parser.add_argument("--xlsx", action="store_true", help="Also export the grouped workbook; requires the reporting extra")
    args = parser.parse_args()
    output = run(args)
    if args.xlsx:
        from .spreadsheet import export_workbook
        export_workbook(read_json(output / "batch_summary_grouped_fingerprints.json"),
                        output / "batch_summary_grouped_fingerprints.xlsx")
    print(output, flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
