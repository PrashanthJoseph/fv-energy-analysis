"""Reclassify saved protonated E3 probes with explicit H-bond geometry.

No ligand construction, H addition, scoring or minimization is repeated. The
new OH proton may be reselected from its twelve already saved orientations.
The chosen single orientation is persisted so contacts match the exported PDB.
"""
from __future__ import annotations

import argparse
from collections import Counter
from concurrent.futures import ProcessPoolExecutor
from datetime import datetime, timezone
from hashlib import sha256
import copy
import json
from pathlib import Path
import shutil
import time
import traceback

import numpy as np

from .e3_discrimination import digest, read_json, write_flat_csv
from .e3_io import replace_hydroxyl_hydrogen
from .paths import native_path
from .report import write_json


def refresh_one(task):
    from .e3_steric import analyze_e3_steric
    from .e3_reports import render_e3_markdown
    output = Path(task["output"])
    output.mkdir(parents=True, exist_ok=True)
    artifact = output / "e3_discrimination.json"
    try:
        original_path = Path(task["artifact"])
        original_hash = digest(original_path)
        if task.get("resume") and artifact.is_file():
            previous = read_json(artifact)
            probe = native_path(previous.get("probe_pdb", ""))
            if (previous.get("status") == "complete" and previous.get("source_e3_artifact_sha256") == original_hash
                    and previous.get("refresh_identity") == task["identity"]
                    and probe.is_file() and digest(probe) == previous.get("probe_pdb_sha256")):
                (output / "automatic_summary.md").write_text(render_e3_markdown(previous), encoding="utf-8")
                previous["reused"] = True
                return previous
        original = read_json(original_path)
        if original.get("status") != "complete" or original.get("candidate") != task["candidate"]:
            raise ValueError("Saved E3 artifact status/candidate differs from the requested case")
        source_probe = native_path(original["probe_pdb"])
        probe_bytes = source_probe.read_bytes()
        if sha256(probe_bytes).hexdigest() != original["probe_pdb_sha256"]:
            raise ValueError("Saved protonated E3 probe hash differs from its preparation record")
        samples = copy.deepcopy(original["hydrogen_sampling"]["orientations"])
        if len(samples) != 12 or original["chemistry"]["ligand_hydrogen_count"] != 24:
            raise ValueError("Expected the saved twelve-orientation full-H E3 preparation")
        # Assess exactly the coordinates that can be persisted in a PDB,
        # preventing rounding from moving an angle across a classification gate.
        for sample in samples:
            sample["hydrogen_xyz_A"] = np.round(sample["hydrogen_xyz_A"], 3).tolist()
        measured = analyze_e3_steric(source_probe, hydrogen_positions=samples,
                                     overlap_warning_A=task["overlap_warning_A"])
        selected_xyz = measured.get("selected_hydrogen_xyz_A")
        if selected_xyz is None:
            raise ValueError("A protonated E3 assessment must identify its single evaluated OH orientation")
        probe = output / "e3_probe.pdb"
        changed = replace_hydroxyl_hydrogen(probe_bytes, selected_xyz)
        probe.write_bytes(changed)
        result = copy.deepcopy(original)
        result.update(measured)
        result.update({"probe_pdb": str(probe), "probe_pdb_sha256": digest(probe),
                       "source_e3_artifact": str(original_path), "source_e3_artifact_sha256": original_hash,
                       "source_e3_probe": str(source_probe), "source_e3_probe_sha256": original["probe_pdb_sha256"],
                       "refresh_identity": task["identity"], "reused": False})
        result["provenance"].update({"pdb": str(probe), "pdb_sha256": digest(probe)})
        chemistry = result["chemistry"]
        chemistry.update({"e3_pdb": str(probe), "e3_pdb_sha256": digest(probe),
                          "added_h_xyz_A": selected_xyz,
                          "added_h_placement": "one joint H-bond-aware orientation reselected from the twelve saved OH torsions"})
        for atom in chemistry.get("ligand_atoms", []):
            if atom["name"] == "H25":
                atom["xyz_A"] = selected_xyz
        sampling = result["hydrogen_sampling"]
        sampling["original_steric_selected_orientation_index"] = sampling["selected_orientation_index"]
        sampling.update({"selected_orientation_index": measured["selected_hydrogen_orientation_index"],
                         "added_h_xyz_A": selected_xyz,
                         "selection_basis": "single joint H-bond geometry assessment of all oxygen overlaps"})
        for key in ("selected_nonpolar_overlap_A", "selected_persisted_nonpolar_overlap_A", "selected_persisted_contacts"):
            sampling.pop(key, None)
        # This independent diagnostic is unchanged: same twelve H positions,
        # protein and non-H-bond-compatible atom types as the saved preparation.
        result["summary"]["e3_oh_h_best_sampled_protein_overlap_A"] = original["summary"].get(
            "e3_oh_h_best_sampled_protein_overlap_A")
        write_json(artifact, result)
        (output / "automatic_summary.md").write_text(render_e3_markdown(result), encoding="utf-8")
        return result
    except Exception as error:
        result = {"schema_version": 1, "candidate": task["candidate"], "status": "failed",
                  "error": f"{type(error).__name__}: {error}", "traceback": traceback.format_exc(),
                  "source_e3_artifact": task["artifact"]}
        write_json(artifact, result)
        return result


def batch_report(manifest, results):
    from .e3_reports import render_e3_markdown
    assessments = Counter(r.get("summary", {}).get("e3_steric_assessment") for r in results if r["status"] == "complete")
    lines = ["# E3 steric screen: H-bond geometry checked", "",
             f"Completed {manifest['status_counts'].get('complete', 0):,} mutant assessments; "
             f"{manifest['status_counts'].get('excluded_input_rejection', 0)} prior chemistry exclusions; "
             f"{manifest['status_counts'].get('failed', 0)} new failures.", "",
             "Every large added-oxygen overlap is checked for compatible donor/acceptor chemistry and actual "
             "D-H...A geometry. A single OH-hydrogen orientation must support the complete contact set; "
             "contacts from different orientations are not pooled into an artificial clearance.", "",
             "| Fixed-pose assessment | Mutants |", "|---|---:|"]
    for key in ("steric_overlap_detected", "hbond_supported_overlaps_only", "no_large_overlap_in_fixed_pose", "unassessed_polar_geometry"):
        lines.append(f"| {key.replace('_', ' ')} | {assessments[key]} |")
    lines += ["", "TRUE means at least one large overlap lacks a supported H-bond exemption in the selected joint pose. "
              "FALSE means all overlaps were assessed and none meets the warning rule. An unassessed geometry is "
              "blank/null, with its precise reason retained. No value establishes binding specificity.", "",
              "The geometry gates follow PLIP's published implementation (donor-acceptor distance <4.1 A, "
              "D-H...A angle >100 degrees), with a project safeguard requiring H...A >=1.2 A. "
              "The added-O overlap warning remains the project Bondi-sphere threshold >=0.4 A. "
              "This is a geometry screen, not a native PLIP/DS run or a Rosetta energy calculation.", "",
              "Saved protonated E3 probes and twelve existing OH-H torsions were reused. Only the new OH proton "
              "may be reselected from those stored orientations. Protein and all shared heavy atoms remain fixed. "
              "No OptH, ligand construction, repacking, minimization or rescoring was repeated.", ""]
    lines += [render_e3_markdown(manifest["parent"])]
    return "\n".join(lines)


def run(args):
    from .e3_reports import append_e3_to_payload, e3_columns, flatten_e3_result
    started = time.monotonic()
    source, output = args.source.resolve(), args.output.resolve()
    if source == output or source in output.parents:
        raise ValueError("Use a new corrected output folder outside the saved E3 source")
    output.mkdir(parents=True, exist_ok=True)
    if (output / "e3_batch_manifest.json").exists() and not args.resume:
        raise FileExistsError("Use --resume or a fresh output folder")
    old_payload_path = source / "batch_summary_grouped_fingerprints.json"
    old_payload_hash = digest(old_payload_path)
    old_payload = read_json(old_payload_path)
    original_source = native_path(old_payload["e3_discrimination_protocol"]["source_combined_batch"])
    base_path = original_source / "batch_summary_grouped_fingerprints.json"
    base_hash = digest(base_path)
    base = read_json(base_path)
    identity = {"scope": "fixed_pocket_steric_only", "large_overlap_A": args.overlap_warning,
                "source_e3_payload_sha256": old_payload_hash,
                "code_sha256": {name: digest(Path(__file__).with_name(name))
                                for name in ("e3_steric.py", "e3_hbond.py")}}
    excluded, tasks = [], []
    for row in base["rows"]:
        if row["status"] != "complete":
            excluded.append({"candidate": row["candidate"], "status": "excluded_input_rejection", "error": row.get("error")})
        else:
            stem = Path(row["candidate"]).stem
            tasks.append({"candidate": row["candidate"], "artifact": str(source / "cases" / stem / "e3_discrimination.json"),
                          "output": str(output / "cases" / stem), "identity": identity,
                          "overlap_warning_A": args.overlap_warning, "resume": args.resume})
    parent = refresh_one({"candidate": "parent_m7", "artifact": str(source / "parent/e3_discrimination.json"),
                          "output": str(output / "parent"), "identity": identity,
                          "overlap_warning_A": args.overlap_warning, "resume": args.resume})
    if parent["status"] != "complete":
        raise ValueError(f"Parent refresh failed: {parent['error']}")
    print(json.dumps({"stage": "Hbond_geometry_refresh", "complexes": len(tasks), "workers": args.jobs}), flush=True)
    results = []
    with ProcessPoolExecutor(max_workers=args.jobs) as pool:
        for result in pool.map(refresh_one, tasks, chunksize=8):
            results.append(result)
            if len(results) % 500 == 0:
                print(json.dumps({"completed": len(results), "total": len(tasks)}), flush=True)
    if digest(old_payload_path) != old_payload_hash or digest(base_path) != base_hash:
        raise ValueError("A source report changed during the geometry refresh")
    by_name = {r["candidate"]: r for r in results + excluded}
    counts = dict(Counter(r["status"] for r in results + excluded))
    manifest = {"status": "complete" if not counts.get("failed") else "complete_with_failures",
                "completed_utc": datetime.now(timezone.utc).isoformat(), "candidate_count": len(by_name),
                "source": str(source), "source_payload_sha256": old_payload_hash,
                "protocol": identity, "status_counts": counts, "parent": parent,
                "source_outputs_modified": False, "elapsed_seconds": round(time.monotonic()-started, 2)}
    write_json(output / "e3_batch_manifest.json", manifest)
    (output / "e3_summary.md").write_text(batch_report(manifest, results), encoding="utf-8")
    write_flat_csv(output / "e3_discrimination_summary.csv", {
        "columns": [{"key": "candidate"}, *e3_columns()],
        "rows": [{"candidate": name, **flatten_e3_result(r)} for name, r in by_name.items()]})
    combined = append_e3_to_payload(base, by_name, parent)
    combined["e3_discrimination_protocol"].update({
        "source_combined_batch": str(original_source), "source_payload_sha256": base_hash,
        "source_e3_batch": str(source), "source_e3_payload_sha256": old_payload_hash,
        "measurement_manifest": str(output / "e3_batch_manifest.json"), "overlap_warning_A": args.overlap_warning})
    for old, new in zip([base["parent"], *base["rows"]], [combined["parent"], *combined["rows"]]):
        if any(old.get(c["key"]) != new.get(c["key"]) for c in base["columns"]):
            raise ValueError("A non-E3 metric changed during geometry refresh")
    write_json(output / "batch_summary_grouped_fingerprints.json", combined)
    csv_path = output / "batch_summary_grouped_fingerprints.csv"
    write_flat_csv(csv_path, combined)
    shutil.copyfile(original_source / "candidate_sources.csv", output / "candidate_sources.csv")
    combined_manifest = read_json(original_source / "combined_manifest.json")
    combined_manifest.update({"columns": len(combined["columns"]), "csv": str(csv_path), "csv_sha256": digest(csv_path),
                              "e3_discrimination": str(output / "e3_batch_manifest.json"),
                              "interaction_fingerprints": {"grouped_summary_xlsx": str(output / "batch_summary_grouped_fingerprints.xlsx")}})
    write_json(output / "combined_manifest.json", combined_manifest)
    print(json.dumps({"status_counts": counts, "output": str(output), "elapsed_seconds": manifest["elapsed_seconds"]}), flush=True)
    if counts.get("failed"):
        raise RuntimeError("E3 geometry refresh failed for some cases; detailed errors retained")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--jobs", type=int, default=3)
    parser.add_argument("--overlap-warning", type=float, default=.4)
    parser.add_argument("--resume", action="store_true")
    args = parser.parse_args()
    if args.jobs < 1 or args.overlap_warning <= 0:
        parser.error("Workers and overlap warning must be positive")
    run(args)


if __name__ == "__main__":
    main()
