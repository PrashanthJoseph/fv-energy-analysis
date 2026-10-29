"""Fixed-pose estriol 16alpha-OH steric screen, separate from E2 recognition.

Read existing prepared complexes and append an optional report section without
rescoring, hydrogen optimization, repacking, minimization or moving input atoms.
"""
from __future__ import annotations

import argparse
import copy
from concurrent.futures import ProcessPoolExecutor
from collections import Counter
from datetime import datetime, timezone
from hashlib import sha256
import csv
import json
from pathlib import Path
import shutil
import time
import traceback

import numpy as np

from .paths import native_path
from .e3_io import replace_hydroxyl_hydrogen
from .report import write_json


def digest(path):
    return sha256(Path(path).read_bytes()).hexdigest()


def read_json(path):
    return json.loads(Path(path).read_text(encoding="utf-8-sig"))


def resolve_saved_complexes(source):
    """Resolve the actual saved analysis PDB, including architecture refreshes."""
    source = Path(source)
    payload_path = source / "batch_summary_grouped_fingerprints.json"
    payload = read_json(payload_path)
    with (source / "candidate_sources.csv").open(encoding="utf-8-sig", newline="") as stream:
        origins = list(csv.DictReader(stream))
    lookup = {r["candidate"]: native_path(r["source"]) for r in origins}
    if len(lookup) != len(origins) or set(lookup) != {r["candidate"] for r in payload["rows"]}:
        raise ValueError("Candidate provenance must uniquely cover the combined dataset")
    tasks, excluded, parent = [], [], None
    for row in payload["rows"]:
        if row["status"] != "complete":
            excluded.append({"candidate": row["candidate"], "status": "excluded_input_rejection",
                             "error": row.get("error"), "analysis_scope": "fixed_pose_steric_only"})
            continue
        folder = lookup[row["candidate"]] / Path(row["candidate"]).stem
        saved = read_json(folder / "estrogen_binding_architecture.json")
        architecture = saved.get("result", saved.get("mutant"))
        if not architecture or architecture.get("e2_architecture_status") != "complete":
            raise ValueError(f"No completed input-associated architecture: {folder}")
        provenance = architecture["e2_provenance"]
        tasks.append({"candidate": row["candidate"], "pdb": str(native_path(provenance["pdb"])),
                      "pdb_sha256": provenance["pdb_sha256"], "source_result": str(folder)})
        if parent is None and saved.get("parent"):
            p = saved["parent"]["e2_provenance"]
            parent = {"candidate": "parent_m7", "pdb": str(native_path(p["pdb"])),
                      "pdb_sha256": p["pdb_sha256"], "source_result": str(folder)}
    if parent is None:
        raise ValueError("No saved parent structure found")
    return payload, tasks, excluded, parent


def measure(task):
    from .e3_chemistry import build_e3_probe
    from .e3_steric import analyze_e3_steric, orient_e3_hydroxyl
    from .e3_reports import render_e3_markdown
    output = Path(task["output"])
    output.mkdir(parents=True, exist_ok=True)
    artifact = output / "e3_discrimination.json"
    try:
        before = digest(task["pdb"])
        if before != task["pdb_sha256"]:
            raise ValueError("Input PDB hash differs from the recorded completed analysis")
        if task.get("resume") and artifact.is_file():
            previous = read_json(artifact)
            if (previous.get("status") == "complete"
                    and previous.get("input_pdb_sha256") == before
                    and previous.get("protocol_identity") == task["protocol_identity"]
                    and Path(previous.get("probe_pdb", "")).is_file()
                    and digest(previous["probe_pdb"]) == previous.get("probe_pdb_sha256")):
                previous["reused"] = True
                (output / "automatic_summary.md").write_text(render_e3_markdown(previous), encoding="utf-8")
                return previous
        probe = output / "e3_probe.pdb"
        chemistry = build_e3_probe(task["pdb"], task["params"], task["template"], output_pdb=probe)
        hydrogen_sampling = orient_e3_hydroxyl(
            probe, added_o_atom=chemistry["added_o_atom"], added_h_atom=chemistry["added_h_atom"],
            c16_atom=chemistry["c16_atom"], n_orientations=12,
        )
        chemistry["initial_probe_sha256"] = chemistry["e3_pdb_sha256"]
        chemistry["e3_pdb_sha256"] = digest(probe)
        chemistry["initial_added_h_xyz_A"] = chemistry["added_h_xyz_A"]
        chemistry["added_h_xyz_A"] = hydrogen_sampling["added_h_xyz_A"]
        chemistry["added_h_placement"] = "selected from twelve fixed-bond steric torsion samples"
        for atom in chemistry.get("ligand_atoms", []):
            if atom["name"] == chemistry["added_h_atom"]:
                atom["xyz_A"] = hydrogen_sampling["added_h_xyz_A"]
        samples = copy.deepcopy(hydrogen_sampling["orientations"])
        for sample in samples:
            sample["hydrogen_xyz_A"] = np.round(sample["hydrogen_xyz_A"], 3).tolist()
        result = analyze_e3_steric(probe, added_o_atom=chemistry["added_o_atom"],
                                   overlap_warning_A=task["overlap_warning_A"], hydrogen_positions=samples)
        selected_xyz = result["selected_hydrogen_xyz_A"]
        probe.write_bytes(replace_hydroxyl_hydrogen(probe.read_bytes(), selected_xyz))
        chemistry["e3_pdb_sha256"] = digest(probe)
        chemistry["added_h_xyz_A"] = selected_xyz
        chemistry["added_h_placement"] = "joint H-bond geometry selection among twelve saved fixed-bond OH torsions"
        for atom in chemistry.get("ligand_atoms", []):
            if atom["name"] == chemistry["added_h_atom"]:
                atom["xyz_A"] = selected_xyz
        hydrogen_sampling["original_steric_selected_orientation_index"] = hydrogen_sampling["selected_orientation_index"]
        hydrogen_sampling["selected_orientation_index"] = result["selected_hydrogen_orientation_index"]
        hydrogen_sampling["added_h_xyz_A"] = selected_xyz
        for key in ("selected_persisted_contacts", "selected_nonpolar_overlap_A", "selected_persisted_nonpolar_overlap_A"):
            hydrogen_sampling.pop(key, None)
        result["provenance"].update({"pdb": str(probe), "pdb_sha256": digest(probe)})
        result["hydrogen_sampling"] = hydrogen_sampling
        result["summary"]["e3_oh_h_best_sampled_protein_overlap_A"] = hydrogen_sampling["minimum_protein_nonpolar_overlap_A"]
        if digest(task["pdb"]) != before:
            raise ValueError("Input structure changed during analysis")
        result.update({"candidate": task["candidate"], "status": "complete", "chemistry": chemistry,
                       "input_pdb": task["pdb"], "input_pdb_sha256": before,
                       "source_result": task.get("source_result"),
                       "probe_pdb": str(probe), "probe_pdb_sha256": digest(probe),
                       "protocol_identity": task["protocol_identity"], "reused": False})
        write_json(artifact, result)
        (output / "automatic_summary.md").write_text(render_e3_markdown(result), encoding="utf-8")
        return result
    except Exception as error:
        result = {"candidate": task["candidate"], "status": "failed",
                  "analysis_scope": "fixed_pose_steric_only", "input_pdb": task["pdb"],
                  "error_type": type(error).__name__, "error": str(error),
                  "traceback": traceback.format_exc(), "protocol_identity": task["protocol_identity"]}
        write_json(artifact, result)
        return result


def write_flat_csv(path, payload):
    with Path(path).open("w", encoding="utf-8-sig", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=[c["key"] for c in payload["columns"]], extrasaction="ignore")
        writer.writeheader()
        writer.writerows(payload["rows"])


def render_batch_report(manifest, results):
    from .e3_reports import render_e3_markdown
    valid = [r for r in results if r["status"] == "complete"]
    assessments = Counter(r["summary"]["e3_steric_assessment"] for r in valid)
    lines = ["# E3 fixed-pose steric screen", "",
             f"Successfully assessed {len(valid):,} mutant poses. "
             f"Prior chemistry exclusions: {manifest['status_counts'].get('excluded_input_rejection', 0)}. "
             f"E3 assessment failures: {manifest['status_counts'].get('failed', 0)}.", "",
             "Every E3 probe contains the correctly stereochemical additional 16alpha-OH and all 24 ligand hydrogens. "
             "Protein coordinates and the twenty shared steroid heavy-atom coordinates remain fixed. "
             "Twelve orientations of the new OH hydrogen are sampled by rotation only. "
             "No Rosetta scoring, hydrogen-energy optimization, repacking or minimization is performed.", "",
             "| Fixed-pose assessment | Mutants |", "|---|---:|"]
    for key in ("steric_overlap_detected", "hbond_supported_overlaps_only", "unassessed_polar_geometry", "no_large_overlap_in_fixed_pose"):
        lines.append(f"| {key.replace('_', ' ')} | {assessments[key]} |")
    lines += ["", "The oxygen warning uses a project heavy-atom overlap threshold of "
              f"{manifest['protocol']['overlap_warning_A']:.2f} A. Each polar overlap receives a chemical-role and "
              "explicit D-H...A geometry check before it can be exempted. A single OH-H orientation must account "
              "for all simultaneous contacts. Missing geometry remains unassessed rather than FALSE. "
              "This borrows the idea of non-hydrogen-bond "
              "overlap from MolProbity; it is not the PROBE all-atom algorithm or a validated estriol exclusion cutoff. "
              "A warning indicates an obstruction in this specific pose. It does not prove E3 rejection; "
              "a warning-free pose does not establish E3 binding.", "",
              "The [SHBG structural study](https://pmc.ncbi.nlm.nih.gov/articles/PMC305588/) proposes "
              "steric interference from estriol's 16alpha-OH near a peptide bond. Here the corresponding "
              "protein atoms are measured directly rather than assumed to be the same residues in m7.", "",
              "Original E2 scores, ranks and detector metrics are retained unchanged. "
              "Per-case probe PDB, contact JSON and report are in `cases/`; the parent's assessment is in `parent/`.", ""]
    if manifest.get("parent"):
        lines += [render_e3_markdown(manifest["parent"])]
    return "\n".join(lines)


def run(args):
    from .e3_reports import append_e3_to_payload, flatten_e3_result, e3_columns
    started = time.monotonic()
    output = args.output.resolve()
    output.mkdir(parents=True, exist_ok=True)
    if (output / "e3_batch_manifest.json").exists() and not args.resume:
        raise FileExistsError("Choose a new output folder or use --resume")
    params, template = args.params.resolve(), args.template.resolve()
    protocol = {
        "scope": "fixed_pose_steric_only", "extra_substituent": "16alpha-OH",
        "params_sha256": digest(params), "template_sha256": digest(template),
        "overlap_warning_A": args.overlap_warning,
        "shared_steroid_and_protein_atoms_fixed": True, "minimization": False,
        "repacking": False, "hydrogen_optimization": False, "energetics_assessed": False,
        "ligand_hydrogens": "all 24 explicit; neutral estriol", "new_oh_h_torsion_samples": 12,
        "implementation_sha256": {name: digest(Path(__file__).with_name(name)) for name in
                                  ("e3_chemistry.py", "e3_steric.py", "e3_hbond.py")},
    }
    payload, excluded, parent = None, [], None
    source = args.source.resolve() if args.source else None
    if source:
        if source == output or source in output.parents:
            raise ValueError("Write E3 results outside the completed source batch")
        payload, tasks, excluded, parent = resolve_saved_complexes(source)
        scientific = read_json(source / "combined_manifest.json").get("scientific_protocol", {})
        for key in ("params_sha256", "template_sha256"):
            if scientific.get(key) and scientific[key] != protocol[key]:
                raise ValueError(f"E3 chemical input differs from the saved scientific protocol: {key}")
        input_hash = digest(source / "batch_summary_grouped_fingerprints.json")
    else:
        inputs = [args.pdb.resolve()] if args.pdb else sorted(args.input_folder.resolve().glob(args.glob))
        if not inputs:
            raise ValueError("No input PDBs found")
        if len({p.name for p in inputs}) != len(inputs):
            raise ValueError("Input PDB basenames must be unique")
        tasks = [{"candidate": p.name, "pdb": str(p), "pdb_sha256": digest(p)} for p in inputs]
        if args.parent:
            p = args.parent.resolve()
            parent = {"candidate": "parent_m7", "pdb": str(p), "pdb_sha256": digest(p)}
        input_hash = None
    if args.limit is not None:
        tasks = tasks[:args.limit]
        if payload:
            selected = {r["candidate"] for r in tasks}
            excluded.extend({"candidate": r["candidate"], "status": "not_assessed",
                             "analysis_scope": "fixed_pose_steric_only"}
                            for r in payload["rows"] if r["status"] == "complete" and r["candidate"] not in selected)
    for task in tasks + ([parent] if parent else []):
        task.update({"params": str(params), "template": str(template),
                     "overlap_warning_A": args.overlap_warning, "protocol_identity": protocol,
                     "resume": args.resume, "output": str(output / "cases" / Path(task["candidate"]).stem)})
    if parent:
        parent["output"] = str(output / "parent")
    parent_result = measure(parent) if parent else None
    if parent_result and parent_result["status"] != "complete":
        raise ValueError(f"Parent E3 probe failed: {parent_result['error']}")
    print(json.dumps({"stage": "fixed_pose_E3_sterics", "complexes": len(tasks), "workers": args.jobs}), flush=True)
    results = []
    with ProcessPoolExecutor(max_workers=args.jobs) as pool:
        for result in pool.map(measure, tasks, chunksize=8):
            results.append(result)
            if len(results) % 250 == 0:
                print(json.dumps({"completed": len(results), "total": len(tasks)}), flush=True)
    by_name = {r["candidate"]: r for r in results + excluded}
    counts = dict(Counter(r["status"] for r in results + excluded))
    manifest = {"status": "complete" if not counts.get("failed") else "complete_with_failures",
                "completed_utc": datetime.now(timezone.utc).isoformat(), "protocol": protocol,
                "source": str(source) if source else None, "source_payload_sha256": input_hash,
                "candidate_count": len(results) + len(excluded), "status_counts": counts,
                "elapsed_seconds": round(time.monotonic()-started, 2),
                "parent": parent_result, "source_outputs_modified": False}
    write_json(output / "e3_batch_manifest.json", manifest)
    (output / "e3_summary.md").write_text(render_batch_report(manifest, results), encoding="utf-8")
    rows = [{"candidate": name, **flatten_e3_result(result)} for name, result in by_name.items()]
    columns = [{"key": "candidate"}, *e3_columns()]
    write_flat_csv(output / "e3_discrimination_summary.csv", {"columns": columns, "rows": rows})
    if payload:
        if digest(source / "batch_summary_grouped_fingerprints.json") != input_hash:
            raise ValueError("Source combined dataset changed during the E3 calculation")
        combined = append_e3_to_payload(payload, by_name, parent_result)
        combined["e3_discrimination_protocol"].update({
            "source_combined_batch": str(source), "source_payload_sha256": input_hash,
            "measurement_manifest": str(output / "e3_batch_manifest.json"),
            "overlap_warning_A": args.overlap_warning,
        })
        keys = [c["key"] for c in payload["columns"]]
        assert len(combined["rows"]) == len(payload["rows"])
        for before, after in zip(payload["rows"], combined["rows"]):
            if any(before.get(key) != after.get(key) for key in keys):
                raise ValueError(f"Unrelated saved metric changed: {before['candidate']}")
        write_json(output / "batch_summary_grouped_fingerprints.json", combined)
        csv_path = output / "batch_summary_grouped_fingerprints.csv"
        write_flat_csv(csv_path, combined)
        shutil.copyfile(source / "candidate_sources.csv", output / "candidate_sources.csv")
        combined_manifest = read_json(source / "combined_manifest.json")
        combined_manifest.update({"columns": len(combined["columns"]), "csv": str(csv_path),
                                  "csv_sha256": digest(csv_path), "e3_discrimination": str(output / "e3_batch_manifest.json"),
                                  "interaction_fingerprints": {"grouped_summary_xlsx": str(output / "batch_summary_grouped_fingerprints.xlsx")}})
        write_json(output / "combined_manifest.json", combined_manifest)
    print(json.dumps({"output": str(output), "status_counts": counts,
                      "elapsed_seconds": manifest["elapsed_seconds"]}), flush=True)
    if counts.get("failed"):
        raise RuntimeError("Some E3 probes failed; see their detailed saved errors")
    return manifest


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    inputs = parser.add_mutually_exclusive_group(required=True)
    inputs.add_argument("--source", type=Path, help="Completed combined dataset with candidate_sources.csv")
    inputs.add_argument("--pdb", type=Path)
    inputs.add_argument("--input-folder", type=Path)
    parser.add_argument("--parent", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--params", type=Path, required=True)
    parser.add_argument("--template", type=Path, required=True)
    parser.add_argument("--glob", default="*.pdb")
    parser.add_argument("--jobs", type=int, default=3)
    parser.add_argument("--overlap-warning", type=float, default=.4)
    parser.add_argument("--limit", type=int)
    parser.add_argument("--resume", action="store_true")
    args = parser.parse_args()
    if args.jobs < 1 or args.overlap_warning <= 0 or args.limit is not None and args.limit < 1:
        parser.error("Worker count, overlap warning and optional limit must be positive")
    run(args)


if __name__ == "__main__":
    main()
