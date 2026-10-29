"""One-pass whole-complex hydrogen optimization and saved-batch rescoring.

Uses the exact historical RosettaLigand protocol of the source rescore. Original
files and fingerprints are never modified. The optimized PDB is reloaded before
authoritative scoring, so the exported coordinates and scores agree.
"""
from __future__ import annotations

import argparse
from concurrent.futures import ProcessPoolExecutor, as_completed
import csv
from datetime import datetime, timezone
from hashlib import sha256
from importlib import metadata, util
import json
import math
import multiprocessing
from pathlib import Path
from time import perf_counter
import traceback

from .hydrogen_benchmark import _hydrogen_task, _import_audit, _snapshot, compare_snapshots
from .paths import native_path
from .report import write_csv
from .rescore_saved import _read_score, _save

_ENGINE = _SCORE = None
INTERFACE_TERMS = ("fa_atr", "fa_rep", "fa_sol", "fa_elec", "fa_pair", "hbond_bb_sc", "hbond_sc")


def _json(path):
    return json.loads(Path(path).read_text(encoding="utf-8-sig"))


def _digest(path):
    return sha256(Path(path).read_bytes()).hexdigest()


def read_source(source: Path, expected_count: int | None = None):
    manifest = _json(source / "rescore_manifest.json")
    if manifest.get("status") != "complete" or manifest.get("failure_count") != 0:
        raise ValueError("Source must be a complete, failure-free saved rescore")
    protocol = manifest["protocol"]
    if protocol.get("weights") != "ligand" or not protocol.get("legacy_pre_talaris"):
        raise ValueError("Expected the original historical ligand/pre-Talaris scoring protocol")
    with (source / "rosettaligand_rescore.csv").open(newline="", encoding="utf-8-sig") as handle:
        rows = list(csv.DictReader(handle))
    names = [row["candidate"] for row in rows]
    if len(set(names)) != len(names):
        raise ValueError("Duplicate candidate names in source")
    if len(rows) != manifest["completed_count"] or (expected_count is not None and len(rows) != expected_count):
        raise ValueError("Source candidate count differs from expected completed count")
    return manifest, rows


def reuse_completed(path: Path, identity: dict):
    if not path.exists():
        return None
    result = _json(path)
    if result.get("status") != "complete" or result.get("identity") != identity:
        raise ValueError(f"Refusing incompatible cached result: {path}")
    if not result.get("artifact_sha256"):
        raise ValueError(f"Cached result has no artifact hashes: {path}")
    for artifact, digest in result["artifact_sha256"].items():
        if not Path(artifact).is_file() or _digest(artifact) != digest:
            raise ValueError(f"Cached artifact is missing or modified: {artifact}")
    return result


def make_summary_row(candidate, substitutions, before, after, parent_before,
                     parent_after, optimized_pdb, optimization):
    key = "interface_delta_X"
    row = {
        "candidate": candidate, "substitutions": substitutions, "status": "complete",
        "optimized_pdb": str(optimized_pdb), "weights": "ligand",
        "interface_parent_REU": parent_after[key], "interface_mutant_REU": after[key],
        "interface_delta_REU": after[key] - parent_after[key],
        "before_interface_parent_REU": parent_before[key],
        "before_interface_mutant_REU": before[key],
        "before_interface_delta_REU": before[key] - parent_before[key],
        "hydrogen_optimization_interface_change_REU": after[key] - before[key],
        "parent_bound_total_REU": parent_after["bound_total"],
        "mutant_bound_total_REU": after["bound_total"],
        "before_parent_bound_total_REU": parent_before["bound_total"],
        "before_mutant_bound_total_REU": before["bound_total"],
    }
    for term in INTERFACE_TERMS:
        metric = "if_X_" + term
        row.update({f"parent_{term}_REU": parent_after[metric],
                    f"mutant_{term}_REU": after[metric],
                    f"delta_{term}_REU": after[metric] - parent_after[metric]})
    row.update({"selected_trial": optimization["selected_trial"],
                "optimization_seconds": optimization["optimization_seconds"],
                "max_heavy_displacement_A": optimization["max_heavy_displacement_A"],
                "moved_hydrogen_count": optimization["moved_hydrogen_count"]})
    return row


def _init_worker(params, weights="ligand", legacy_pre_talaris=None):
    global _ENGINE, _SCORE
    from .pyrosetta_score import _initialize, _scorefunction, resolve_legacy_pre_talaris
    _ENGINE = _initialize([params], legacy=resolve_legacy_pre_talaris(weights, legacy_pre_talaris))
    _ENGINE.rosetta.basic.options.set_boolean_option("packing:no_optH", True)
    _ENGINE.rosetta.basic.options.set_boolean_option("packing:flip_HNQ", False)
    _SCORE = _scorefunction(_ENGINE, weights)


def _score_metrics(pose, chain, output):
    from .pyrosetta_score import write_interface_metrics
    write_interface_metrics(pose, _SCORE, chain, output, Path(output).stem)
    values = _read_score(output)
    if not all(math.isfinite(value) for value in values.values()):
        raise ValueError("Non-finite score encountered")
    if not math.isclose(sum(v for k, v in values.items() if k.startswith(f"if_{chain}_")),
                        values[f"interface_delta_{chain}"], abs_tol=5e-5):
        raise ValueError("Interface decomposition does not reconstruct the total")
    return values


def _validate_import(pdb, pose):
    audit = _import_audit(pdb, pose)
    if audit["input_atoms_absent_from_pose"] or audit["max_shared_atom_displacement_A"] > 1e-6:
        raise ValueError(f"Input atoms were discarded or moved during import: {audit}")
    info = pose.pdb_info()
    virtual = {f"{info.chain(i)}:{info.number(i)}:{pose.residue(i).atom_name(j).strip()}"
               for i in range(1, pose.size() + 1)
               for j in range(1, pose.residue(i).natoms() + 1)
               if pose.residue(i).atom_type(j).is_virtual()}
    if set(audit["pose_atoms_absent_from_input"]) - virtual:
        raise ValueError(f"Unexpected nonvirtual atoms added on import: {audit}")
    return audit


def _run_one(record_path, output, pdb_output, protocol):
    started = perf_counter()
    record_path, output, pdb_output = Path(record_path), Path(output), Path(pdb_output)
    source_record = _json(record_path)
    pdb = native_path(source_record["pdb"])
    if source_record["identity"]["protocol"] != protocol["source_protocol"]:
        raise ValueError("Source case scoring protocol differs from batch")
    if _digest(pdb) != source_record["identity"]["pdb_sha256"]:
        raise ValueError(f"Original prepared PDB was changed: {pdb}")
    identity = {"protocol": protocol, "pdb_sha256": source_record["identity"]["pdb_sha256"]}
    cached = reuse_completed(output / "result.json", identity)
    if cached is not None:
        return cached
    output.mkdir(parents=True, exist_ok=True)
    pose = _ENGINE.pose_from_file(str(pdb))
    import_audit = _validate_import(pdb, pose)
    snapshot = _snapshot(pose)
    if not all(math.isfinite(v) for atom in snapshot["atoms"].values() for v in atom["xyz"]):
        raise ValueError("Non-finite input coordinates")
    chain = protocol["source_protocol"]["ligand_chain"]
    before = _score_metrics(pose, chain, output / "before_interface.sc")
    for key, value in source_record["metrics"].items():
        if not math.isclose(before[key], value, abs_tol=1e-4):
            raise ValueError(f"Saved baseline score differs for {key}: {before[key]} vs {value}")
    initial_pose = pose.clone()
    _ENGINE.rosetta.numeric.random.rg().set_seed(protocol["seed"])
    opt_started = perf_counter()
    task, supported = _hydrogen_task(_ENGINE, pose)
    _SCORE(pose)
    _ENGINE.rosetta.core.pack.pack_rotamers(pose, _SCORE, task)
    optimization_seconds = perf_counter() - opt_started
    trial_snapshot = _snapshot(pose)
    if not all(math.isfinite(v) for atom in trial_snapshot["atoms"].values() for v in atom["xyz"]):
        raise ValueError("Non-finite optimized coordinates")
    trial_audit = compare_snapshots(snapshot, trial_snapshot)
    trial_bound = float(_SCORE(pose))
    if not math.isfinite(trial_bound):
        raise ValueError("Non-finite optimized bound score")
    selected_trial = 1 if trial_bound < before["bound_total"] else 0
    if not selected_trial:
        pose = initial_pose
    selected_audit = compare_snapshots(snapshot, _snapshot(pose))
    pdb_output.parent.mkdir(parents=True, exist_ok=True)
    temporary = pdb_output.with_name(pdb_output.stem + ".tmp.pdb")
    pose.dump_pdb(str(temporary))
    temporary.replace(pdb_output)
    # Score exactly the persisted coordinates, including PDB rounding of H atoms.
    persisted = _ENGINE.pose_from_file(str(pdb_output))
    persisted_audit = compare_snapshots(snapshot, _snapshot(persisted))
    after = _score_metrics(persisted, chain, output / "interface.sc")
    from .hbonds import write_pose_hbonds
    write_pose_hbonds(persisted, output / "hbonds.csv", scorefxn=_SCORE)
    with (output / "hbonds.csv").open(newline="") as handle:
        hbonds = list(csv.DictReader(handle))
    ligand_hbonds = [hb for hb in hbonds if hb["donor_residue"].startswith(chain + ":")
                     != hb["acceptor_residue"].startswith(chain + ":")]
    weighted_hb = sum(float(hb["hbond_weighted_energy_REU"]) for hb in hbonds)
    from .pyrosetta_score import _active_terms, _weighted_totals
    bound_terms = _weighted_totals(persisted, _SCORE, _active_terms(_SCORE))
    if not math.isclose(weighted_hb, sum(v for k, v in bound_terms.items() if k.startswith("hbond_")), abs_tol=1e-4):
        raise ValueError("H-bond records do not reconstruct the scored energy")
    optimization = {"selected_trial": selected_trial, "optimization_seconds": optimization_seconds,
                    "trial_bound_total_REU": trial_bound,
                    "max_heavy_displacement_A": persisted_audit["max_heavy_displacement_A"],
                    "moved_hydrogen_count": persisted_audit["moved_hydrogen_count"],
                    "selected_audit_before_serialization": selected_audit,
                    "persisted_audit": persisted_audit, "trial_audit": trial_audit,
                    "supported_residues": supported}
    result = {"status": "complete", "identity": identity, "source_result": str(record_path),
              "input_pdb": str(pdb), "optimized_pdb": str(pdb_output), "before_metrics": before,
              "metrics": after, "before_ligand_hbond_count": source_record["ligand_hbond_count"],
              "ligand_hbond_count": len(ligand_hbonds), "optimization": optimization,
              "import_audit": import_audit, "elapsed_seconds": perf_counter() - started,
              "artifact_sha256": {str(path): _digest(path) for path in
                                   (pdb_output, output / "before_interface.sc", output / "interface.sc", output / "hbonds.csv")}}
    _save(output / "result.json", result)
    return result


def run(source: Path, output: Path, jobs=3, expected_count=None, seed=20261022):
    started = perf_counter()
    source, output = source.resolve(), output.resolve()
    if jobs < 1 or source == output or source in output.parents:
        raise ValueError("Use positive jobs and a separate sibling output directory")
    manifest, rows = read_source(source, expected_count)
    source_batch = native_path(manifest["source_batch"])
    first_workflow = _json(native_path(rows[0]["source_result_folder"]) / "workflow_manifest.json")
    params = native_path(first_workflow["params"]).resolve()
    if _digest(params) != manifest["protocol"]["params_sha256"]:
        raise ValueError("Ligand parameter hash differs from source")
    if metadata.version("pyrosetta") != manifest["protocol"]["pyrosetta_version"]:
        raise ValueError("PyRosetta version differs from source")
    database = Path(util.find_spec("pyrosetta").origin).parent / "database"
    weights = database / "scoring/weights/ligand.wts"
    if _digest(weights) != manifest["protocol"]["weights_sha256"]:
        raise ValueError("Installed ligand.wts differs from source")
    protocol = {
        "version": 1, "source_protocol": manifest["protocol"], "seed": seed,
        "passes_per_complex": 1, "selection": "lowest bound energy of original and one optimized trial",
        "optimization": "whole-complex supported proton chis; fixed heavy atoms and chemical states",
        "score_persisted_pdb": True,
        "implementation_sha256": {name: _digest(Path(__file__).with_name(name)) for name in
                                   ("hydrogen_rescore.py", "hydrogen_benchmark.py", "pyrosetta_score.py", "hbonds.py")},
    }
    output.mkdir(parents=True, exist_ok=True)
    manifest_path = output / "hydrogen_rescore_manifest.json"
    if manifest_path.exists() and _json(manifest_path)["protocol"] != protocol:
        raise ValueError("Existing output belongs to a different optimization protocol")
    progress = {"status": "running", "source_rescore": str(source), "source_batch": str(source_batch),
                "candidate_count": len(rows), "jobs": jobs, "protocol": protocol,
                "started_utc": datetime.now(timezone.utc).isoformat(), "completed_count": 0, "failure_count": 0,
                "fingerprints_regenerated": False}
    _save(manifest_path, progress)
    structures = output / "hydrogen_optimized_structures"
    results, failures = {}, []
    context = multiprocessing.get_context("spawn")
    with ProcessPoolExecutor(max_workers=jobs, mp_context=context,
                             initializer=_init_worker, initargs=(str(params),)) as pool:
        parent = pool.submit(_run_one, str(source / "parent/result.json"), str(output / "parent"),
                             str(structures / "parent_m7.pdb"), protocol).result()
        print(json.dumps({"parent": "complete", "interface_REU": parent["metrics"]["interface_delta_X"]}), flush=True)
        # Validate one mutant end-to-end before dispatching the full set. Its
        # completed result is part of the final batch, not discarded smoke work.
        first = rows[0]
        first_name = first["candidate"]
        results[first_name] = pool.submit(
            _run_one, str(source / Path(first_name).stem / "result.json"),
            str(output / "cases" / Path(first_name).stem), str(structures / first_name), protocol,
        ).result()
        futures = {pool.submit(_run_one, str(source / Path(row["candidate"]).stem / "result.json"),
                               str(output / "cases" / Path(row["candidate"]).stem),
                               str(structures / row["candidate"]), protocol): row for row in rows[1:]}
        with (output / "progress.jsonl").open("a", encoding="utf-8") as log:
            log.write(json.dumps({"candidate": first_name, "status": "complete", "preflight": True}) + "\n")
            log.flush()
            for future in as_completed(futures):
                row = futures[future]
                try:
                    results[row["candidate"]] = future.result()
                    event = {"candidate": row["candidate"], "status": "complete"}
                except Exception as error:
                    event = {"candidate": row["candidate"], "status": "failed", "error": str(error),
                             "traceback": traceback.format_exc()}
                    failures.append(event)
                log.write(json.dumps(event) + "\n")
                log.flush()
                count = len(results) + len(failures)
                if count % 100 == 0 or count == len(rows):
                    progress.update(completed_count=len(results), failure_count=len(failures),
                                    elapsed_seconds=perf_counter() - started)
                    _save(manifest_path, progress)
                    print(json.dumps({"processed": count, "completed": len(results), "failed": len(failures)}), flush=True)
    table = []
    for source_row in rows:
        name = source_row["candidate"]
        if name not in results:
            table.append({"candidate": name, "substitutions": source_row.get("substitutions", ""), "status": "failed"})
            continue
        result = results[name]
        row = make_summary_row(name, source_row.get("substitutions", ""), result["before_metrics"], result["metrics"],
                               parent["before_metrics"], parent["metrics"], result["optimized_pdb"], result["optimization"])
        row.update(parent_ligand_hbonds=parent["ligand_hbond_count"], mutant_ligand_hbonds=result["ligand_hbond_count"],
                   before_parent_ligand_hbonds=parent["before_ligand_hbond_count"],
                   before_mutant_ligand_hbonds=result["before_ligand_hbond_count"],
                   parent_optimized_pdb=parent["optimized_pdb"],
                   optimized_pdb_sha256=result["artifact_sha256"][result["optimized_pdb"]])
        table.append(row)
    table.sort(key=lambda row: (row["status"] != "complete", row.get("interface_delta_REU", math.inf)))
    for rank, row in enumerate(table, 1):
        row["rank_by_interface_delta"] = rank if row["status"] == "complete" else ""
    csv_path = output / "rosettaligand_hydrogen_optimized_rescore.csv"
    write_csv(csv_path, table)
    _save(output / "failures.json", failures)
    progress.update(status="complete_with_failures" if failures else "complete", completed_count=len(results),
                    failure_count=len(failures), summary_csv=str(csv_path), structures_folder=str(structures),
                    parent_optimized_pdb=parent["optimized_pdb"], elapsed_seconds=perf_counter() - started,
                    finished_utc=datetime.now(timezone.utc).isoformat())
    _save(manifest_path, progress)
    print(json.dumps({key: progress[key] for key in
                      ("status", "completed_count", "failure_count", "summary_csv", "elapsed_seconds")}), flush=True)
    return 1 if failures else 0


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("source", type=Path)
    parser.add_argument("output", type=Path)
    parser.add_argument("--jobs", type=int, default=3)
    parser.add_argument("--expected-count", type=int)
    parser.add_argument("--seed", type=int, default=20261022)
    args = parser.parse_args()
    try:
        return run(args.source, args.output, args.jobs, args.expected_count, args.seed)
    except Exception as error:
        path = args.output / "hydrogen_rescore_manifest.json"
        if path.exists():
            progress = _json(path)
            owned = (progress.get("source_rescore") == str(args.source.resolve())
                     and progress.get("protocol", {}).get("seed") == args.seed
                     and progress.get("protocol", {}).get("implementation_sha256", {}).get("hydrogen_rescore.py")
                     == _digest(__file__))
            if owned:
                progress.update(status="failed", error=str(error), traceback=traceback.format_exc())
                _save(path, progress)
        raise


if __name__ == "__main__":
    raise SystemExit(main())
