"""Prepare a fresh relaxed folder and optimize hydrogens once before analysis.

Reuses the validated preparation and whole-complex proton-chi implementation.
No interaction detector is run here. Heavy atoms and chemical states are fixed;
the selected, persisted PDBs are the inputs to the normal full batch module.
"""
from __future__ import annotations

import argparse
from concurrent.futures import ProcessPoolExecutor, as_completed
from datetime import datetime, timezone
from hashlib import sha256
from importlib import metadata, util
import json
import math
import multiprocessing
from pathlib import Path
import shutil
from time import perf_counter
import traceback

from .hydrogen_benchmark import _hydrogen_task, _snapshot, compare_snapshots
from .hydrogen_rescore import _score_metrics, _validate_import, reuse_completed
from .paths import native_path
from .rescore_saved import _save
from .resources import resource_path

_ENGINE = _SCORE = _PREPARER = None


def digest(path):
    return sha256(Path(path).read_bytes()).hexdigest()


def read_json(path):
    return json.loads(Path(path).read_text(encoding="utf-8-sig"))


def list_inputs(folder: Path, parent: Path, recursive=False, pattern="*.pdb") -> list[Path]:
    files = sorted((folder.rglob(pattern) if recursive else folder.glob(pattern)))
    files = [path.resolve() for path in files if path.resolve() != parent.resolve()]
    if not files:
        raise ValueError("Input folder contains no PDB complexes")
    names = [path.name for path in files]
    if len(set(names)) != len(names) or "parent_m7.pdb" in names:
        raise ValueError("Candidate basenames must be unique and must not be parent_m7.pdb")
    return files


def reuse_case(path: Path, identity: dict):
    """Bind resume to raw input/settings and verify every saved output hash."""
    if not path.exists():
        return None
    record = read_json(path)
    if any(record.get("identity", {}).get(key) != value for key, value in identity.items()):
        raise ValueError(f"Refusing incompatible cached result: {path}")
    return reuse_completed(path, record["identity"])


def validate_reusable_parent(record: dict, protocol: dict, prepared: Path):
    """Reuse only the same prepared parent and the exact one-pass protocol."""
    previous = record["identity"]["protocol"]
    old, new = previous["source_protocol"], protocol["source_protocol"]
    keys = ("weights", "legacy_pre_talaris", "params_sha256", "weights_sha256",
            "pyrosetta_version", "ligand_chain")
    if (record.get("status") != "complete" or any(old.get(key) != new.get(key) for key in keys)
            or previous.get("passes_per_complex") != 1
            or previous.get("seed") != protocol["seed"]
            or previous.get("selection") != protocol["selection"]
            or record["identity"].get("pdb_sha256") != digest(prepared)):
        raise ValueError("Existing parent OptH result has different inputs or settings")
    implementations = previous.get("implementation_sha256", {})
    for name in ("hydrogen_benchmark.py", "pyrosetta_score.py", "hbonds.py"):
        if implementations.get(name) != digest(Path(__file__).with_name(name)):
            raise ValueError(f"Existing parent OptH implementation differs: {name}")
    for artifact, expected in record.get("artifact_sha256", {}).items():
        if digest(native_path(artifact)) != expected:
            raise ValueError(f"Existing parent OptH artifact changed: {artifact}")
    optimized = record["optimized_pdb"]
    if record.get("artifact_sha256", {}).get(optimized) != digest(native_path(optimized)):
        raise ValueError("Existing parent optimized PDB is not hash-validated")


def _init_worker(params, weights="ligand", legacy_pre_talaris=None):
    global _ENGINE, _SCORE, _PREPARER
    from . import hydrogen_rescore
    hydrogen_rescore._init_worker(params, weights, legacy_pre_talaris)
    _ENGINE, _SCORE = hydrogen_rescore._ENGINE, hydrogen_rescore._SCORE
    spec = util.spec_from_file_location("fv_prepare_ligand", resource_path("scripts/prepare_ligand_complexes.py"))
    _PREPARER = util.module_from_spec(spec)
    spec.loader.exec_module(_PREPARER)


def _prepare(source, folder, ligand_chain, params, template):
    from .validation import validate_ligand_against_params
    folder.mkdir(parents=True, exist_ok=True)
    prepared = folder / "prepared.pdb"
    _PREPARER.rename_ligand_atoms(source, prepared, ligand_chain, "X", template, params)
    validation = validate_ligand_against_params(prepared, params, "X", template=template)
    _save(folder / "ligand_validation.json", validation)
    if not validation["valid"]:
        raise ValueError(f"Ligand validation failed; see {folder / 'ligand_validation.json'}")
    return prepared


def _run_one(source, output, pdb_output, ligand_chain, params, template, protocol, parent_record=None):
    source, output, pdb_output = Path(source), Path(output), Path(pdb_output)
    started = perf_counter()
    identity = {"protocol": protocol, "source_sha256": digest(source), "input_ligand_chain": ligand_chain}
    cached = reuse_case(output / "result.json", identity)
    if cached is not None:
        return cached
    prepared = _prepare(source, output, ligand_chain, Path(params), Path(template))
    identity["pdb_sha256"] = digest(prepared)
    # Retain the stable raw-source key for resume, and validate the prepared
    # file through artifact hashes rather than requiring it before reuse.
    if parent_record:
        previous = read_json(parent_record)
        validate_reusable_parent(previous, protocol, prepared)
        pdb_output.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(native_path(previous["optimized_pdb"]), pdb_output)
        result = {**previous, "identity": identity, "input_pdb": str(prepared),
                  "optimized_pdb": str(pdb_output), "source_pdb": str(source),
                  "reused_parent_result": str(parent_record), "elapsed_seconds": perf_counter()-started,
                  "artifact_sha256": {str(pdb_output): digest(pdb_output), str(prepared): digest(prepared)}}
        _save(output / "result.json", result)
        return result
    pose = _ENGINE.pose_from_file(str(prepared))
    imported = _validate_import(prepared, pose)
    original = _snapshot(pose)
    before = _score_metrics(pose, "X", output / "before_interface.sc")
    baseline = pose.clone()
    _ENGINE.rosetta.numeric.random.rg().set_seed(protocol["seed"])
    task, supported = _hydrogen_task(_ENGINE, pose)
    _SCORE(pose)
    opt_started = perf_counter()
    _ENGINE.rosetta.core.pack.pack_rotamers(pose, _SCORE, task)
    optimization_seconds = perf_counter()-opt_started
    trial = _snapshot(pose)
    if not all(math.isfinite(value) for atom in trial["atoms"].values() for value in atom["xyz"]):
        raise ValueError("Non-finite optimized coordinates")
    trial_audit = compare_snapshots(original, trial)
    trial_bound = float(_SCORE(pose))
    if not math.isfinite(trial_bound):
        raise ValueError("Non-finite optimized bound energy")
    selected = int(trial_bound < before["bound_total"])
    if not selected:
        pose = baseline
    pdb_output.parent.mkdir(parents=True, exist_ok=True)
    temporary = pdb_output.with_suffix(".tmp.pdb")
    pose.dump_pdb(str(temporary))
    temporary.replace(pdb_output)
    persisted = _ENGINE.pose_from_file(str(pdb_output))
    audit = compare_snapshots(original, _snapshot(persisted))
    after = _score_metrics(persisted, "X", output / "interface.sc")
    result = {"status": "complete", "identity": identity, "source_pdb": str(source),
              "input_pdb": str(prepared), "optimized_pdb": str(pdb_output),
              "before_metrics": before, "metrics": after, "import_audit": imported,
              "optimization": {"selected_trial": selected, "optimization_seconds": optimization_seconds,
                               "trial_bound_total_REU": trial_bound, "trial_audit": trial_audit,
                               "persisted_audit": audit, "supported_residues": supported,
                               "max_heavy_displacement_A": audit["max_heavy_displacement_A"],
                               "moved_hydrogen_count": audit["moved_hydrogen_count"]},
              "elapsed_seconds": perf_counter()-started,
              "artifact_sha256": {str(path): digest(path) for path in
                                   (prepared, pdb_output, output/"before_interface.sc", output/"interface.sc")}}
    _save(output / "result.json", result)
    return result


def _protocol(params: Path, template: Path, seed: int, weights="ligand", legacy_pre_talaris=None):
    from .pyrosetta_score import scoring_weights_identity
    database = Path(util.find_spec("pyrosetta").origin).parent / "database"
    scoring = scoring_weights_identity(weights, legacy_pre_talaris, database=database)
    return {"version": 1, "seed": seed, "passes_per_complex": 1,
            "selection": "lowest bound energy of original and one optimized trial",
            "optimization": "whole-complex supported proton chis; fixed heavy atoms and chemical states",
            "score_persisted_pdb": True,
            "source_protocol": {"weights": scoring["weights"], "legacy_pre_talaris": scoring["legacy_pre_talaris"],
                                "params_sha256": digest(params), "template_sha256": digest(template),
                                "weights_sha256": scoring["weights_sha256"],
                                "pyrosetta_version": metadata.version("pyrosetta"), "ligand_chain": "X"},
            "implementation_sha256": {name: digest(Path(__file__).with_name(name)) for name in
                                      ("hydrogen_batch.py", "hydrogen_benchmark.py", "hydrogen_rescore.py",
                                       "pyrosetta_score.py", "hbonds.py", "ligand_chemistry.py", "validation.py")},
            "preparation_implementation_sha256": digest(resource_path("scripts/prepare_ligand_complexes.py"))}


def run(args):
    started = perf_counter()
    parent, folder, output, params, template = [Path(value).resolve() for value in
                                              (args.parent, args.input_folder, args.output, args.params, args.template)]
    if args.jobs < 1 or output == folder or output in folder.parents or folder in output.parents:
        raise ValueError("Use positive jobs and a separate output folder")
    for path in (parent, params, template):
        if not path.is_file():
            raise FileNotFoundError(path)
    candidates = list_inputs(folder, parent, args.recursive, args.glob)
    if args.expected_count is not None and len(candidates) != args.expected_count:
        raise ValueError(f"Expected {args.expected_count} candidates, found {len(candidates)}")
    weights = getattr(args, "weights", "ligand")
    legacy_pre_talaris = getattr(args, "legacy_pre_talaris", None)
    protocol = _protocol(params, template, args.seed, weights, legacy_pre_talaris)
    output.mkdir(parents=True, exist_ok=True)
    manifest_path = output/"hydrogen_batch_manifest.json"
    if manifest_path.exists() and read_json(manifest_path)["protocol"] != protocol:
        raise ValueError("Existing output uses a different OptH protocol")
    identity = {"parent_sha256": digest(parent), "inputs": {path.name: digest(path) for path in candidates}}
    if manifest_path.exists():
        old_inputs = read_json(manifest_path).get("input_identity", {})
        if (old_inputs.get("parent_sha256") != identity["parent_sha256"]
                or any(identity["inputs"].get(name) != value for name, value in old_inputs.get("inputs", {}).items())):
            raise ValueError("Existing output has changed parent/candidate inputs or omitted previously selected candidates")
    manifest = {"status": "running", "started_utc": datetime.now(timezone.utc).isoformat(),
                "parent": str(parent), "input_folder": str(folder), "candidate_count": len(candidates),
                "input_identity": identity, "protocol": protocol, "jobs": args.jobs,
                "completed_count": 0, "failure_count": 0, "cases": []}
    _save(manifest_path, manifest)
    structures = output/"hydrogen_optimized_structures"
    failures, results = [], []
    try:
        with ProcessPoolExecutor(max_workers=args.jobs, mp_context=multiprocessing.get_context("spawn"),
                                 initializer=_init_worker, initargs=(str(params), weights, legacy_pre_talaris)) as pool:
            parent_result = pool.submit(_run_one, str(parent), str(output/"parent"), str(structures/"parent_m7.pdb"),
                                        args.parent_ligand_chain, str(params), str(template), protocol,
                                        str(Path(args.reuse_parent_result).resolve()) if args.reuse_parent_result else None).result()
            manifest["parent_optimized_pdb"] = parent_result["optimized_pdb"]
            _save(manifest_path, manifest)
            futures = {pool.submit(_run_one, str(path), str(output/"cases"/path.stem), str(structures/path.name),
                                   args.ligand_chain, str(params), str(template), protocol): path for path in candidates}
            with (output/"progress.jsonl").open("a", encoding="utf-8") as log:
                for future in as_completed(futures):
                    path = futures[future]
                    try:
                        result = future.result()
                        event = {"candidate": path.name, "status": "complete", "optimized_pdb": result["optimized_pdb"]}
                        results.append(event)
                    except Exception as error:
                        event = {"candidate": path.name, "status": "failed", "error_type": type(error).__name__,
                                 "error": str(error), "traceback": traceback.format_exc(),
                                 "diagnostic_folder": str(output/"cases"/path.stem)}
                        failures.append(event)
                        (output/"cases"/path.stem).mkdir(parents=True, exist_ok=True)
                        _save(output/"cases"/path.stem/"failure.json", event)
                    log.write(json.dumps(event)+"\n")
                    log.flush()
                    count = len(results)+len(failures)
                    if count % 100 == 0 or count == len(candidates):
                        manifest.update(completed_count=len(results), failure_count=len(failures),
                                        elapsed_seconds=perf_counter()-started)
                        _save(manifest_path, manifest)
        manifest.update(status="complete_with_failures" if failures else "complete", completed_count=len(results),
                        failure_count=len(failures), cases=sorted(results+failures, key=lambda r:r["candidate"]),
                        structures_folder=str(structures), elapsed_seconds=perf_counter()-started,
                        finished_utc=datetime.now(timezone.utc).isoformat())
        _save(output/"failures.json", failures)
        _save(manifest_path, manifest)
        print(json.dumps({key: manifest[key] for key in ("status","completed_count","failure_count","structures_folder")}))
        return 1 if failures else 0
    except BaseException as error:
        manifest.update(status="failed", error=str(error), traceback=traceback.format_exc())
        _save(manifest_path, manifest)
        raise


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for option in ("parent", "input-folder", "output", "params", "template"):
        parser.add_argument("--"+option, type=Path, required=True)
    parser.add_argument("--parent-ligand-chain", default="d")
    parser.add_argument("--ligand-chain", default="d")
    parser.add_argument("--reuse-parent-result", type=Path)
    parser.add_argument("--jobs", type=int, default=3)
    parser.add_argument("--seed", type=int, default=20261022)
    parser.add_argument("--weights", default="ligand", help="Rosetta built-in weight set or custom .wts file")
    parser.add_argument("--legacy-pre-talaris", action=argparse.BooleanOptionalAction, default=None,
                        help="Override historical scoring defaults; automatic for ligand/ligand_soft_rep")
    parser.add_argument("--expected-count", type=int)
    parser.add_argument("--recursive", action="store_true")
    parser.add_argument("--glob", default="*.pdb")
    return run(parser.parse_args())


if __name__ == "__main__":
    raise SystemExit(main())
