"""Benchmark whole-complex proton-chi optimization on a saved parent/mutant pair.

This is a separate diagnostic, not a modification of production batch results.
All residues are eligible for supported hydrogen rotamers; heavy atoms, chemical
states and histidine tautomers must remain unchanged. No separated-state packing
or heavy-atom minimization is performed.
"""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
from hashlib import sha256
from importlib.metadata import version
import json
import math
from pathlib import Path
from statistics import median
from time import perf_counter

from .paths import native_path


def compare_snapshots(before: dict, after: dict, tolerance: float = 1e-6) -> dict:
    """Reject any chemical-state or heavy-atom changes, and describe H motion."""
    if before["residue_types"] != after["residue_types"]:
        raise ValueError("Residue types/protonation/tautomers changed")
    if before["atoms"].keys() != after["atoms"].keys():
        raise ValueError("Atom identities/counts changed")
    heavy, hydrogen, moved, moved_residues = [], [], 0, set()
    for key, original in before["atoms"].items():
        updated = after["atoms"][key]
        if original["is_hydrogen"] != updated["is_hydrogen"]:
            raise ValueError("Atom element classification changed")
        displacement = math.dist(original["xyz"], updated["xyz"])
        if original["is_hydrogen"]:
            hydrogen.append(displacement)
            if displacement > tolerance:
                moved += 1
                moved_residues.add(original["residue"])
        else:
            heavy.append(displacement)
    if max(heavy, default=0.0) > tolerance:
        raise ValueError(f"Heavy atoms moved by {max(heavy):.8g} A")
    return {
        "heavy_atom_count": len(heavy), "hydrogen_atom_count": len(hydrogen),
        "max_heavy_displacement_A": max(heavy, default=0.0),
        "max_hydrogen_displacement_A": max(hydrogen, default=0.0),
        "moved_hydrogen_count": moved,
        "moved_hydrogen_residues": sorted(moved_residues),
    }


def _snapshot(pose) -> dict:
    info, atoms, types = pose.pdb_info(), {}, []
    for index in range(1, pose.size() + 1):
        residue = pose.residue(index)
        types.append(str(residue.type().name()))
        label = f"{info.chain(index)}:{info.number(index)}:{residue.name3()}"
        for atom in range(1, residue.natoms() + 1):
            xyz = residue.xyz(atom)
            atoms[index, residue.atom_name(atom).strip()] = {
                "xyz": (float(xyz.x), float(xyz.y), float(xyz.z)),
                "is_hydrogen": bool(residue.atom_is_hydrogen(atom)),
                "residue": label,
            }
    return {"residue_types": types, "atoms": atoms}


def _import_audit(path: Path, pose) -> dict:
    """Separate any import effects from the measured optimization effects."""
    original = {}
    for line in path.read_text().splitlines():
        if line.startswith(("ATOM  ", "HETATM")):
            key = (line[21], int(line[22:26]), line[12:16].strip())
            original[key] = (float(line[30:38]), float(line[38:46]), float(line[46:54]))
    info, displacements, missing = pose.pdb_info(), [], []
    for (index, name), atom in _snapshot(pose)["atoms"].items():
        key = (str(info.chain(index)), int(info.number(index)), name)
        if key in original:
            displacements.append(math.dist(original.pop(key), atom["xyz"]))
        else:
            missing.append(":".join(map(str, key)))
    return {
        "shared_atom_count": len(displacements),
        "max_shared_atom_displacement_A": max(displacements, default=0.0),
        "pose_atoms_absent_from_input": missing,
        "input_atoms_absent_from_pose": [":".join(map(str, key)) for key in original],
    }


def _hydrogen_task(pyrosetta, pose):
    task = pyrosetta.rosetta.core.pack.task.TaskFactory.create_packer_task(pose)
    task.or_optimize_h_mode(True)
    task.or_include_current(True)
    supported = []
    for index in range(1, pose.size() + 1):
        residue_task = task.nonconst_residue_task(index)
        residue_task.or_fix_his_tautomer(True)
        if residue_task.flip_HNQ():
            raise ValueError("Heavy-atom H/N/Q flips must be disabled")
        if pose.residue_type(index).n_proton_chi():
            supported.append({
                "pose_index": index,
                "residue": f"{pose.pdb_info().chain(index)}:{pose.pdb_info().number(index)}:{pose.residue(index).name3()}",
                "proton_chi_count": int(pose.residue_type(index).n_proton_chi()),
            })
    return task, supported


def _diagnostics(pyrosetta, pose, scorefxn, chain: str) -> dict:
    """Report the Ser35/O17 geometry and weighted residue-pair energies."""
    from .pyrosetta_score import _active_terms, _term_name

    info = pose.pdb_info()
    ligand = [i for i in range(1, pose.size() + 1) if info.chain(i) == chain]
    serine = [i for i in range(1, pose.size() + 1)
              if info.chain(i) == "H" and info.number(i) == 35]
    if len(ligand) != 1 or len(serine) != 1:
        raise ValueError("This diagnostic expects one ligand and canonical H:35")
    ser, lig = pose.residue(serine[0]), pose.residue(ligand[0])
    if ser.name3() != "SER":
        raise ValueError("H:35 must be serine for this benchmark diagnostic")
    energies = pyrosetta.rosetta.core.scoring.EMapVector()
    scorefxn.eval_ci_2b(ser, lig, pose, energies)
    scorefxn.eval_cd_2b(ser, lig, pose, energies)
    weighted = {_term_name(term): float(energies[term]) * float(scorefxn.weights()[term])
                for term in _active_terms(scorefxn)}
    def distance(first, second):
        return float(ser.xyz(first).distance(lig.xyz(second)))
    return {
        "ser35_ligand_weighted_pair_terms_REU": weighted,
        "ser35_ligand_pair_total_REU": sum(weighted.values()),
        "ser35_HG_ligand_H2_distance_A": distance("HG", "H2"),
        "ser35_OG_ligand_O2_distance_A": distance("OG", "O2"),
        "ser35_HG_ligand_O2_distance_A": distance("HG", "O2"),
        "ser35_OG_ligand_H2_distance_A": distance("OG", "H2"),
        "ligand_proton_chis_degrees": {
            str(chi): float(lig.chi(chi)) for chi in range(1, lig.nchi() + 1)
            if lig.type().is_proton_chi(chi)
        },
    }


def _measure(pyrosetta, pose, scorefxn, chain: str, path: Path) -> dict:
    from .pyrosetta_score import write_interface_metrics

    started = perf_counter()
    write_interface_metrics(pose, scorefxn, chain, path, path.stem)
    lines = path.read_text().splitlines()
    values = {key: float(value) for key, value in
              zip(lines[0].split()[1:-1], lines[1].split()[1:-1])}
    return {"metrics": values, "diagnostics": _diagnostics(pyrosetta, pose, scorefxn, chain),
            "scoring_seconds": perf_counter() - started}


def _summary_markdown(report: dict) -> str:
    prefix = f"if_{report['protocol']['ligand_chain']}_"
    text = ["# Whole-complex hydrogen-only benchmark", "",
            "All supported proton-chi orientations were sampled in the intact complex. "
            "Heavy atoms and chemical states were fixed; no separated-state optimization was performed.", "",
            f"Score function: `{report['protocol']['weights']}`; "
            f"initialization: {report['initialization_seconds']:.3f} s; "
            f"independent repeats per complex: {report['protocol']['repeats']}.", "",
            "| Complex | Median optimization (s) | Interface before → selected (REU) | "
            "Interface repulsion before → selected (REU) | Ser35 H–H before → selected (Å) |",
            "| --- | ---: | ---: | ---: | ---: |"]
    for name, case in report["cases"].items():
        before, after = case["before"], case["selected"]
        text.append(
            f"| {name} | {median(t['optimization_seconds'] for t in case['trials']):.3f} | "
            f"{before['metrics']['interface_delta_X']:.4f} → {after['metrics']['interface_delta_X']:.4f} | "
            f"{before['metrics'][prefix+'fa_rep']:.4f} → {after['metrics'][prefix+'fa_rep']:.4f} | "
            f"{before['diagnostics']['ser35_HG_ligand_H2_distance_A']:.3f} → "
            f"{after['diagnostics']['ser35_HG_ligand_H2_distance_A']:.3f} |")
    parent, mutant = report["cases"]["parent"], report["cases"]["mutant"]
    before_delta = mutant["before"]["metrics"]["interface_delta_X"] - parent["before"]["metrics"]["interface_delta_X"]
    after_delta = mutant["selected"]["metrics"]["interface_delta_X"] - parent["selected"]["metrics"]["interface_delta_X"]
    text += ["", f"Mutant − parent interface difference: **{before_delta:+.4f} → {after_delta:+.4f} REU**.", "",
             "Selection uses the lowest whole-complex bound score, not the most favorable interface score. "
             "The original is retained if no trial improves its bound score. Every accepted trial must "
             "preserve all heavy coordinates, residue types, atom identities and hydrogen counts.", "",
             "These are diagnostic results only. Production workbooks and detector fingerprints were not updated.", ""]
    return "\n".join(text)


def run_benchmark(rescore_root: Path, candidate: str, params: Path, output: Path,
                  repeats: int = 3, seed: int = 20261022) -> dict:
    if repeats < 1:
        raise ValueError("At least one repeat is required")
    if output.exists():
        raise FileExistsError(f"Refusing to overwrite benchmark output: {output}")
    records = {
        "parent": json.loads((rescore_root / "parent/result.json").read_text()),
        "mutant": json.loads((rescore_root / Path(candidate).stem / "result.json").read_text()),
    }
    protocol = records["parent"]["identity"]["protocol"]
    if records["mutant"]["identity"]["protocol"] != protocol:
        raise ValueError("Parent and mutant have different scoring protocols")
    if sha256(params.read_bytes()).hexdigest() != protocol["params_sha256"]:
        raise ValueError("Ligand parameter hash differs from saved analysis")
    for record in records.values():
        if sha256(native_path(record["pdb"]).read_bytes()).hexdigest() != record["identity"]["pdb_sha256"]:
            raise ValueError("Saved input structure was modified")

    output.mkdir(parents=True)
    started = perf_counter()
    from .pyrosetta_score import _initialize, _legacy_weights, _scorefunction
    pyrosetta = _initialize([str(params)], legacy=_legacy_weights(protocol["weights"]))
    pyrosetta.rosetta.basic.options.set_boolean_option("packing:no_optH", True)
    pyrosetta.rosetta.basic.options.set_boolean_option("packing:flip_HNQ", False)
    scorefxn = _scorefunction(pyrosetta, protocol["weights"])
    report = {
        "status": "running", "started_utc": datetime.now(timezone.utc).isoformat(),
        "initialization_seconds": perf_counter() - started,
        "protocol": {"weights": protocol["weights"], "legacy_pre_talaris": protocol["legacy_pre_talaris"],
                     "pyrosetta_version": version("pyrosetta"), "params_sha256": protocol["params_sha256"],
                     "repeats": repeats, "first_seed": seed, "ligand_chain": protocol["ligand_chain"],
                     "scope": "whole-complex supported proton chis; fixed heavy atoms and chemical states",
                     "optimization": "optimize-H pack_rotamers; no minimization or HNQ flips",
                     "implementation_sha256": sha256(Path(__file__).read_bytes()).hexdigest()},
        "cases": {},
    }
    if report["protocol"]["pyrosetta_version"] != protocol["pyrosetta_version"]:
        raise ValueError("PyRosetta version differs from saved scoring")
    for name, record in records.items():
        case_path = output / name
        case_path.mkdir()
        load_start = perf_counter()
        original = pyrosetta.pose_from_file(str(native_path(record["pdb"])))
        load_seconds = perf_counter() - load_start
        baseline_snapshot = _snapshot(original)
        before = _measure(pyrosetta, original, scorefxn, protocol["ligand_chain"], case_path / "before.sc")
        for key, value in record["metrics"].items():
            if not math.isclose(before["metrics"][key], value, abs_tol=1e-4):
                raise ValueError(f"Baseline scoring does not reproduce saved {name} {key}")
        case = {"input_pdb": record["pdb"], "input_sha256": record["identity"]["pdb_sha256"],
                "load_seconds": load_seconds, "import_audit": _import_audit(native_path(record["pdb"]), original),
                "before": before, "trials": []}
        report["cases"][name] = case
        best, best_measure, best_trial = original, before, "original"
        for repeat in range(repeats):
            pose = original.clone()
            pyrosetta.rosetta.numeric.random.rg().set_seed(seed + repeat)
            trial_start = perf_counter()
            task, supported = _hydrogen_task(pyrosetta, pose)
            scorefxn(pose)
            pyrosetta.rosetta.core.pack.pack_rotamers(pose, scorefxn, task)
            optimization_seconds = perf_counter() - trial_start
            audit = compare_snapshots(baseline_snapshot, _snapshot(pose))
            after = _measure(pyrosetta, pose, scorefxn, protocol["ligand_chain"],
                             case_path / f"trial_{repeat + 1}.sc")
            trial = {"repeat": repeat + 1, "seed": seed + repeat,
                     "optimization_seconds": optimization_seconds, "supported_residues": supported,
                     "audit": audit, **after}
            case["trials"].append(trial)
            write_start = perf_counter()
            pose.dump_pdb(str(case_path / f"trial_{repeat + 1}.pdb"))
            trial["pdb_write_seconds"] = perf_counter() - write_start
            if after["metrics"]["bound_total"] < best_measure["metrics"]["bound_total"]:
                best, best_measure, best_trial = pose, after, repeat + 1
            print(json.dumps({"case": name, "repeat": repeat + 1,
                              "optimization_seconds": optimization_seconds,
                              "interface_REU": after["metrics"]["interface_delta_X"],
                              "fa_rep_REU": after["metrics"]["if_X_fa_rep"],
                              "max_heavy_displacement_A": audit["max_heavy_displacement_A"]}), flush=True)
        best.dump_pdb(str(case_path / "selected.pdb"))
        case["selected_trial"] = best_trial
        case["selected"] = best_measure
        (output / "benchmark.json").write_text(json.dumps(report, indent=2) + "\n")
    report["status"] = "complete"
    report["elapsed_seconds"] = perf_counter() - started
    report["finished_utc"] = datetime.now(timezone.utc).isoformat()
    (output / "benchmark.json").write_text(json.dumps(report, indent=2) + "\n")
    (output / "benchmark_summary.md").write_text(_summary_markdown(report), encoding="utf-8")
    print(json.dumps({"status": "complete", "output": str(output),
                      "elapsed_seconds": report["elapsed_seconds"]}), flush=True)
    return report


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--rescore-root", type=Path, required=True)
    parser.add_argument("--candidate", required=True)
    parser.add_argument("--params", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--repeats", type=int, default=3)
    parser.add_argument("--seed", type=int, default=20261022)
    args = parser.parse_args()
    run_benchmark(args.rescore_root, args.candidate, args.params, args.output, args.repeats, args.seed)


if __name__ == "__main__":
    main()
