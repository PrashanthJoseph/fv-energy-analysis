"""Rebuild mutation-network reports from saved scoring-only calculations.

The source batch is read-only. Geometry annotations and substitutions are reused,
but pair energies, hydrogen bonds, energy-presence events and quality warnings
come exclusively from the new scorefunction. No scoring engine is initialized.
"""

from __future__ import annotations

import argparse
import csv
from hashlib import sha256
import json
from pathlib import Path

from .batch_workflow import _candidate_network_summary
from .hbonds import aggregate_hbond_csv
from .network import mutation_network_changes
from .paths import native_path
from .report import merge_by_pair, outer_compare, write_csv, write_json
from .rosetta import read_residue_energy_breakdown
from .structure import read_pdb


def _json(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8-sig"))


def _csv(path: Path) -> list[dict]:
    with path.open(newline="", encoding="utf-8-sig") as handle:
        return list(csv.DictReader(handle))


def _geometry(path: Path, mutated: set[str]) -> list[dict]:
    """Use only the saved geometric table, never a previous merged comparison."""
    rows = _csv(path)
    for row in rows:
        labels = set(str(row.get("interaction_types", "")).split(";"))
        if labels & {"energy_only", "polymer_conformational", "rosetta_hbond"}:
            raise ValueError(f"Saved geometry contains score-dependent labels: {path}")
    return _mutation_pairs(rows, mutated)


def _mutation_pairs(rows: list[dict], mutated: set[str]) -> list[dict]:
    return [row for row in rows if mutated.intersection(
        (str(row.get("residue1", "")), str(row.get("residue2", ""))))]


def _context(rows: list[dict], fv: set[str], ligand: set[str]) -> list[dict]:
    for row in rows:
        chains = [str(row.get(f"residue{index}", "")).split(":", 1)[0]
                  for index in (1, 2)]
        row["pair_scope"] = ("aa_to_ligand" if any(c in ligand for c in chains)
                             else "aa_to_aa" if all(c in fv for c in chains)
                             else "fv_to_other")
    return rows


def _validated_result(folder: Path, protocol: dict) -> dict:
    saved = _json(folder / "result.json")
    if saved.get("status") != "complete" or saved["identity"]["protocol"] != protocol:
        raise ValueError(f"Incomplete or inconsistent scoring protocol: {folder}")
    for name, expected in saved["artifact_sha256"].items():
        if sha256((folder / name).read_bytes()).hexdigest() != expected:
            raise ValueError(f"Scoring artifact hash mismatch: {folder / name}")
    pdb = native_path(saved["pdb"])
    if sha256(pdb.read_bytes()).hexdigest() != saved["identity"]["pdb_sha256"]:
        raise ValueError(f"Scored coordinate hash mismatch: {pdb}")
    return saved


def comparison_warning(parent: dict, mutant: dict, ligand_chain: str) -> bool:
    """Same existing heuristic, recalculated within the new scoring function."""
    return (abs(mutant["bound_total"] - parent["bound_total"]) > 100.0
            or mutant.get(f"if_{ligand_chain}_fa_rep", 0.0)
            - parent.get(f"if_{ligand_chain}_fa_rep", 0.0) > 10.0)


def export_network(source: Path, rescore: Path,
                   candidates: list[str] | None = None) -> dict[str, dict]:
    """Save fresh mutation-centred CSV reports and return workbook metrics.

    Candidate keys include ``.pdb``. ``candidates`` optionally limits work for a
    representative validation. Outputs are written only beneath ``rescore``.
    """
    source, rescore = Path(source).resolve(), Path(rescore).resolve()
    if source == rescore or source in rescore.parents:
        raise ValueError("Network export must be outside the source batch")
    rows = _csv(rescore / "rosettaligand_rescore.csv")
    wanted = set(candidates) if candidates is not None else None
    if wanted is not None:
        rows = [row for row in rows if row["candidate"] in wanted]
        missing = wanted - {row["candidate"] for row in rows}
        if missing:
            raise ValueError(f"Candidates absent from rescore: {sorted(missing)}")
    if not rows:
        raise ValueError("No completed candidates to export")
    protocol = _json(rescore / "rescore_manifest.json")["protocol"]
    parent_saved = _validated_result(rescore / "parent", protocol)
    parent_energy = read_residue_energy_breakdown(
        rescore / "parent/pair_energies.sc", read_pdb(native_path(parent_saved["pdb"])))
    parent_hbonds = aggregate_hbond_csv(rescore / "parent/hbonds.csv")
    metrics, aggregate = {}, []
    for row in rows:
        candidate = row["candidate"]
        source_case = native_path(row["source_result_folder"]).resolve()
        if source_case.parent != source:
            raise ValueError(f"Source case lies outside requested source batch: {source_case}")
        case = rescore / Path(candidate).stem
        saved = _validated_result(case, protocol)
        summary = _json(source_case / "final/summary.json")
        substitutions = _csv(source_case / "final/substitutions.csv")
        mutated = {str(mutation[key]) for mutation in substitutions
                   for key in ("parent_residue", "mutant_residue")}
        fv, ligand = set(summary["fv_chains"]), set(summary["ligand_chains"])
        if len(ligand) != 1 or next(iter(ligand)) != protocol["ligand_chain"]:
            raise ValueError(f"Ligand-chain mismatch: {candidate}")
        mutant_energy = read_residue_energy_breakdown(
            case / "pair_energies.sc", read_pdb(native_path(saved["pdb"])))
        p_merged = merge_by_pair(
            _mutation_pairs(parent_energy, mutated),
            _geometry(source_case / "final/interactions_parent.csv", mutated),
            _mutation_pairs(parent_hbonds, mutated))
        m_merged = merge_by_pair(
            _mutation_pairs(mutant_energy, mutated),
            _geometry(source_case / "final/interactions_mutant.csv", mutated),
            _mutation_pairs(aggregate_hbond_csv(case / "hbonds.csv"), mutated))
        comparison = outer_compare(_context(p_merged, fv, ligand), _context(m_merged, fv, ligand))
        details, summaries, type_energies = mutation_network_changes(
            comparison, substitutions,
            energy_presence_cutoff=float(summary.get("network_energy_presence_cutoff_REU", 0.1)),
            energy_change_cutoff=float(summary.get("network_energy_change_cutoff_REU", 0.25)))
        final = case / "final"
        final.mkdir(exist_ok=True)
        write_csv(final / "mutation_partner_impacts.csv", comparison)
        write_csv(final / "mutation_network_changes.csv", details)
        write_csv(final / "mutation_network_summary.csv", summaries)
        write_csv(final / "mutation_network_type_energy_summary.csv", type_energies)
        compact, _ = _candidate_network_summary(candidate, case, details)
        compact["comparison_quality_warning"] = comparison_warning(
            parent_saved["metrics"], saved["metrics"], protocol["ligand_chain"])
        compact["scorefunction"] = protocol["weights"]
        metrics[candidate] = compact
        aggregate.append({"candidate": candidate, **compact})
        write_json(final / "network_metrics.json", compact)
    # A limited check must never replace a previously exported complete table.
    suffix = "_subset" if candidates is not None else ""
    write_csv(rescore / f"rescore_network_summary{suffix}.csv", aggregate)
    write_json(rescore / f"rescore_network_summary{suffix}.json", metrics)
    return metrics


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, required=True)
    parser.add_argument("--rescore", type=Path, required=True)
    parser.add_argument("--candidate", action="append")
    args = parser.parse_args(argv)
    result = export_network(args.source, args.rescore, args.candidate)
    print(json.dumps({"completed": len(result), "output": str(args.rescore)}))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
