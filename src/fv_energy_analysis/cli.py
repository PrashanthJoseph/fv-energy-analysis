from __future__ import annotations

import argparse
import math
from pathlib import Path

from .arpeggio import read_arpeggio
from .analysis_stages import add_analysis_arguments, analysis_status
from .chemistry import substitution_chemistry
from .geometry import BACKBONE, structural_rmsds
from .hbonds import aggregate_hbond_csv, extract_pyrosetta_hbonds
from .network import mutation_network_changes
from .numbering import annotate_fv
from .pyrosetta_score import score_and_extract_hbonds, score_with_pyrosetta
from .report import merge_by_pair, outer_compare, write_csv, write_json
from .rosetta import read_residue_energy_breakdown, read_scorefile_metrics, run_residue_energy_breakdown
from .rotamers import compare_rotamers
from .structure import THREE_TO_ONE, read_pdb, write_transformed_pdb


def _parse_chains(text: str) -> list[str]:
    values = [item.strip() for item in text.split(",") if item.strip()]
    if not values:
        raise argparse.ArgumentTypeError("At least one chain ID is required")
    return values


def _validated_pairs(parent, mutant, chains: list[str]):
    pairs = {}
    for chain in chains:
        p = parent.chain_residues(chain)
        m = mutant.chain_residues(chain)
        if not p or not m:
            raise ValueError(f"Chain {chain!r} must occur as a protein chain in both structures")
        if len(p) != len(m):
            raise ValueError(f"Chain {chain}: parent has {len(p)} residues but mutant has {len(m)}")
        p_positions = [(res.key.number, res.key.insertion) for res in p]
        m_positions = [(res.key.number, res.key.insertion) for res in m]
        if p_positions != m_positions:
            raise ValueError(
                f"Chain {chain}: PDB residue numbers/insertion codes differ. Structures must have "
                "one-to-one residue identity by position, not merely equal lengths."
            )
        pairs[chain] = (p, m)
    return pairs


ROTAMER_FIELDS = (
    "parent_rotamer_bin", "mutant_rotamer_bin", "parent_chi_degrees",
    "mutant_chi_degrees", "n_comparable_chi", "same_chi1_bin",
    "same_common_chi_bins", "same_full_rotamer_bin",
    "rotamer_method", "parent_chi_complete", "mutant_chi_complete", "rotamer_comparison_note",
)


def _mutation_rows(
    chain_pairs,
    p_annotations,
    m_annotations,
    rmsd_by_region: dict[str, float],
    rotamer_rmsd_cutoff: float,
    perform_rotamers: bool = True,
):
    rows = []
    for chain, (p_residues, m_residues) in chain_pairs.items():
        for sequence_index, (p_res, m_res) in enumerate(zip(p_residues, m_residues), start=1):
            p_aa, m_aa = THREE_TO_ONE[p_res.name], THREE_TO_ONE[m_res.name]
            if p_aa == m_aa:
                continue
            p_region = p_annotations[chain].region_by_key.get(p_res.key, "UNNUMBERED")
            m_region = m_annotations[chain].region_by_key.get(m_res.key, "UNNUMBERED")
            p_number = p_annotations[chain].number_by_key.get(p_res.key, "")
            m_number = m_annotations[chain].number_by_key.get(m_res.key, "")
            row = {
                "chain": chain,
                "sequence_index": sequence_index,
                "parent_residue": p_res.key.label(),
                "mutant_residue": m_res.key.label(),
                "parent_aa": p_aa,
                "mutant_aa": m_aa,
                "substitution": f"{p_aa}{sequence_index}{m_aa}",
                "parent_antibody_number": p_number,
                "mutant_antibody_number": m_number,
                "parent_region": p_region,
                "mutant_region": m_region,
                "same_numbered_region": p_region == m_region,
            }
            row.update(substitution_chemistry(p_aa, m_aa))
            region_label = f"{chain}_{p_region}"
            cdr_rmsd = rmsd_by_region.get(region_label) if p_region.startswith("CDR") else None
            row["cdr_global_fit_rmsd_A"] = cdr_rmsd if cdr_rmsd is not None else ""
            gate_region = (f"{chain}_CDR3_ANCHORS"
                           if p_region == "CDR3" and p_annotations[chain].chain_type == "H"
                           else region_label if p_region.startswith("CDR") else "")
            gate_rmsd = rmsd_by_region.get(gate_region) if gate_region else None
            row["rotamer_gate_region"] = gate_region
            row["rotamer_gate_rmsd_A"] = gate_rmsd if gate_rmsd is not None else ""
            row["rotamer_rmsd_cutoff_A"] = rotamer_rmsd_cutoff
            if not perform_rotamers:
                row["rotamer_analysis_status"] = "not_run"
                row.update({field: "" for field in ROTAMER_FIELDS})
            elif p_region.startswith("CDR") and (gate_rmsd is None or not math.isfinite(gate_rmsd)):
                row["rotamer_analysis_status"] = "skipped_cdr_rmsd_unavailable"
                row.update({field: "" for field in ROTAMER_FIELDS})
            elif gate_rmsd is not None and gate_rmsd > rotamer_rmsd_cutoff:
                row["rotamer_analysis_status"] = "skipped_cdr_rmsd_above_cutoff"
                row.update({field: "" for field in ROTAMER_FIELDS})
            else:
                row["rotamer_analysis_status"] = "performed"
                row.update(compare_rotamers(p_res, m_res))
            rows.append(row)
    return rows


def _annotate_pair_context(rows, fv_chains: set[str], ligand_chains: set[str]):
    for row in rows:
        labels = [str(row.get("residue1", "")), str(row.get("residue2", ""))]
        chains = [label.split(":", 1)[0] if ":" in label else "" for label in labels]
        row["pair_scope"] = (
            "aa_to_ligand" if any(c in ligand_chains for c in chains)
            else "aa_to_aa" if all(c in fv_chains for c in chains)
            else "fv_to_other"
        )
    return rows


def compare(args: argparse.Namespace) -> int:
    output = Path(args.output)
    output.mkdir(parents=True, exist_ok=True)
    parent = read_pdb(args.parent)
    mutant = read_pdb(args.mutant)
    chain_pairs = _validated_pairs(parent, mutant, args.fv_chains)

    p_annotations = annotate_fv(parent, args.fv_chains, args.scheme, args.cdr_definition)
    m_annotations = annotate_fv(mutant, args.fv_chains, args.scheme, args.cdr_definition)
    for chain, (p_residues, m_residues) in chain_pairs.items():
        for p_res, m_res in zip(p_residues, m_residues):
            p_number = p_annotations[chain].number_by_key.get(p_res.key)
            m_number = m_annotations[chain].number_by_key.get(m_res.key)
            p_region = p_annotations[chain].region_by_key.get(p_res.key)
            m_region = m_annotations[chain].region_by_key.get(m_res.key)
            if (p_number, p_region) != (m_number, m_region):
                raise ValueError(
                    f"ANARCII mapping differs at chain {chain} residue {p_res.key.number}{p_res.key.insertion}: "
                    f"parent={(p_number, p_region)}, mutant={(m_number, m_region)}"
                )
    p_regions = {key: region for ann in p_annotations.values() for key, region in ann.region_by_key.items()}

    transform, rmsd_rows = structural_rmsds(
        parent, mutant, chain_pairs, p_regions, fit_region=args.fit_region,
        atom_names=("CA",) if args.rmsd_atoms == "ca" else BACKBONE,
        heavy_chain_ids=tuple(chain for chain, ann in p_annotations.items() if ann.chain_type == "H"),
        hcdr3_tip_sequence=getattr(args, "hcdr3_tip_sequence", "YYGSS"),
    )
    write_csv(output / "structure_rmsd.csv", rmsd_rows)
    write_transformed_pdb(mutant, transform.rotation, transform.translation, output / "mutant_aligned_to_parent.pdb")

    rmsd_by_region = {
        str(row["region"]): float(row["rmsd_after_global_fit_A"])
        for row in rmsd_rows
        if row["rmsd_after_global_fit_A"] != ""
    }
    mutation_rows = _mutation_rows(
        chain_pairs, p_annotations, m_annotations, rmsd_by_region,
        args.rotamer_rmsd_cutoff,
        perform_rotamers=not getattr(args, "skip_rotamer_analysis", False),
    )
    write_csv(output / "substitutions.csv", mutation_rows)

    fv_set, ligand_set = set(args.fv_chains), set(args.ligand_chains)
    p_energy = read_residue_energy_breakdown(args.parent_rosetta, parent) if args.parent_rosetta else []
    m_energy = read_residue_energy_breakdown(args.mutant_rosetta, mutant) if args.mutant_rosetta else []
    p_contacts = read_arpeggio(args.parent_arpeggio, fv_set) if args.parent_arpeggio and not getattr(args, "skip_arpeggio", False) else []
    m_contacts = read_arpeggio(args.mutant_arpeggio, fv_set) if args.mutant_arpeggio and not getattr(args, "skip_arpeggio", False) else []
    p_hbonds = aggregate_hbond_csv(args.parent_hbonds) if args.parent_hbonds else []
    m_hbonds = aggregate_hbond_csv(args.mutant_hbonds) if args.mutant_hbonds else []
    _annotate_pair_context(p_energy, fv_set, ligand_set)
    _annotate_pair_context(m_energy, fv_set, ligand_set)
    _annotate_pair_context(p_contacts, fv_set, ligand_set)
    _annotate_pair_context(m_contacts, fv_set, ligand_set)
    _annotate_pair_context(p_hbonds, fv_set, ligand_set)
    _annotate_pair_context(m_hbonds, fv_set, ligand_set)

    if p_energy:
        write_csv(output / "rosetta_pairs_parent.csv", p_energy)
    if m_energy:
        write_csv(output / "rosetta_pairs_mutant.csv", m_energy)
    if p_contacts:
        write_csv(output / "interactions_parent.csv", p_contacts)
    if m_contacts:
        write_csv(output / "interactions_mutant.csv", m_contacts)

    p_merged = merge_by_pair(p_energy, p_contacts, p_hbonds)
    m_merged = merge_by_pair(m_energy, m_contacts, m_hbonds)
    comparison = outer_compare(p_merged, m_merged) if (p_merged or m_merged) else []
    if comparison:
        write_csv(output / "pair_comparison.csv", comparison)
        mutated_labels = {str(row["parent_residue"]) for row in mutation_rows} | {str(row["mutant_residue"]) for row in mutation_rows}
        mutation_impacts = [
            row for row in comparison
            if any(str(row.get(col, "")) in mutated_labels for col in (
                "parent_residue1", "parent_residue2", "mutant_residue1", "mutant_residue2"
            ))
        ]
        write_csv(output / "mutation_partner_impacts.csv", mutation_impacts)
        if not getattr(args, "skip_network_analysis", False):
            network_details, network_summary, network_type_energy = mutation_network_changes(
                comparison,
                mutation_rows,
                energy_presence_cutoff=args.network_energy_presence_cutoff,
                energy_change_cutoff=args.network_energy_change_cutoff,
            )
            write_csv(output / "mutation_network_changes.csv", network_details)
            write_csv(output / "mutation_network_summary.csv", network_summary)
            write_csv(output / "mutation_network_type_energy_summary.csv", network_type_energy)

    if args.parent_interface or args.mutant_interface:
        if not (args.parent_interface and args.mutant_interface):
            raise ValueError("Provide both --parent-interface and --mutant-interface")
        p_metrics = read_scorefile_metrics(args.parent_interface)
        m_metrics = read_scorefile_metrics(args.mutant_interface)
        metric_rows = []
        for name in sorted(set(p_metrics) | set(m_metrics)):
            pv, mv = p_metrics.get(name, float("nan")), m_metrics.get(name, float("nan"))
            metric_rows.append({
                "metric": name,
                "parent": pv,
                "mutant": mv,
                "delta_mutant_minus_parent": mv - pv if not (math.isnan(pv) or math.isnan(mv)) else "",
            })
        write_csv(output / "interface_metrics_comparison.csv", metric_rows)

    stages = analysis_status(args)
    if not (args.parent_arpeggio or args.mutant_arpeggio):
        stages["arpeggio"] = "not_run"
    # The comparison CLI consumes saved pair/contact tables; detector execution
    # and E2 geometry belong to its surrounding scientific workflow.
    stages["plip"] = stages["e2_architecture"] = "not_run"
    summary = {
        "analysis_status": stages,
        "parent": str(Path(args.parent).resolve()),
        "mutant": str(Path(args.mutant).resolve()),
        "fv_chains": args.fv_chains,
        "ligand_chains": args.ligand_chains,
        "numbering_scheme": args.scheme,
        "cdr_definition": args.cdr_definition,
        "numbering_backend": sorted({ann.backend for ann in p_annotations.values()}),
        "fit_region": args.fit_region,
        "rmsd_atoms": args.rmsd_atoms,
        "hcdr3_tip_reference_sequence": getattr(args, "hcdr3_tip_sequence", "YYGSS"),
        "mutation_count": len(mutation_rows),
        "rotamer_rmsd_cutoff_A": args.rotamer_rmsd_cutoff,
        "network_energy_presence_cutoff_REU": args.network_energy_presence_cutoff,
        "network_energy_change_cutoff_REU": args.network_energy_change_cutoff,
        "warnings": [
            "Pair-strength labels are heuristic; use numeric Rosetta REU values for analysis.",
            "Rotamer fields are coarse canonical-chi torsion bins, not Dunbrack assignments; full-bin identity is not assigned across different amino-acid types.",
            "CDR and framework RMSD_after_global_fit values all use the same global transform; local_fit_rmsd is diagnostic only.",
        ],
    }
    write_json(output / "summary.json", summary)
    return 0


def run_rosetta(args: argparse.Namespace) -> int:
    run_residue_energy_breakdown(
        args.executable, args.pdb, args.output, params=args.params,
        weights=args.weights, database=args.database,
    )
    return 0


def extract_hbonds(args: argparse.Namespace) -> int:
    extract_pyrosetta_hbonds(args.pdb, args.output, params=args.params, weights=args.weights,
                            **({"legacy_pre_talaris": args.legacy_pre_talaris} if getattr(args, "legacy_pre_talaris", None) is not None else {}))
    return 0


def score_pyrosetta(args: argparse.Namespace) -> int:
    score_options = {} if getattr(args, "legacy_pre_talaris", None) is None else {"legacy_pre_talaris": args.legacy_pre_talaris}
    if args.hbond_output:
        if not args.interface_output or not args.ligand_chain:
            raise ValueError('--hbond-output requires --interface-output and --ligand-chain')
        score_and_extract_hbonds(
            args.pdb, args.pair_output, args.interface_output, args.hbond_output,
            args.ligand_chain, params=args.params, weights=args.weights, **score_options,
        )
        return 0
    score_with_pyrosetta(
        args.pdb,
        args.pair_output,
        interface_output=args.interface_output,
        ligand_chain=args.ligand_chain,
        params=args.params,
        weights=args.weights, **score_options,
    )
    return 0


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog="fv-energy", description=__doc__)
    sub = parser.add_subparsers(dest="command", required=True)
    comp = sub.add_parser("compare", help="Compare parent and mutant Fv-ligand PDB structures")
    comp.add_argument("--parent", required=True)
    comp.add_argument("--mutant", required=True)
    comp.add_argument("--fv-chains", type=_parse_chains, required=True, help="Comma-separated Fv chains, e.g. H,L")
    comp.add_argument("--ligand-chains", type=_parse_chains, required=True, help="Comma-separated ligand chains, e.g. X")
    comp.add_argument("--scheme", default="imgt", choices=["imgt", "chothia", "kabat", "aho"])
    comp.add_argument("--cdr-definition", default="imgt", choices=["imgt", "chothia", "kabat", "north"])
    comp.add_argument("--fit-region", default="fv", choices=["fv", "framework"])
    comp.add_argument("--rmsd-atoms", default="backbone", choices=["backbone", "ca"])
    comp.add_argument("--rotamer-rmsd-cutoff", type=float, default=3.0)
    comp.add_argument("--hcdr3-tip-sequence", default="YYGSS", help="Unique parent HCDR3 tip motif; mutant positions are inherited from the parent")
    comp.add_argument("--network-energy-presence-cutoff", type=float, default=0.1)
    comp.add_argument("--network-energy-change-cutoff", type=float, default=0.25)
    comp.add_argument("--parent-rosetta", help="Parent residue_energy_breakdown score table")
    comp.add_argument("--mutant-rosetta", help="Mutant residue_energy_breakdown score table")
    comp.add_argument("--parent-arpeggio", help="Parent PDBe Arpeggio JSON")
    comp.add_argument("--mutant-arpeggio", help="Mutant PDBe Arpeggio JSON")
    comp.add_argument("--parent-hbonds", help="Parent CSV from fv-energy extract-hbonds")
    comp.add_argument("--mutant-hbonds", help="Mutant CSV from fv-energy extract-hbonds")
    comp.add_argument("--parent-interface", help="Parent InterfaceScoreCalculator Rosetta scorefile")
    comp.add_argument("--mutant-interface", help="Mutant InterfaceScoreCalculator Rosetta scorefile")
    comp.add_argument("--output", required=True)
    add_analysis_arguments(comp)
    comp.set_defaults(func=compare)

    ros = sub.add_parser("run-rosetta", help="Run Rosetta residue_energy_breakdown")
    ros.add_argument("--executable", required=True)
    ros.add_argument("--pdb", required=True)
    ros.add_argument("--output", required=True)
    ros.add_argument("--params", action="append", default=[])
    ros.add_argument("--weights", default="ref2015")
    ros.add_argument("--database")
    ros.set_defaults(func=run_rosetta)

    hb = sub.add_parser("extract-hbonds", help="Extract individual scored H-bonds using PyRosetta")
    hb.add_argument("--pdb", required=True)
    hb.add_argument("--output", required=True)
    hb.add_argument("--params", action="append", default=[])
    hb.add_argument("--weights", default="ref2015")
    hb.add_argument("--legacy-pre-talaris", action=argparse.BooleanOptionalAction, default=None)
    hb.set_defaults(func=extract_hbonds)

    pyscore = sub.add_parser(
        "score-pyrosetta",
        help="Generate residue-pair and ligand-interface score tables directly with PyRosetta",
    )
    pyscore.add_argument("--pdb", required=True)
    pyscore.add_argument("--pair-output", required=True)
    pyscore.add_argument("--interface-output")
    pyscore.add_argument("--hbond-output", help="Also extract H-bonds from the same loaded pose")
    pyscore.add_argument("--ligand-chain")
    pyscore.add_argument("--params", action="append", default=[])
    pyscore.add_argument("--weights", default="ref2015")
    pyscore.add_argument("--legacy-pre-talaris", action=argparse.BooleanOptionalAction, default=None)
    pyscore.set_defaults(func=score_pyrosetta)
    return parser


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    try:
        return int(args.func(args))
    except (ValueError, RuntimeError, OSError) as exc:
        parser.error(str(exc))
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
