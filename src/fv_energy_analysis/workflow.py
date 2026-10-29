"""One-command parent-versus-mutant Fv/ligand analysis workflow.

Run inside the activated WSL environment, for example::

    python -m fv_energy_analysis.workflow --parent parent.pdb --mutant mutant.pdb \
        --output results/mutant_vs_parent

The module owns every stage and writes ``workflow.log`` plus
``workflow_manifest.json`` so callers do not need to monitor intermediate jobs.
"""

from __future__ import annotations

import argparse
import csv
import json
import os
import shlex
import shutil
import subprocess
import sys
import time
from hashlib import sha256
from datetime import datetime, timezone
from pathlib import Path

from .estrogen_architecture import analyze_estrogen_architecture
from .plip import summarize_plip_xml
from .report import write_csv
from .resources import resource_path
from .arpeggio import run_arpeggio, mutation_neighbourhood_selection, filter_arpeggio_contacts
from .analysis_stages import add_analysis_arguments, analysis_flags, analysis_status, record_comparison_status, scoring_protocol, omitted_e2_architecture
from .validation import validate_ligand_against_params


def _chains(text: str) -> str:
    values = [value.strip() for value in text.split(",") if value.strip()]
    if not values:
        raise argparse.ArgumentTypeError("At least one chain is required")
    return ",".join(values)


def _chain_map(text: str) -> dict[str, str]:
    if not text:
        return {}
    mapping: dict[str, str] = {}
    for item in text.split(","):
        try:
            source, target = (value.strip() for value in item.split(":", 1))
        except ValueError as exc:
            raise argparse.ArgumentTypeError(
                "Chain maps must use SOURCE:TARGET, for example H:L,L:H"
            ) from exc
        if not source or not target:
            raise argparse.ArgumentTypeError("Chain map names cannot be empty")
        mapping[source] = target
    return mapping


def _convert_structure_input(
    source: Path,
    destination: Path,
    chain_map: dict[str, str],
    ligand_chain: str,
) -> Path:
    """Convert mmCIF to PDB and/or remap chains for downstream Rosetta tools."""
    if source.suffix.lower() not in {".cif", ".mmcif"}:
        if not chain_map:
            return source
        # A fixed-column rewrite preserves atom serials and CONECT records.
        # Rewriting a PDB through Gemmi discards CONECT, which can make RDKit
        # infer spurious ligand bonds from close coordinates after relaxation.
        lines: list[str] = []
        for line in source.read_text(encoding="utf-8").splitlines(keepends=True):
            if line.startswith(("ATOM  ", "HETATM", "ANISOU", "TER   ")) and len(line) > 21:
                source_chain = line[21:22].strip()
                target_chain = chain_map.get(source_chain, source_chain)
                if len(target_chain) > 1:
                    raise ValueError(f"PDB chain IDs must be one character: {target_chain!r}")
                line = f"{line[:21]}{target_chain or ' '}{line[22:]}"
            lines.append(line)
        destination.write_text("".join(lines), encoding="utf-8")
        return destination
    try:
        import gemmi
    except ImportError as exc:
        raise RuntimeError("gemmi is required for CIF input or chain remapping") from exc
    structure = gemmi.read_structure(str(source))
    for model in structure:
        chains = list(model)
        original_names = [chain.name for chain in chains]
        for chain, original_name in zip(chains, original_names):
            chain.name = chain_map.get(original_name, original_name)
            if chain.name == ligand_chain:
                for residue in chain:
                    if len(residue.name) > 3:
                        residue.name = "LIG"
    structure.write_pdb(str(destination))
    return destination


class Runner:
    def __init__(self, log_path: Path):
        self.log_path = log_path
        self.commands: list[list[str]] = []

    def run(
        self,
        command: list[str],
        *,
        cwd: Path | None = None,
        env_overrides: dict[str, str] | None = None,
    ) -> None:
        command = [str(value) for value in command]
        self.commands.append(command)
        with self.log_path.open("a", encoding="utf-8") as log:
            log.write(f"\n[{datetime.now(timezone.utc).isoformat()}] $ {shlex.join(command)}\n")
            log.flush()
            completed = subprocess.run(
                command,
                cwd=str(cwd) if cwd else None,
                env={
                    **os.environ,
                    "PYTHONHASHSEED": os.environ.get("PYTHONHASHSEED", "0"),
                    **(env_overrides or {}),
                },
                stdout=log,
                stderr=subprocess.STDOUT,
                text=True,
                check=False,
            )
        if completed.returncode:
            raise RuntimeError(
                f"Command failed with exit code {completed.returncode}: {shlex.join(command)}. "
                f"See {self.log_path}"
            )


def _read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def _number(row: dict[str, str], key: str) -> float:
    value = row.get(key, "")
    return float(value) if value not in {"", None} else 0.0


def _residue_contact_present(row: dict[str, str], prefix: str) -> bool:
    """Collapse atom annotations to one structural residue-pair observation."""
    return (
        _number(row, f"{prefix}_total_contact_count") > 0
        or _number(row, f"{prefix}_count_rosetta_hbond") > 0
    )


def _mutation_interpretations(
    substitutions: list[dict[str, str]],
    mutation_effects: list[dict[str, object]],
    network_summaries: list[dict[str, str]],
    network_changes: list[dict[str, str]],
    interface_delta: float,
) -> list[dict[str, object]]:
    """Create conservative pose-level evidence statements for each mutation."""
    effects_by_mutation = {str(row["substitution"]): row for row in mutation_effects}
    summaries_by_mutation = {row["mutation"]: row for row in network_summaries}
    interpretations: list[dict[str, object]] = []

    def pair_evidence(row: dict[str, str]) -> str:
        partner = row.get("mutant_partner") or row.get("parent_partner") or "unknown partner"
        old_types = row.get("parent_interaction_types", "none")
        new_types = row.get("mutant_interaction_types", "none")
        return (
            f"{partner} ({old_types} -> {new_types}; "
            f"delta {_number(row, 'delta_total'):+.3f} REU)"
        )

    for substitution in substitutions:
        mutation = substitution["substitution"]
        effect = effects_by_mutation[mutation]
        network = summaries_by_mutation.get(mutation, {})
        local_delta = _number(network, "sum_nonbonded_pair_energy_delta_REU" if "sum_nonbonded_pair_energy_delta_REU" in network else "sum_pair_energy_delta_REU")
        ligand_delta = float(effect["direct_ligand_delta_REU"])
        parent_contacts = float(effect["parent_contact_annotations"])
        mutant_contacts = float(effect["mutant_contact_annotations"])
        has_ligand_edge = (
            parent_contacts > 0 or mutant_contacts > 0
            or abs(float(effect["direct_ligand_parent_REU"])) >= 0.1
            or abs(float(effect["direct_ligand_mutant_REU"])) >= 0.1
        )

        if has_ligand_edge and ligand_delta <= -0.25 and local_delta <= 0.5:
            assessment = "Supported individually in this pose"
        elif has_ligand_edge and ligand_delta <= -0.25 and local_delta > 0.5:
            assessment = "Mixed: ligand interaction improves, local network worsens"
        elif has_ligand_edge and ligand_delta >= 0.25 and local_delta <= -0.5:
            assessment = "Mixed: local network improves, direct ligand interaction worsens"
        elif has_ligand_edge and ligand_delta >= 0.25 and local_delta >= 0.5:
            assessment = "Not supported individually in this pose"
        elif not has_ligand_edge and local_delta <= -0.5:
            assessment = "Supported indirectly through the local network"
        elif not has_ligand_edge and local_delta >= 0.5:
            assessment = "Not supported individually; possible combination-only role"
        else:
            assessment = "Context-dependent / energetically near-neutral"

        mutation_rows = [row for row in network_changes if row.get("mutation") == mutation]
        favorable_rows = sorted(
            (row for row in mutation_rows if _number(row, "delta_total") <= -0.25),
            key=lambda row: _number(row, "delta_total"),
        )
        unfavorable_rows = sorted(
            (row for row in mutation_rows if _number(row, "delta_total") >= 0.25),
            key=lambda row: _number(row, "delta_total"), reverse=True,
        )
        support = "; ".join(pair_evidence(row) for row in favorable_rows[:2]) or "No pair improved by at least 0.25 REU"
        opposition = "; ".join(pair_evidence(row) for row in unfavorable_rows[:2]) or "No pair worsened by at least 0.25 REU"

        same_class = substitution.get("same_primary_chemical_class") == "True"
        chemistry = (
            f"Conservative class change: {substitution['parent_chemical_class']} -> "
            f"{substitution['mutant_chemical_class']}"
            if same_class else
            f"Non-conservative change: {substitution['parent_chemical_class']} -> "
            f"{substitution['mutant_chemical_class']}"
        )
        rotamer_status = substitution.get("rotamer_analysis_status", "")
        if rotamer_status != "performed":
            rotamer = f"Not interpreted ({rotamer_status or 'not available'})"
        elif substitution.get("same_full_rotamer_bin") == "True":
            rotamer = "Same coarse canonical-chi bins (not a Dunbrack assignment)"
        elif substitution.get("same_common_chi_bins") == "True":
            rotamer = "Same bins for comparable chi angles; full bins are not directly comparable"
        elif int(float(substitution.get("n_comparable_chi") or 0)) == 0:
            rotamer = "No comparable chi angles (one residue lacks the corresponding side-chain torsion)"
        else:
            rotamer = "Different coarse canonical-chi bins (not a Dunbrack assignment)"

        ligand_statement = (
            f"{float(effect['direct_ligand_parent_REU']):.3f} -> "
            f"{float(effect['direct_ligand_mutant_REU']):.3f} REU "
            f"(delta {ligand_delta:+.3f}); contact annotations "
            f"{parent_contacts:.0f} -> {mutant_contacts:.0f}"
            if has_ligand_edge else "No direct ligand edge detected"
        )
        context = (
            f"The complete multi-mutant interface changes by {interface_delta:+.3f} REU. "
            "This combined score is favorable but is non-additive and cannot be assigned to this mutation alone."
            if interface_delta < 0 else
            f"The complete multi-mutant interface changes by {interface_delta:+.3f} REU; "
            "the combined pose does not provide favorable global support."
        )
        interpretations.append({
            "mutation": mutation,
            "antibody_number": substitution["parent_antibody_number"],
            "region": substitution["parent_region"],
            "assessment": assessment,
            "direct_ligand_evidence": ligand_statement,
            "local_network_delta_REU": local_delta,
            "strongest_supporting_pair_changes": support,
            "strongest_opposing_pair_changes": opposition,
            "chemistry_interpretation": chemistry,
            "rotamer_interpretation": rotamer,
            "cdr_global_fit_rmsd_A": _number(substitution, "cdr_global_fit_rmsd_A"),
            "combination_context": context,
            "evidence_limit": "Single relaxed pose; not a causal single-mutant or ensemble free-energy calculation",
        })
    return interpretations


def _write_automatic_summary(
    output: Path,
    estrogen_architecture: dict[str, dict[str, object]] | None = None,
    plip_summary: dict[str, dict[str, object]] | None = None,
) -> None:
    final = output / "final"
    settings = json.loads((final / "summary.json").read_text(encoding="utf-8"))
    metric_rows = _read_csv(final / "interface_metrics_comparison.csv")
    metrics = {row["metric"]: row for row in metric_rows}
    rmsd_rows = _read_csv(final / "structure_rmsd.csv")
    rmsds = {row["region"]: row for row in rmsd_rows}
    substitutions = _read_csv(final / "substitutions.csv")
    pairs = _read_csv(final / "pair_comparison.csv")
    network_changes_path = final / "mutation_network_changes.csv"
    network_summary_path = final / "mutation_network_summary.csv"
    network_type_energy_path = final / "mutation_network_type_energy_summary.csv"
    network_changes = _read_csv(network_changes_path) if network_changes_path.exists() else []
    network_summaries = _read_csv(network_summary_path) if network_summary_path.exists() else []
    network_type_energies = (
        _read_csv(network_type_energy_path) if network_type_energy_path.exists() else []
    )
    if settings.get("analysis_status", {}).get("network") == "not_run":
        network_changes = network_summaries = network_type_energies = []
    ligand_validation_path = output / "ligand_atom_validation.json"
    ligand_validation = (
        json.loads(ligand_validation_path.read_text(encoding="utf-8"))
        if ligand_validation_path.exists() else {}
    )
    ligand_pairs = [
        row for row in pairs
        if row.get("parent_pair_scope") == "aa_to_ligand"
        or row.get("mutant_pair_scope") == "aa_to_ligand"
    ]
    interface_names = [key for key in metrics if key.startswith("interface_delta_")]
    if len(interface_names) != 1:
        raise ValueError("Expected one ligand-interface metric; found " + str(interface_names))
    interface_name = interface_names[0]
    ligand_chain = interface_name.removeprefix("interface_delta_")
    interface = metrics[interface_name]
    mutation_effects = []
    for substitution in substitutions:
        label = substitution["parent_residue"]
        partner_rows = [
            row for row in pairs
            if label in {
                row.get("parent_residue1"), row.get("parent_residue2"),
                row.get("mutant_residue1"), row.get("mutant_residue2"),
            }
        ]
        direct = [
            row for row in partner_rows
            if row.get("parent_pair_scope") == "aa_to_ligand"
            or row.get("mutant_pair_scope") == "aa_to_ligand"
        ]
        mutation_effects.append({
            "substitution": substitution["substitution"],
            "antibody_number": substitution["parent_antibody_number"],
            "all_partner_delta_REU": sum(_number(row, "delta_total") for row in partner_rows),
            "direct_ligand_parent_REU": sum(_number(row, "parent_total") for row in direct),
            "direct_ligand_mutant_REU": sum(_number(row, "mutant_total") for row in direct),
            "direct_ligand_delta_REU": sum(_number(row, "delta_total") for row in direct),
            "parent_contact_annotations": sum(
                _number(row, "parent_total_contact_count") for row in direct
            ),
            "mutant_contact_annotations": sum(
                _number(row, "mutant_total_contact_count") for row in direct
            ),
            "parent_ligand_residue_pair_count": sum(
                _residue_contact_present(row, "parent") for row in direct
            ),
            "mutant_ligand_residue_pair_count": sum(
                _residue_contact_present(row, "mutant") for row in direct
            ),
            "parent_interaction_types": sorted({
                name for row in direct for name in row.get("parent_interaction_types", "").split(";")
                if name and name != "0"
            }),
            "mutant_interaction_types": sorted({
                name for row in direct for name in row.get("mutant_interaction_types", "").split(";")
                if name and name != "0"
            }),
        })

    bound_total_delta = float(metrics["bound_total"]["delta_mutant_minus_parent"])
    interface_rep_delta = float(metrics.get(f"if_{ligand_chain}_fa_rep", {}).get("delta_mutant_minus_parent", 0.0))
    comparison_quality_warning = abs(bound_total_delta) > 100.0 or interface_rep_delta > 10.0
    mutation_interpretations = [] if settings.get("analysis_status", {}).get("network") == "not_run" else _mutation_interpretations(
        substitutions,
        mutation_effects,
        network_summaries,
        network_changes,
        float(interface["delta_mutant_minus_parent"]),
    )
    for row in mutation_interpretations:
        row["comparison_quality"] = (
            "LOW: large total-energy/repulsive mismatch; matched relaxation is required"
            if comparison_quality_warning else "No large total-energy or interface-repulsion mismatch detected"
        )
    write_csv(final / "mutation_interpretation_summary.csv", mutation_interpretations)

    top_ligand_changes = sorted(
        ligand_pairs, key=lambda row: abs(_number(row, "delta_total")), reverse=True
    )[:15]
    ligand_hbonds: dict[str, list[dict[str, str]] | None] = {}
    for label in ("parent", "mutant"):
        path = output / f"{label}_hbonds.csv"
        ligand_hbonds[label] = [
            row for row in _read_csv(path)
            if row.get("donor_residue", "").startswith(f"{ligand_chain}:")
            or row.get("acceptor_residue", "").startswith(f"{ligand_chain}:")
        ] if path.exists() else None

    summary = {
        "analysis_status": settings.get("analysis_status", {}),
        "network_annotation_sources": [] if settings.get("analysis_status", {}).get("network") == "not_run" else ["Rosetta pair energies", "Rosetta H-bonds"] + ([] if settings.get("analysis_status", {}).get("arpeggio") == "not_run" else ["Arpeggio contacts"]),
        "ligand_chain": ligand_chain,
        "ligand_resname": ((plip_summary or {}).get("parent") or {}).get("ligand_resname"),
        "interface_parent_REU": float(interface["parent"]),
        "interface_mutant_REU": float(interface["mutant"]),
        "interface_delta_mutant_minus_parent_REU": float(interface["delta_mutant_minus_parent"]),
        "bound_total_delta_mutant_minus_parent_REU": bound_total_delta,
        "interface_rep_delta_mutant_minus_parent_REU": interface_rep_delta,
        "comparison_quality_warning": comparison_quality_warning,
        "global_fv_rmsd_A": float(rmsds["GLOBAL_FV"]["rmsd_after_global_fit_A"]),
        "framework_rmsd_A": float(rmsds["FRAMEWORK"]["rmsd_after_global_fit_A"]),
        "h_cdr3_rmsd_A": _number(rmsds.get("H_CDR3", {}), "rmsd_after_global_fit_A"),
        "h_cdr3_anchors_rmsd_A": float(rmsds['H_CDR3_ANCHORS']['rmsd_after_global_fit_A']) if rmsds.get('H_CDR3_ANCHORS', {}).get('rmsd_after_global_fit_A') not in ('',None) else None,
        "h_cdr3_tip_rmsd_A": float(rmsds['H_CDR3_TIP']['rmsd_after_global_fit_A']) if rmsds.get('H_CDR3_TIP', {}).get('rmsd_after_global_fit_A') not in ('',None) else None,
        "parent_ligand_contact_annotations": sum(
            _number(row, "parent_total_contact_count") for row in ligand_pairs
        ),
        "mutant_ligand_contact_annotations": sum(
            _number(row, "mutant_total_contact_count") for row in ligand_pairs
        ),
        "parent_ligand_residue_pair_count": sum(
            _residue_contact_present(row, "parent") for row in ligand_pairs
        ),
        "mutant_ligand_residue_pair_count": sum(
            _residue_contact_present(row, "mutant") for row in ligand_pairs
        ),
        "parent_detected_ligand_hbonds": len(ligand_hbonds["parent"]) if ligand_hbonds["parent"] is not None else None,
        "mutant_detected_ligand_hbonds": len(ligand_hbonds["mutant"]) if ligand_hbonds["mutant"] is not None else None,
        "parent_arpeggio_strict_ligand_hbond_annotations": sum(
            _number(row, "parent_count_hbond") for row in ligand_pairs
        ),
        "mutant_arpeggio_strict_ligand_hbond_annotations": sum(
            _number(row, "mutant_count_hbond") for row in ligand_pairs
        ),
        "substitution_count": len(substitutions),
        "substitutions": substitutions,
        "mutation_effects": mutation_effects,
        "rmsds": rmsd_rows,
        "interface_terms": metric_rows,
        "ligand_hbonds": ligand_hbonds,
        "mutation_network_summary": network_summaries,
        "mutation_network_changes": network_changes,
        "mutation_network_type_energy_summary": network_type_energies,
        "mutation_interpretations": mutation_interpretations,
        "ligand_atom_validation": ligand_validation,
        "estrogen_binding_architecture": estrogen_architecture or {},
        "plip_ligand_interactions": plip_summary or {},
    }
    if settings.get("analysis_status", {}).get("arpeggio") == "not_run":
        summary["parent_arpeggio_strict_ligand_hbond_annotations"] = None
        summary["mutant_arpeggio_strict_ligand_hbond_annotations"] = None
    (output / "automatic_summary.json").write_text(
        json.dumps(summary, indent=2) + "\n", encoding="utf-8"
    )

    direction = "more favorable" if summary["interface_delta_mutant_minus_parent_REU"] < 0 else "less favorable"
    lines = [
        "# Automated parent-versus-mutant summary",
        "",
        f"The mutant interface is **{direction}** in this pose comparison: "
        f"{summary['interface_parent_REU']:.3f} to {summary['interface_mutant_REU']:.3f} REU "
        f"(delta mutant-parent {summary['interface_delta_mutant_minus_parent_REU']:+.3f} REU).",
        "",
        f"Global Fv RMSD: **{summary['global_fv_rmsd_A']:.3f} A**.  "
        f"Framework RMSD: **{summary['framework_rmsd_A']:.3f} A**.",
        "",
        f"Ligand-contacting residue pairs: {summary['parent_ligand_residue_pair_count']:.0f} -> "
        f"{summary['mutant_ligand_residue_pair_count']:.0f}.  Atom-level ligand contact annotations: "
        f"{summary['parent_ligand_contact_annotations']:.0f} -> "
        f"{summary['mutant_ligand_contact_annotations']:.0f}.  Rosetta-scored ligand H-bonds: "
        f"{summary['parent_detected_ligand_hbonds'] if summary['parent_detected_ligand_hbonds'] is not None else 'NaN'} -> "
        f"{summary['mutant_detected_ligand_hbonds'] if summary['mutant_detected_ligand_hbonds'] is not None else 'NaN'}.  Arpeggio strict geometric H-bond "
        f"annotations: {summary['parent_arpeggio_strict_ligand_hbond_annotations'] if summary['parent_arpeggio_strict_ligand_hbond_annotations'] is not None else 'not run'} -> "
        f"{summary['mutant_arpeggio_strict_ligand_hbond_annotations'] if summary['mutant_arpeggio_strict_ligand_hbond_annotations'] is not None else 'not run'}.",
        "",
    ]
    omitted = [name for name, status in settings.get("analysis_status", {}).items() if status == "not_run"]
    if omitted:
        lines += ["Analyses not run by configuration: " + ", ".join(omitted) + ".", ""]
    if estrogen_architecture:
        parent_e2 = estrogen_architecture["parent"]
        mutant_e2 = estrogen_architecture["mutant"]
        lines += [
            "## Conserved E2-binding architecture", "",
            "The four features use DS conventional O3/O17 H-bonds, calibrated DS C18 packing, "
            "and the paper-reference-calibrated two-sided A-ring enclosure geometry. "
            "All four must occur in the same pose. DS A-ring pi-mediated contacts "
            "(stacked/T-shaped pi-pi, pi-cation, pi-anion) are a separate descriptor, not a pass requirement.", "",
            "DS-dependent criteria are unavailable until the separate Discovery Studio stage supplies its results.", "",
            "| Pose | Score | DS O3 H-bond | DS O17 H-bond | DS C18 packing | Steric sandwich | DS pi-mediated residues | All four features |",
            "|---|---:|---|---|---|---|---|---|",
            f"| Parent | {parent_e2['e2_conserved_architecture_score_0_to_4']}/4 | "
            f"{parent_e2['e2_o3_hbond']} | {parent_e2['e2_o17_hbond']} | "
            f"{parent_e2['e2_c18_hydrophobic_contact']} | {parent_e2['e2_geometric_sandwich']} | "
            f"{parent_e2.get('e2_a_ring_pi_contact_residues', 'unavailable')} | "
            f"{parent_e2['e2_strict_screen_pass'] if parent_e2['e2_strict_screen_pass'] is not None else 'unavailable'} |",
            f"| Mutant | {mutant_e2['e2_conserved_architecture_score_0_to_4']}/4 | "
            f"{mutant_e2['e2_o3_hbond']} | {mutant_e2['e2_o17_hbond']} | "
            f"{mutant_e2['e2_c18_hydrophobic_contact']} | {mutant_e2['e2_geometric_sandwich']} | "
            f"{mutant_e2.get('e2_a_ring_pi_contact_residues', 'unavailable')} | "
            f"{mutant_e2['e2_strict_screen_pass'] if mutant_e2['e2_strict_screen_pass'] is not None else 'unavailable'} |", "",
        ]
    if plip_summary:
        parent_plip = plip_summary["parent"]
        mutant_plip = plip_summary["mutant"]
        parent_counts = parent_plip["interaction_counts"]
        mutant_counts = mutant_plip["interaction_counts"]
        lines += [
            "## PLIP independent ligand-interaction cross-check", "",
            "PLIP was run on the same prepared, explicitly protonated structures with "
            "`--nohydro`. Its definitions are independent of Arpeggio and Rosetta; "
            "agreement is stronger evidence, while disagreement should be inspected geometrically.", "",
            "| Interaction | Parent | Mutant |",
            "|---|---:|---:|",
        ]
        for name in parent_counts:
            label = name.replace("_", " ").title()
            lines.append(f"| {label} | {parent_counts[name]} | {mutant_counts[name]} |")
        lines.append("")
    if comparison_quality_warning:
        lines += [
            "> **Comparability warning:** the structures have a very large bound-total energy "
            f"difference ({bound_total_delta:+.3f} REU) and/or interface-repulsion change "
            f"({interface_rep_delta:+.3f} REU). This usually indicates an unrelaxed or clashing "
            "pose. Do not attribute the apparent penalties causally to the substitutions until parent "
            "and mutant have undergone the same ligand-aware constrained relaxation protocol.",
            "",
        ]
    if ligand_validation:
        parent_validation = ligand_validation.get("parent", {})
        mutant_validation = ligand_validation.get("mutant", {})
        lines += [
            f"Ligand atom validation: parent **{parent_validation.get('status', 'unknown')}** "
            f"({parent_validation.get('explicit_hydrogen_count', 0)}/"
            f"{parent_validation.get('expected_hydrogen_count', 0)} explicit H); mutant "
            f"**{mutant_validation.get('status', 'unknown')}** "
            f"({mutant_validation.get('explicit_hydrogen_count', 0)}/"
            f"{mutant_validation.get('expected_hydrogen_count', 0)} explicit H).",
            "",
        ]
    lines += [
        "## Interface energy decomposition",
        "",
        "| Weighted Rosetta term | Parent | Mutant | Delta mutant-parent |",
        "|---|---:|---:|---:|",
    ]
    preferred_terms = [
        "interface_delta_X", "bound_total", "if_X_fa_atr", "if_X_fa_rep",
        "if_X_fa_sol", "if_X_fa_elec", "if_X_hbond_sc", "if_X_hbond_bb_sc",
        "if_X_hbond_lr_bb", "if_X_hbond_sr_bb", "if_X_lk_ball_wtd",
    ]
    preferred_terms = [name.replace("interface_delta_X", interface_name).replace("if_X_", f"if_{ligand_chain}_") for name in preferred_terms]
    for name in preferred_terms:
        row = metrics.get(name)
        if row:
            lines.append(
                f"| `{name}` | {float(row['parent']):.3f} | {float(row['mutant']):.3f} | "
                f"{float(row['delta_mutant_minus_parent']):+.3f} |"
            )
    lines += [
        "",
        "## Structural RMSDs",
        "",
        "All `RMSD after global fit` values use the same global Fv transform. The local-fit value "
        "reports internal shape change after fitting that region independently.",
        f"HCDR3 is also split using the parent tip sequence `{settings.get('hcdr3_tip_reference_sequence', 'unrecorded')}`; anchors are all other HCDR3 residues. Both use the same global transform. Anchor-only HCDR3 RMSD is the rotamer-analysis gate.",
        "",
        "| Region | Residues | RMSD after global fit (A) | Local-fit RMSD (A) |",
        "|---|---:|---:|---:|",
    ]
    for row in rmsd_rows:
        global_value = f"{float(row['rmsd_after_global_fit_A']):.3f}" if row['rmsd_after_global_fit_A'] else "NA"
        local_value = f"{float(row['local_fit_rmsd_A']):.3f}" if row['local_fit_rmsd_A'] else "NA"
        lines.append(
            f"| {row['region']} | {row['n_residues']} | "
            f"{global_value} | {local_value} |"
        )

    lines += [
        "",
        "## Substitutions, chemistry, and rotamers",
        "",
        "CDR rotamers are compared only when the gate RMSD after the global fit is at or below "
        f"the configured cutoff ({settings.get('rotamer_rmsd_cutoff_A', 'unrecorded')} A). HCDR3 uses anchor-only RMSD; other CDRs use their "
        "overall RMSD. Framework substitutions have no CDR RMSD gate. Missing anchor RMSD skips "
        "HCDR3 comparisons without falling back to whole-HCDR3 RMSD.",
        "",
        "| Substitution | Antibody number | Region | CDR RMSD | Rotamer gate region | Gate RMSD | Parent class | Mutant class | Same class | Rotamer status | Parent bin | Mutant bin | Same full bin |",
        "|---|---|---|---:|---|---:|---|---|---|---|---|---|---|",
    ]
    for row in substitutions:
        full_rotamer = row.get("same_full_rotamer_bin") or "N/A"
        cdr_rmsd = (
            f"{float(row['cdr_global_fit_rmsd_A']):.3f}"
            if row.get("cdr_global_fit_rmsd_A") else "N/A"
        )
        lines.append(
            f"| {row['substitution']} | {row['parent_antibody_number']} | {row['parent_region']} | "
            f"{cdr_rmsd} | {row.get('rotamer_gate_region') or 'N/A'} | "
            f"{row.get('rotamer_gate_rmsd_A', 'N/A')} | "
            f"{row['parent_chemical_class']} | {row['mutant_chemical_class']} | "
            f"{row['same_primary_chemical_class']} | {row['rotamer_analysis_status']} | "
            f"{row['parent_rotamer_bin'] or 'N/A'} | "
            f"{row['mutant_rotamer_bin']} | {full_rotamer} |"
        )

    lines += [
        "",
        "## Per-mutation energetic and contact effects",
        "",
        "`All-partner delta` sums residue-pair energy changes involving that position. Direct-ligand "
        "columns include only its pair with the ligand. `Residue pair present` is a 0/1 "
        "structural-contact indicator; `atom annotations` counts overlapping atom-level labels. "
        "Negative deltas are favorable.",
        "",
        "| Substitution | All-partner delta | Ligand parent | Ligand mutant | Ligand delta | Residue pair present parent->mutant | Atom annotations parent->mutant | Interaction types parent -> mutant |",
        "|---|---:|---:|---:|---:|---:|---:|---|",
    ]
    for row in mutation_effects:
        p_types = ", ".join(row["parent_interaction_types"]) or "none"
        m_types = ", ".join(row["mutant_interaction_types"]) or "none"
        lines.append(
            f"| {row['substitution']} ({row['antibody_number']}) | "
            f"{row['all_partner_delta_REU']:+.3f} | {row['direct_ligand_parent_REU']:.3f} | "
            f"{row['direct_ligand_mutant_REU']:.3f} | {row['direct_ligand_delta_REU']:+.3f} | "
            f"{row['parent_ligand_residue_pair_count']:.0f}->{row['mutant_ligand_residue_pair_count']:.0f} | "
            f"{row['parent_contact_annotations']:.0f}->{row['mutant_contact_annotations']:.0f} | "
            f"{p_types} -> {m_types} |"
        )

    lines += [
        "",
        "## Mutation-by-mutation evidence interpretation",
        "",
        "These are pose-level evidence assessments, not experimental proof. Negative energy deltas "
        "are favorable. `Local network` sums all scored residue-pair changes involving the mutated "
        "position, whereas `direct ligand` includes only that residue-ligand pair. Because all "
        f"{len(substitutions)} mutations occur in the same model, their effects are coupled and are not additive.",
        "",
        "| Mutation | Assessment | Direct ligand evidence | Local network delta | Chemistry and rotamer |",
        "|---|---|---|---:|---|",
    ]
    for row in mutation_interpretations:
        lines.append(
            f"| {row['mutation']} ({row['antibody_number']}) | **{row['assessment']}** | "
            f"{row['direct_ligand_evidence']} | {float(row['local_network_delta_REU']):+.3f} | "
            f"{row['chemistry_interpretation']}; {row['rotamer_interpretation']} |"
        )

    for row in mutation_interpretations:
        lines += [
            "",
            f"### Interpretation: {row['mutation']} ({row['antibody_number']})",
            "",
            f"- **Assessment:** {row['assessment']}.",
            f"- **Best supporting changes:** {row['strongest_supporting_pair_changes']}.",
            f"- **Main opposing changes:** {row['strongest_opposing_pair_changes']}.",
            f"- **Direct ligand:** {row['direct_ligand_evidence']}.",
            f"- **Structure/chemistry:** CDR global-fit RMSD {float(row['cdr_global_fit_rmsd_A']):.3f} A; "
            f"{row['chemistry_interpretation']}; {row['rotamer_interpretation']}.",
            f"- **Combination context:** {row['combination_context']}",
        ]

    lines += [
        "",
        "## Mutation-centred interaction-network changes",
        "",
        "A pair is included when it has a typed structural contact, a detected H-bond, or a "
        "pair-energy magnitude above the configured presence cutoff. `Formed` and `broken` describe "
        "whole residue-pair edges; the type columns show changes within retained edges.",
        "",
        "| Mutation | Partners | Formed pairs | Broken pairs | Retained changed | Retained stable | Formed types | Broken types | Sum pair-energy delta |",
        "|---|---:|---:|---:|---:|---:|---|---|---:|",
    ]
    for row in network_summaries:
        lines.append(
            f"| {row['mutation']} ({row['antibody_number']}) | {row['partner_count']} | "
            f"{row['formed_pairs']} | {row['broken_pairs']} | {row['retained_changed_pairs']} | "
            f"{row['retained_stable_pairs']} | {row['formed_interaction_types']} | "
            f"{row['broken_interaction_types']} | {float(row['sum_pair_energy_delta_REU']):+.3f} |"
        )

    lines += [
        "",
        "### Explicit residue pairs that formed or broke",
        "",
        "A **formed pair** means the residue pair did not meet the network-presence rule in the "
        "parent but does in the mutant. A **broken pair** means the reverse. `Structural contact/H-bond` "
        "means a contact type or Rosetta H-bond was detected. `Energy threshold only` means no geometric "
        f"contact type was assigned, but the absolute Rosetta pair energy crossed the configured {settings.get('network_energy_presence_cutoff_REU', 'unrecorded')} REU "
        "presence threshold; this is an energetic edge, not proof that a new physical contact formed.",
        "",
        "| Mutation | Residue pair | Event | Evidence | Scope | Old types | New types | Old energy/strength | New energy/strength | Delta REU | Pair charge old->new |",
        "|---|---|---|---|---|---|---|---:|---:|---:|---:|",
    ]
    formed_broken = [
        row for row in network_changes if row.get("network_status") in {"formed", "broken"}
    ]
    formed_broken.sort(key=lambda row: (
        str(row.get("mutation", "")), str(row.get("network_status", "")),
        str(row.get("mutant_partner") or row.get("parent_partner") or ""),
    ))
    for row in formed_broken:
        mutated = row.get("mutant_mutated_residue") or row.get("parent_mutated_residue") or "-"
        partner = row.get("mutant_partner") or row.get("parent_partner") or "-"
        basis = (
            "Energy threshold only"
            if row.get("event_basis") == "energy_threshold_only"
            else "Structural contact/H-bond"
        )
        old_strength = (
            f"{_number(row, 'parent_total'):.3f} / absent by network rule"
            if row["network_status"] == "formed"
            else f"{_number(row, 'parent_total'):.3f} / {row.get('parent_pair_energy_strength') or 'absent'}"
        )
        new_strength = (
            f"{_number(row, 'mutant_total'):.3f} / absent by network rule"
            if row["network_status"] == "broken"
            else f"{_number(row, 'mutant_total'):.3f} / {row.get('mutant_pair_energy_strength') or 'absent'}"
        )
        lines.append(
            f"| {row['mutation']} ({row['antibody_number']}) | {mutated} <-> {partner} | "
            f"**{row['network_status']} pair** | {basis} | {row.get('pair_scope', '')} | "
            f"{row['parent_interaction_types']} | {row['mutant_interaction_types']} | "
            f"{old_strength} | {new_strength} | {_number(row, 'delta_total'):+.3f} | "
            f"{_number(row, 'parent_pair_nominal_charge'):+.0f}->"
            f"{_number(row, 'mutant_pair_nominal_charge'):+.0f} |"
        )
    if not formed_broken:
        lines.append("| - | - | No formed or broken pairs | - | - | - | - | - | - | - | - |")

    lines += [
        "",
        "### Secondary aggregate by formed/broken interaction type",
        "",
        "This compact table names the residue pairs contributing to every formed or broken interaction "
        "type, including type changes within pairs retained in both structures. These are associated whole-pair Rosetta energies, not "
        "energies uniquely attributable to a geometric contact type. A multi-annotated pair therefore "
        "appears under each relevant type.",
        "",
        "| Mutation | Charge parent->mutant | Event/type | Contributing residue pairs | Pair count | Associated energy sum | Parent total | Mutant total | Delta total |",
        "|---|---:|---|---|---:|---:|---:|---:|---:|",
    ]
    for row in network_type_energies:
        lines.append(
            f"| {row['mutation']} ({row['antibody_number']}) | "
            f"{float(row['parent_mutation_nominal_charge']):+.0f}->"
            f"{float(row['mutant_mutation_nominal_charge']):+.0f} | "
            f"{row['event']} / {row['interaction_type']} | "
            f"{row.get('contributing_pairs') or 'see detailed pair table'} | {row['pair_count']} | "
            f"{float(row['sum_corresponding_pair_energy_REU']):.3f} | "
            f"{float(row['sum_parent_total_REU']):.3f} | "
            f"{float(row['sum_mutant_total_REU']):.3f} | "
            f"{float(row['sum_delta_total_REU']):+.3f} |"
        )

    for substitution in substitutions:
        mutation = substitution["substitution"]
        changed = [
            row for row in network_changes
            if row.get("mutation") == mutation and row.get("network_status") != "retained_stable"
        ]
        changed.sort(key=lambda row: abs(_number(row, "delta_total")), reverse=True)
        lines += ["", f"### {mutation} ({substitution['parent_antibody_number']})", ""]
        if not changed:
            lines += ["No formed, broken, or materially changed interaction edges.", ""]
            continue
        lines += [
            "| Partner | Pair charge old->new | Scope | Status | Old types | New types | Broken types | Formed types | Type/count changes | Old strength/REU | New strength/REU | Delta REU | Delta atr | Delta rep | Delta elec | Delta sol | Delta H-bond |",
            "|---|---:|---|---|---|---|---|---|---|---:|---:|---:|---:|---:|---:|---:|---:|",
        ]
        for row in changed:
            partner = row.get("mutant_partner") or row.get("parent_partner") or "-"
            old_strength = f"{row.get('parent_pair_energy_strength') or 'absent'} / {_number(row, 'parent_total'):.3f}"
            new_strength = f"{row.get('mutant_pair_energy_strength') or 'absent'} / {_number(row, 'mutant_total'):.3f}"
            delta_hbond = sum(
                _number(row, f"delta_{term}")
                for term in ("hbond_sc", "hbond_bb_sc", "hbond_sr_bb", "hbond_lr_bb")
            )
            lines.append(
                f"| {partner} | {_number(row, 'parent_pair_nominal_charge'):+.0f}->"
                f"{_number(row, 'mutant_pair_nominal_charge'):+.0f} | "
                f"{row.get('pair_scope', '')} | {row['network_status']} | "
                f"{row['parent_interaction_types']} | {row['mutant_interaction_types']} | "
                f"{row['broken_interaction_types']} | {row['formed_interaction_types']} | "
                f"{row['type_count_changes']} | {old_strength} | {new_strength} | "
                f"{_number(row, 'delta_total'):+.3f} | {_number(row, 'delta_fa_atr'):+.3f} | "
                f"{_number(row, 'delta_fa_rep'):+.3f} | {_number(row, 'delta_fa_elec'):+.3f} | "
                f"{_number(row, 'delta_fa_sol'):+.3f} | {delta_hbond:+.3f} |"
            )

    lines += ["", "## Rosetta-scored protein-ligand hydrogen bonds", ""]
    for label in ("parent", "mutant"):
        lines += [f"### {label.capitalize()}", ""]
        rows = ligand_hbonds[label]
        if rows is None:
            lines += ["Unavailable: Rosetta H-bond output missing.", ""]
            continue
        if not rows:
            lines += ["None detected.", ""]
            continue
        lines += [
            "| Donor | Donor atom | Acceptor | Acceptor atom | Detector energy |",
            "|---|---|---|---|---:|",
        ]
        for row in rows:
            lines.append(
                f"| {row['donor_residue']} {row['donor_residue_name']} | "
                f"{row['donor_atom']}/{row['donor_hydrogen']} | "
                f"{row['acceptor_residue']} {row['acceptor_residue_name']} | "
                f"{row['acceptor_atom']} | {float(row['hbond_energy']):.3f} |"
            )
        lines.append("")

    lines += [
        "## Largest ligand-pair energy changes",
        "",
        "| Parent residue | Mutant residue | Parent REU | Mutant REU | Delta | Residue pair present parent->mutant | Atom annotations parent->mutant | Interaction types parent -> mutant |",
        "|---|---|---:|---:|---:|---:|---:|---|",
    ]
    for row in top_ligand_changes:
        lines.append(
            f"| {row.get('parent_residue1') or row.get('parent_residue2') or '-'} | "
            f"{row.get('mutant_residue1') or row.get('mutant_residue2') or '-'} | "
            f"{_number(row, 'parent_total'):.3f} | {_number(row, 'mutant_total'):.3f} | "
            f"{_number(row, 'delta_total'):+.3f} | "
            f"{int(_residue_contact_present(row, 'parent'))}->{int(_residue_contact_present(row, 'mutant'))} | "
            f"{_number(row, 'parent_total_contact_count'):.0f}->{_number(row, 'mutant_total_contact_count'):.0f} | "
            f"{row.get('parent_interaction_types') or 'none'} -> "
            f"{row.get('mutant_interaction_types') or 'none'} |"
        )
    lines += [
        "",
        "## Interpretation limits",
        "",
        "- REU are Rosetta energy units, not experimental binding free energies.",
        "- Contact annotations can overlap at the atom level and are not unique residue contacts.",
        "- Detector H-bonds should be interpreted alongside the weighted Rosetta H-bond terms; a "
        "geometrically detected bond can be excluded or weighted differently by the score function.",
        "- Discovery Studio can report broader potential H-bond annotations, including carbon-H bonds "
        "and contacts accepted from heavy-atom geometry. These counts are not directly equivalent to "
        "Rosetta-scored H-bonds or Arpeggio's strict `hbond` classification.",
        "- Charges are nominal side-chain formal charges near neutral pH (D/E=-1, K/R=+1, "
        "histidine=0/titratable, E2L=0), not partial atomic charges or a pKa calculation.",
        "- Replicate relaxation or ensemble scoring is recommended before prioritizing small deltas.",
    ]
    (output / "automatic_summary.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


def run_workflow(args: argparse.Namespace) -> Path:
    started = time.time()
    output = Path(args.output).expanduser().resolve()
    if output.exists() and any(output.iterdir()):
        raise FileExistsError(f"Output directory is not empty; preserve existing results and choose a new directory: {output}")
    output.mkdir(parents=True, exist_ok=True)
    prepared = output / "prepared"
    prepared.mkdir(exist_ok=True)
    final = output / "final"
    arpeggio_parent = output / "arpeggio_parent"
    arpeggio_mutant = output / "arpeggio_mutant"
    plip_parent = output / "plip_parent"
    plip_mutant = output / "plip_mutant"
    arpeggio_parent.mkdir(exist_ok=True)
    arpeggio_mutant.mkdir(exist_ok=True)
    plip_parent.mkdir(exist_ok=True)
    plip_mutant.mkdir(exist_ok=True)
    log_path = output / "workflow.log"
    log_path.write_text("", encoding="utf-8")
    runner = Runner(log_path)

    prepare_script = resource_path("scripts/prepare_ligand_complexes.py")
    cif_script = resource_path("scripts/pdb_to_arpeggio_cif.py")
    parent_source = Path(args.parent).expanduser().resolve()
    mutant_source = Path(args.mutant).expanduser().resolve()
    parent_ligand_chain = args.parent_ligand_chain or args.ligand_chain
    mutant_ligand_chain = args.mutant_ligand_chain or args.ligand_chain
    same_complex = (
        parent_source == mutant_source
        and args.parent_chain_map == args.mutant_chain_map
        and parent_ligand_chain == mutant_ligand_chain
    )
    if not parent_source.is_file() or not mutant_source.is_file():
        raise FileNotFoundError(f"Missing parent or mutant structure: {parent_source}; {mutant_source}")
    parent = _convert_structure_input(
        parent_source, prepared / "parent_input.pdb",
        args.parent_chain_map, parent_ligand_chain,
    )
    if parent != prepared / "parent_input.pdb":
        shutil.copyfile(parent, prepared / "parent_input.pdb")
        parent = prepared / "parent_input.pdb"
    mutant = parent if same_complex else _convert_structure_input(
        mutant_source, prepared / "mutant_input.pdb",
        args.mutant_chain_map, mutant_ligand_chain,
    )
    if not same_complex and mutant != prepared / "mutant_input.pdb":
        shutil.copyfile(mutant, prepared / "mutant_input.pdb")
        mutant = prepared / "mutant_input.pdb"

    manifest = {
        "status": "running",
        "started_utc": datetime.now(timezone.utc).isoformat(),
        "parent": str(parent_source),
        "mutant": str(mutant_source),
        "mutant_input_sha256": sha256(mutant_source.read_bytes()).hexdigest(),
        "parent_analysis_input": str(parent),
        "mutant_analysis_input": str(mutant),
        "ligand_resname": args.ligand_resname,
        "output": str(output),
        "analysis_status": analysis_status(args, pending=True),
    }
    manifest_path = output / "workflow_manifest.json"
    manifest_path.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")

    try:
        from .ligand_chemistry import reference_smiles_from_source
        ligand_source = getattr(args, 'ligand_source', None)
        if bool(args.params) != bool(args.ligand_template_pdb):
            raise ValueError("Supply both --params and --ligand-template-pdb, or neither")
        if args.params:
            params = Path(args.params).expanduser().resolve()
            template = Path(args.ligand_template_pdb).expanduser().resolve()
        else:
            if not ligand_source:
                raise ValueError('New parameter generation requires --ligand-source SDF/MOL2')
            from .ligand_parameters import prepare_parameters
            definition = prepare_parameters(Path(ligand_source), prepared / 'ligand_definition', args.ligand_resname,
                                            args.rosetta_ligand_chain, args.molfile_to_params)
            params, template = Path(definition['params']), Path(definition['template_pdb'])

        ligand_smiles = getattr(args, 'ligand_smiles', None) or reference_smiles_from_source(
            ligand_source, params, template
        )
        chemistry_flags = ['--ligand-smiles', ligand_smiles] if ligand_smiles else []

        runner.run([
            sys.executable, str(prepare_script), "--reference", str(parent),
            "--complex", str(parent),
            "--generated-ligand-pdb", str(template),
            "--params", str(params),
            "--ligand-chain", parent_ligand_chain,
            "--output-ligand-chain", args.rosetta_ligand_chain,
            "--output", str(prepared),
            *chemistry_flags,
        ])
        if not same_complex:
            runner.run([
                sys.executable, str(prepare_script), "--reference", str(parent),
                "--complex", str(mutant),
                "--generated-ligand-pdb", str(template),
                "--params", str(params),
                "--ligand-chain", mutant_ligand_chain,
                "--output-ligand-chain", args.rosetta_ligand_chain,
                "--output", str(prepared),
                *chemistry_flags,
            ])
        parent_prepared = prepared / f"{parent.stem}_rosetta.pdb"
        mutant_prepared = prepared / f"{mutant.stem}_rosetta.pdb"

        parent_validation = validate_ligand_against_params(
            parent_prepared, params, args.rosetta_ligand_chain, template=template,
        )
        ligand_validation = {
            "parent": parent_validation,
            "mutant": parent_validation if same_complex else validate_ligand_against_params(
                mutant_prepared, params, args.rosetta_ligand_chain, template=template,
            ),
        }
        (output / "ligand_atom_validation.json").write_text(
            json.dumps(ligand_validation, indent=2) + "\n", encoding="utf-8"
        )
        invalid = [label for label, result in ligand_validation.items() if not result["valid"]]
        if invalid:
            raise ValueError(
                f"Ligand atom validation failed for {', '.join(invalid)}; "
                f"see {output / 'ligand_atom_validation.json'}"
            )

        common_score = ["--params", str(params), "--ligand-chain", args.rosetta_ligand_chain, "--weights", args.weights]
        if getattr(args, "legacy_pre_talaris", None) is not None:
            common_score.append("--legacy-pre-talaris" if args.legacy_pre_talaris else "--no-legacy-pre-talaris")
        for label, pdb in (("parent", parent_prepared), ("mutant", mutant_prepared)):
            if same_complex and label == "mutant":
                for suffix in ("pair_energies.sc", "interface.sc", "hbonds.csv"):
                    shutil.copyfile(output / f"parent_{suffix}", output / f"mutant_{suffix}")
                continue
            runner.run([
                sys.executable, "-m", "fv_energy_analysis.cli", "score-pyrosetta",
                "--pdb", str(pdb), "--pair-output", str(output / f"{label}_pair_energies.sc"),
                "--interface-output", str(output / f"{label}_interface.sc"),
                "--hbond-output", str(output / f"{label}_hbonds.csv"), *common_score,
            ])

        if not getattr(args, "skip_arpeggio", False):
            parent_cif, mutant_cif = prepared / "parent_arpeggio.cif", prepared / "mutant_arpeggio.cif"
            runner.run([sys.executable, str(cif_script), str(parent_prepared), str(parent_cif)])
            if same_complex:
                shutil.copyfile(parent_cif, mutant_cif)
            else:
                runner.run([sys.executable, str(cif_script), str(mutant_prepared), str(mutant_cif)])
            arpeggio = shutil.which("pdbe-arpeggio")
            if not arpeggio:
                raise RuntimeError("pdbe-arpeggio is not installed in the active environment")
            selection = None
            if args.arpeggio_scope == "neighbourhood":
                selection = mutation_neighbourhood_selection(
                    parent_prepared, mutant_prepared, args.fv_chains.split(","),
                    args.rosetta_ligand_chain, args.arpeggio_neighbourhood_radius,
                )
            manifest["arpeggio_selection"] = selection or {"scope": "full"}
            selection_options = {"selectors": selection["selectors"]} if selection else {}
            parent_json = run_arpeggio(runner, arpeggio, arpeggio_parent, parent_cif, **selection_options)
            if same_complex:
                mutant_json = arpeggio_mutant / "mutant_arpeggio.json"
                shutil.copyfile(parent_json, mutant_json)
            else:
                mutant_json = run_arpeggio(
                    runner, arpeggio, arpeggio_mutant, mutant_cif, **selection_options
                )
            if selection:
                for contacts in (parent_json, mutant_json):
                    filter_arpeggio_contacts(contacts, contacts, selection["selected_residues"])

        else:
            parent_json = mutant_json = None
            manifest["arpeggio_selection"] = {"scope": "not_run"}

        if not getattr(args, "skip_plip", False):
            plip = shutil.which("plip")
            if not plip:
                raise RuntimeError(
                    "PLIP is not installed in the active environment; run scripts/install_plip_wsl.sh once"
                )
            for label, pdb, directory in (
                ("parent", parent_prepared, plip_parent),
                ("mutant", mutant_prepared, plip_mutant),
            ):
                if same_complex and label == "mutant":
                    for suffix in ("xml", "txt"):
                        source = plip_parent / f"parent_plip.{suffix}"
                        if source.is_file():
                            shutil.copyfile(source, directory / f"mutant_plip.{suffix}")
                    continue
                runner.run([
                    plip, "-f", str(pdb), "-o", str(directory), "-x", "-t",
                    "--nohydro", "--name", f"{label}_plip",
                ])
            parent_plip = summarize_plip_xml(plip_parent / "parent_plip.xml", args.ligand_resname)
            plip_summary = {
                "parent": parent_plip,
                "mutant": parent_plip if same_complex else summarize_plip_xml(
                    plip_mutant / "mutant_plip.xml", args.ligand_resname
                ),
            }
            (output / "plip_ligand_interactions.json").write_text(
                json.dumps(plip_summary, indent=2) + "\n", encoding="utf-8"
            )

        else:
            plip_summary = None

        runner.run([
            sys.executable, "-m", "fv_energy_analysis.cli", "compare",
            "--parent", str(parent_prepared), "--mutant", str(mutant_prepared),
            "--fv-chains", args.fv_chains, "--ligand-chains", args.rosetta_ligand_chain,
            "--scheme", args.scheme, "--cdr-definition", args.cdr_definition,
            "--fit-region", args.fit_region, "--rmsd-atoms", args.rmsd_atoms,
            "--rotamer-rmsd-cutoff", str(args.rotamer_rmsd_cutoff),
            "--hcdr3-tip-sequence", getattr(args, "hcdr3_tip_sequence", "YYGSS"),
            "--network-energy-presence-cutoff", str(args.network_energy_presence_cutoff),
            "--network-energy-change-cutoff", str(args.network_energy_change_cutoff),
            "--parent-rosetta", str(output / "parent_pair_energies.sc"),
            "--mutant-rosetta", str(output / "mutant_pair_energies.sc"),
            *(["--parent-arpeggio", str(parent_json), "--mutant-arpeggio", str(mutant_json)] if parent_json is not None else []),
            *analysis_flags(args),
            "--parent-hbonds", str(output / "parent_hbonds.csv"),
            "--mutant-hbonds", str(output / "mutant_hbonds.csv"),
            "--parent-interface", str(output / "parent_interface.sc"),
            "--mutant-interface", str(output / "mutant_interface.sc"),
            "--output", str(final),
        ])
        record_comparison_status(final, args)
        fv_chain_set = set(args.fv_chains.split(","))
        if not getattr(args, "skip_e2_architecture", False):
            parent_architecture = analyze_estrogen_architecture(
                parent_prepared, parent_json,
                ligand_resname=args.ligand_resname, fv_chains=fv_chain_set, params_path=params, template_path=template,
            )
            estrogen_architecture = {
                "parent": parent_architecture,
                "mutant": parent_architecture if same_complex else analyze_estrogen_architecture(
                    mutant_prepared, mutant_json,
                    ligand_resname=args.ligand_resname, fv_chains=fv_chain_set,
                    params_path=params, template_path=template,
                ),
            }
        else:
            estrogen_architecture = omitted_e2_architecture()
        (output / "estrogen_binding_architecture.json").write_text(
            json.dumps(estrogen_architecture, indent=2) + "\n", encoding="utf-8"
        )
        _write_automatic_summary(output, None if getattr(args, "skip_e2_architecture", False) else estrogen_architecture, plip_summary)
        manifest.update({
            "status": "complete",
            "completed_utc": datetime.now(timezone.utc).isoformat(),
            "elapsed_seconds": round(time.time() - started, 3),
            "params": str(params),
            "ligand_template_pdb": str(template),
            "parent_artifacts_shared": same_complex,
            "scoring_protocol": scoring_protocol(params, args.weights, getattr(args, "legacy_pre_talaris", None)),
            "analysis_status": analysis_status(args),
            "commands": runner.commands,
            "summary": str(output / "automatic_summary.md"),
        })
    except Exception as exc:
        manifest.update({
            "status": "failed",
            "completed_utc": datetime.now(timezone.utc).isoformat(),
            "elapsed_seconds": round(time.time() - started, 3),
            "error": str(exc),
            "commands": runner.commands,
        })
        raise
    finally:
        manifest_path.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    return output


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--parent", required=True)
    parser.add_argument("--mutant", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--fv-chains", type=_chains, default="H,L")
    parser.add_argument("--ligand-chain", default="d", help="Ligand chain in the input PDBs")
    parser.add_argument("--parent-ligand-chain", help="Override the ligand chain for the parent input")
    parser.add_argument("--mutant-ligand-chain", help="Override the ligand chain for the mutant input")
    parser.add_argument(
        "--parent-chain-map", type=_chain_map, default={},
        help="Optional chain renaming before analysis, e.g. A:H,B:L",
    )
    parser.add_argument(
        "--mutant-chain-map", type=_chain_map, default={},
        help="Optional chain renaming before analysis, e.g. H:L,L:H for swapped CIF chains",
    )
    parser.add_argument("--rosetta-ligand-chain", default="X")
    parser.add_argument("--ligand-resname", default="E2L")
    parser.add_argument("--ligand-smiles", help="Explicit authoritative ligand identity; never inferred from a PDB")
    parser.add_argument("--ligand-source", help="User SDF or MOL2 defining chemical topology and stereochemistry")
    parser.add_argument("--scheme", default="imgt", choices=["imgt", "chothia", "kabat", "aho"])
    parser.add_argument("--cdr-definition", default="imgt", choices=["imgt", "chothia", "kabat", "north"])
    parser.add_argument("--fit-region", default="fv", choices=["fv", "framework"])
    parser.add_argument("--rmsd-atoms", default="backbone", choices=["backbone", "ca"])
    parser.add_argument("--rotamer-rmsd-cutoff", type=float, default=3.0)
    parser.add_argument("--hcdr3-tip-sequence", default="YYGSS")
    parser.add_argument("--network-energy-presence-cutoff", type=float, default=0.1)
    parser.add_argument("--network-energy-change-cutoff", type=float, default=0.25)
    parser.add_argument("--weights", default="ref2015")
    parser.add_argument("--legacy-pre-talaris", action=argparse.BooleanOptionalAction, default=None)
    add_analysis_arguments(parser)
    parser.add_argument("--arpeggio-scope", choices=["neighbourhood", "full"], default="neighbourhood")
    parser.add_argument("--arpeggio-neighbourhood-radius", type=float, default=6.0)
    parser.add_argument("--params", help="Reuse an existing ligand params file")
    parser.add_argument("--ligand-template-pdb", help="Template PDB generated with the supplied params")
    parser.add_argument("--molfile-to-params", help="Path to your licensed Rosetta molfile_to_params.py (or set FV_MOLFILE_TO_PARAMS)")
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    output = run_workflow(args)
    print(f"Analysis complete: {output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
