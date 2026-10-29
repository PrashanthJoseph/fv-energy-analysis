"""Batch parent-versus-mutant analysis with per-complex and master reports."""

from __future__ import annotations

import argparse
import concurrent.futures
import json
import math
import multiprocessing
import os
import re
import shutil
import subprocess
import sys
import time
from datetime import datetime, timezone
from hashlib import sha1, sha256
from pathlib import Path

from .estrogen_architecture import analyze_estrogen_architecture
from .report import write_csv
from .resources import resource_path
from .issues import record_issue
from .analysis_stages import add_analysis_arguments, analysis_flags, analysis_status, validate_analysis_reuse, validate_comparison_reuse
from .workflow import _chain_map, _chains
from .fast_workflow import run_cached_candidate, run_rebased_candidate, validate_scoring_reuse


DEFAULT_GLOBS = ("*_relaxed.pdb", "*_relaxed.cif", "*_relaxed.mmcif")
CDR_REGIONS = ("H_CDR1", "H_CDR2", "H_CDR3", "H_CDR3_ANCHORS", "H_CDR3_TIP", "L_CDR1", "L_CDR2", "L_CDR3")


def _safe_result_name(path: Path) -> str:
    name = re.sub(r"[^A-Za-z0-9._-]+", "_", path.stem).strip("._")
    return name or "complex"


def _discover(input_folder: Path, globs: list[str], recursive: bool) -> list[Path]:
    found: set[Path] = set()
    for pattern in globs or DEFAULT_GLOBS:
        iterator = input_folder.rglob(pattern) if recursive else input_folder.glob(pattern)
        found.update(path.resolve() for path in iterator if path.is_file())
    return sorted(found, key=lambda path: str(path).lower())


def _result_names(candidates: list[Path], input_folder: Path) -> dict[Path, str]:
    base_names = [_safe_result_name(path) for path in candidates]
    duplicate_names = {name for name in base_names if base_names.count(name) > 1}
    result: dict[Path, str] = {}
    for path, base in zip(candidates, base_names):
        if base in duplicate_names:
            relative = str(path.relative_to(input_folder.resolve())).encode("utf-8")
            base = f"{base}_{sha1(relative).hexdigest()[:8]}"
        result[path] = base
    return result


def _completed(result_dir: Path) -> bool:
    manifest = result_dir / "workflow_manifest.json"
    summary = result_dir / "automatic_summary.json"
    if not manifest.exists() or not summary.exists():
        return False
    try:
        return json.loads(manifest.read_text(encoding="utf-8")).get("status") == "complete"
    except (OSError, json.JSONDecodeError):
        return False


def _interaction_types(value: object) -> list[str]:
    return [
        item.strip().split(":", 1)[0]
        for item in str(value or "").split(";")
        if item.strip() not in {"", "none", "0"}
    ]


def _candidate_network_summary(
    candidate: str,
    result_dir: Path,
    changes: list[dict[str, object]],
) -> tuple[dict[str, object], list[dict[str, object]]]:
    """Summarize unique mutation-network pair/type events without double counting."""
    events: dict[tuple[str, str, str], dict[str, object]] = {}
    pair_records: dict[str, dict[str, object]] = {}
    pairs_by_event: dict[str, set[str]] = {"added": set(), "removed": set()}
    for change in changes:
        residue = str(change.get("parent_mutated_residue") or change.get("mutant_mutated_residue") or "")
        partner = str(change.get("parent_partner") or change.get("mutant_partner") or "")
        if not residue or not partner:
            continue
        pair = "<->".join(sorted((residue, partner)))
        pair_records.setdefault(pair, change)
        for event, field in (
            ("added", "formed_interaction_types"),
            ("removed", "broken_interaction_types"),
        ):
            for interaction_type in _interaction_types(change.get(field)):
                events.setdefault((event, interaction_type, pair), change)
                pairs_by_event[event].add(pair)

    type_rows: list[dict[str, object]] = []
    type_counts: dict[str, dict[str, int]] = {"added": {}, "removed": {}}
    for event in ("added", "removed"):
        interaction_types = sorted({key[1] for key in events if key[0] == event})
        for interaction_type in interaction_types:
            selected = {
                pair: change
                for (row_event, row_type, pair), change in events.items()
                if row_event == event and row_type == interaction_type
            }
            type_counts[event][interaction_type] = len(selected)
            type_rows.append({
                "candidate": candidate,
                "result_folder": str(result_dir),
                "event": event,
                "interaction_class": interaction_type,
                "unique_pair_count": len(selected),
                "contributing_pairs": ";".join(sorted(selected)),
                "sum_parent_pair_energy_REU": sum(_number(row.get("parent_total")) for row in selected.values()),
                "sum_mutant_pair_energy_REU": sum(_number(row.get("mutant_total")) for row in selected.values()),
                "sum_pair_energy_delta_REU": sum(_number(row.get("delta_total")) for row in selected.values()),
                "energy_overlap_warning": "Pair energy is repeated across classes when one pair has multiple annotations",
            })

    def pair_energy(event: str, field: str) -> float:
        return sum(_number(pair_records[pair].get(field)) for pair in pairs_by_event[event])

    changed_pairs = pairs_by_event["added"] | pairs_by_event["removed"] | {
        pair for pair, change in pair_records.items() if change.get("network_status") == "retained_changed"
    }
    compact = {
        event: ";".join(f"{name}:{count}" for name, count in sorted(type_counts[event].items())) or "none"
        for event in ("added", "removed")
    }
    return {
        "network_added_pair_type_events": sum(type_counts["added"].values()),
        "network_removed_pair_type_events": sum(type_counts["removed"].values()),
        "network_added_unique_pairs": len(pairs_by_event["added"]),
        "network_removed_unique_pairs": len(pairs_by_event["removed"]),
        "network_added_interaction_classes": compact["added"],
        "network_removed_interaction_classes": compact["removed"],
        "network_added_pairs_parent_energy_REU": pair_energy("added", "parent_total"),
        "network_added_pairs_mutant_energy_REU": pair_energy("added", "mutant_total"),
        "network_added_pairs_energy_delta_REU": pair_energy("added", "delta_total"),
        "network_removed_pairs_parent_energy_REU": pair_energy("removed", "parent_total"),
        "network_removed_pairs_mutant_energy_REU": pair_energy("removed", "mutant_total"),
        "network_removed_pairs_energy_delta_REU": pair_energy("removed", "delta_total"),
        "network_changed_unique_pairs": len(changed_pairs),
        "network_changed_unique_pairs_energy_delta_REU": sum(
            _number(pair_records[pair].get("delta_total")) for pair in changed_pairs
        ),
    }, type_rows


def _hbond_saltbridge_summary(
    candidate: str,
    result_dir: Path,
    changes: list[dict[str, object]],
) -> tuple[dict[str, object], list[dict[str, object]]]:
    """Summarize formed/broken H-bonds and salt bridges by chain scope.

    A residue pair is counted once per event/family/scope, even when both
    Arpeggio (``hbond``) and Rosetta (``rosetta_hbond``) annotate it.
    ``ionic`` is Arpeggio's ionic-contact/salt-bridge annotation.
    """
    scopes = ("intra_chain", "inter_chain", "protein_ligand")
    families = ("hbond", "salt_bridge")
    events = ("formed", "broken")
    grouped: dict[tuple[str, str, str], dict[str, dict[str, object]]] = {}

    def scope_for(change: dict[str, object], residue: str, partner: str) -> str:
        if str(change.get("pair_scope", "")) == "aa_to_ligand":
            return "protein_ligand"
        residue_chain = residue.split(":", 1)[0]
        partner_chain = partner.split(":", 1)[0]
        return "intra_chain" if residue_chain == partner_chain else "inter_chain"

    for change in changes:
        residue = str(change.get("parent_mutated_residue") or change.get("mutant_mutated_residue") or "")
        partner = str(change.get("parent_partner") or change.get("mutant_partner") or "")
        if not residue or not partner:
            continue
        pair = "<->".join(sorted((residue, partner)))
        scope = scope_for(change, residue, partner)
        for event, field in (("formed", "formed_interaction_types"), ("broken", "broken_interaction_types")):
            types = set(_interaction_types(change.get(field)))
            # A detector-label swap is not formation/loss of the interaction
            # family if either detector already identified that family.
            old_types = set(_interaction_types(change.get("parent_interaction_types")))
            new_types = set(_interaction_types(change.get("mutant_interaction_types")))
            family_sets = {"hbond": {"hbond", "rosetta_hbond"}, "salt_bridge": {"ionic"}}
            have_full_types = "parent_interaction_types" in change and "mutant_interaction_types" in change
            selected_families: list[tuple[str, list[str]]] = []
            hbond_sources = sorted(types & {"hbond", "rosetta_hbond"})
            if hbond_sources and (not have_full_types or not (old_types if event == "formed" else new_types) & family_sets["hbond"]):
                selected_families.append(("hbond", hbond_sources))
            if "ionic" in types and (not have_full_types or not (old_types if event == "formed" else new_types) & family_sets["salt_bridge"]):
                selected_families.append(("salt_bridge", ["ionic"]))
            for family, sources in selected_families:
                key = (event, family, scope)
                records = grouped.setdefault(key, {})
                record = records.setdefault(pair, dict(change))
                existing = set(str(record.get("source_interaction_types", "")).split(";")) - {""}
                record["source_interaction_types"] = ";".join(sorted(existing | set(sources)))

    columns: dict[str, object] = {}
    detail_rows: list[dict[str, object]] = []
    for event in events:
        for family in families:
            for scope in scopes:
                records = grouped.get((event, family, scope), {})
                prefix = f"{event}_{scope}_{family}"
                pairs = sorted(records)
                columns[f"{prefix}_pair_count"] = len(pairs)
                columns[f"{prefix}_pairs"] = ";".join(pairs) or "none"
                columns[f"{prefix}_energy_delta_REU"] = sum(
                    _number(record.get("delta_total")) for record in records.values()
                )
                for pair, record in sorted(records.items()):
                    detail_rows.append({
                        "candidate": candidate,
                        "result_folder": str(result_dir),
                        "event": event,
                        "interaction_family": family,
                        "scope": scope,
                        "pair": pair,
                        "mutation": record.get("mutation", ""),
                        "antibody_number": record.get("antibody_number", ""),
                        "source_interaction_types": record.get("source_interaction_types", ""),
                        "parent_pair_energy_REU": record.get("parent_total", ""),
                        "mutant_pair_energy_REU": record.get("mutant_total", ""),
                        "pair_energy_delta_REU": record.get("delta_total", ""),
                        "parent_pair_nominal_charge": record.get("parent_pair_nominal_charge", ""),
                        "mutant_pair_nominal_charge": record.get("mutant_pair_nominal_charge", ""),
                        "delta_pair_nominal_charge": record.get("delta_pair_nominal_charge", ""),
                    })
    return columns, detail_rows


def _summary_row(
    candidate: Path,
    result_dir: Path,
    status: str,
    error: str = "",
) -> tuple[dict[str, object], list[dict[str, object]], list[dict[str, object]], list[dict[str, object]]]:
    row: dict[str, object] = {
        "candidate": candidate.name,
        "candidate_path": str(candidate),
        "result_folder": str(result_dir),
        "status": status,
        "error": error,
    }
    if status != "complete":
        return row, [], [], []
    summary = json.loads((result_dir / "automatic_summary.json").read_text(encoding="utf-8"))
    substitutions = summary.get("substitutions", [])
    parent_architecture = summary.get("estrogen_binding_architecture", {}).get("parent", {})
    mutant_architecture = summary.get("estrogen_binding_architecture", {}).get("mutant", {})
    parent_plip_pi_stacks = summary.get("plip_ligand_interactions", {}).get("parent", {}).get("interaction_counts", {}).get("pi_stacks")
    mutant_plip_pi_stacks = summary.get("plip_ligand_interactions", {}).get("mutant", {}).get("interaction_counts", {}).get("pi_stacks")
    from .fingerprint_counts import _arpeggio_counts, _ligand_resname
    def arpeggio_ring_contact(role):
        if summary.get("analysis_status", {}).get("arpeggio") == "not_run":
            return None
        path = result_dir / f"arpeggio_{role}" / f"{role}_arpeggio.json"
        return _arpeggio_counts(path, _ligand_resname(summary, result_dir))["plane-plane"] > 0 if path.is_file() else None
    parent_arpeggio_pi_stack = arpeggio_ring_contact("parent")
    mutant_arpeggio_pi_stack = arpeggio_ring_contact("mutant")
    parent_plip_pi_stack = _number(parent_plip_pi_stacks) > 0 if parent_plip_pi_stacks is not None else None
    mutant_plip_pi_stack = _number(mutant_plip_pi_stacks) > 0 if mutant_plip_pi_stacks is not None else None
    row.update({
        "substitution_count": summary.get("substitution_count", len(substitutions)),
        "substitutions": ";".join(item.get("substitution", "") for item in substitutions),
        "interface_parent_REU": summary.get("interface_parent_REU"),
        "interface_mutant_REU": summary.get("interface_mutant_REU"),
        "interface_delta_REU": summary.get("interface_delta_mutant_minus_parent_REU"),
        "bound_total_delta_REU": summary.get("bound_total_delta_mutant_minus_parent_REU"),
        "interface_rep_delta_REU": summary.get("interface_rep_delta_mutant_minus_parent_REU"),
        "global_fv_rmsd_A": summary.get("global_fv_rmsd_A"),
        "framework_rmsd_A": summary.get("framework_rmsd_A"),
        "parent_ligand_contacts": summary.get("parent_ligand_contact_annotations"),
        "mutant_ligand_contacts": summary.get("mutant_ligand_contact_annotations"),
        "parent_ligand_residue_pairs": summary.get("parent_ligand_residue_pair_count"),
        "mutant_ligand_residue_pairs": summary.get("mutant_ligand_residue_pair_count"),
        "parent_rosetta_ligand_hbonds": summary.get("parent_detected_ligand_hbonds"),
        "mutant_rosetta_ligand_hbonds": summary.get("mutant_detected_ligand_hbonds"),
        "parent_arpeggio_ligand_hbonds": summary.get("parent_arpeggio_strict_ligand_hbond_annotations"),
        "mutant_arpeggio_ligand_hbonds": summary.get("mutant_arpeggio_strict_ligand_hbond_annotations"),
        "parent_plip_ligand_hbonds": summary.get("plip_ligand_interactions", {}).get("parent", {}).get("interaction_counts", {}).get("hydrogen_bonds"),
        "mutant_plip_ligand_hbonds": summary.get("plip_ligand_interactions", {}).get("mutant", {}).get("interaction_counts", {}).get("hydrogen_bonds"),
        "parent_arpeggio_pi_stack_detected": parent_arpeggio_pi_stack,
        "mutant_arpeggio_pi_stack_detected": mutant_arpeggio_pi_stack,
        "parent_arpeggio_pi_stack_sandwich": None if "e2_pi_stack_source" in parent_architecture else parent_architecture.get("e2_pi_stack_sandwich"),
        "mutant_arpeggio_pi_stack_sandwich": None if "e2_pi_stack_source" in mutant_architecture else mutant_architecture.get("e2_pi_stack_sandwich"),
        "parent_plip_pi_stacks": parent_plip_pi_stacks,
        "mutant_plip_pi_stacks": mutant_plip_pi_stacks,
        "parent_plip_pi_stack_detected": parent_plip_pi_stack,
        "mutant_plip_pi_stack_detected": mutant_plip_pi_stack,
        "parent_arpeggio_plip_pi_stack_consensus": None if parent_arpeggio_pi_stack is None or parent_plip_pi_stack is None else parent_arpeggio_pi_stack and parent_plip_pi_stack,
        "mutant_arpeggio_plip_pi_stack_consensus": None if mutant_arpeggio_pi_stack is None or mutant_plip_pi_stack is None else mutant_arpeggio_pi_stack and mutant_plip_pi_stack,
        "comparison_quality_warning": summary.get("comparison_quality_warning", False),
    })
    ligand_validation = summary.get("ligand_atom_validation", {})
    for label in ("parent", "mutant"):
        validation = ligand_validation.get(label, {})
        row[f"{label}_ligand_validation_status"] = validation.get("status", "not_recorded")
        row[f"{label}_explicit_ligand_hydrogens"] = validation.get("explicit_hydrogen_count", "")
        row[f"{label}_expected_ligand_hydrogens"] = validation.get("expected_hydrogen_count", "")
    rmsds = {item.get("region"): item for item in summary.get("rmsds", [])}
    for region in CDR_REGIONS:
        value = rmsds.get(region, {}).get("rmsd_after_global_fit_A", "")
        row[f"{region.lower()}_rmsd_A"] = value
    interpretations: list[dict[str, object]] = []
    assessments = []
    for item in summary.get("mutation_interpretations", []):
        combined = {"candidate": candidate.name, "result_folder": str(result_dir), **item}
        interpretations.append(combined)
        assessments.append(f"{item.get('mutation')}: {item.get('assessment')}")
    row["mutation_assessments"] = "; ".join(assessments)
    saved_architecture = mutant_architecture
    architecture_pdb = result_dir / "prepared" / f"{candidate.stem}_input_rosetta.pdb"
    architecture_json = result_dir / "arpeggio_mutant" / "mutant_arpeggio.json"
    if saved_architecture:
        row.update(saved_architecture)
        row["e2_metric_status"] = "not_run" if saved_architecture.get("e2_architecture_status") == "not_run" else "complete"
        row["e2_metric_error"] = ""
    elif architecture_pdb.is_file() and architecture_json.is_file():
        try:
            from .fingerprint_reports import _saved_architecture
            architecture = _saved_architecture(result_dir, "mutant", result_dir / "discovery_studio_interactions.csv")
            row.update(architecture)
            row["e2_metric_status"] = "complete"
            row["e2_metric_error"] = ""
            (result_dir / "estrogen_binding_architecture.json").write_text(
                json.dumps(architecture, indent=2) + "\n", encoding="utf-8"
            )
        except Exception as exc:
            row["e2_metric_status"] = "failed"
            row["e2_metric_error"] = str(exc)
    else:
        row["e2_metric_status"] = "not_available"
        row["e2_metric_error"] = "Prepared mutant PDB or mutant Arpeggio JSON is missing"
    network_columns, network_type_rows = _candidate_network_summary(
        candidate.name,
        result_dir,
        summary.get("mutation_network_changes", []),
    )
    row.update(network_columns)
    bond_columns, bond_rows = _hbond_saltbridge_summary(
        candidate.name,
        result_dir,
        summary.get("mutation_network_changes", []),
    )
    row.update(bond_columns)
    stages = summary.get("analysis_status", {})
    if stages.get("network") == "not_run":
        for key in [*network_columns, *bond_columns]:
            row[key] = None
        network_type_rows = bond_rows = []
    elif stages.get("arpeggio") == "not_run":
        for key in bond_columns:
            if "salt_bridge" in key:
                row[key] = None
    row["network_analysis_status"] = stages.get("network", "complete")
    row["network_annotation_sources"] = ";".join(summary.get("network_annotation_sources", []))
    return row, interpretations, network_type_rows, bond_rows


def _number(value: object) -> float:
    try:
        return float(value)
    except (TypeError, ValueError):
        return 0.0


def _write_master_reports(
    output: Path,
    rows: list[dict[str, object]],
    mutation_rows: list[dict[str, object]],
    network_rows: list[dict[str, object]] | None = None,
    bond_rows: list[dict[str, object]] | None = None,
) -> tuple[Path, Path, Path, Path, Path, Path, Path]:
    complete = [row for row in rows if row["status"] == "complete"]
    complete.sort(key=lambda row: _number(row.get("interface_delta_REU")))
    for rank, row in enumerate(complete, 1):
        row["rank_by_interface_delta"] = rank
    screen_order = sorted(
        complete,
        key=lambda row: (
            not bool(row.get("e2_strict_screen_pass")),
            -_number(row.get("e2_conserved_architecture_score_0_to_4")),
            _number(row.get("interface_delta_REU")),
        ),
    )
    for rank, row in enumerate(screen_order, 1):
        row["rank_by_e2_screen"] = rank
    failures = [row for row in rows if row["status"] != "complete"]
    ordered = complete + failures
    summary_csv = output / "batch_summary.csv"
    mutation_csv = output / "batch_mutation_summary.csv"
    network_csv = output / "batch_interaction_network_summary.csv"
    bond_csv = output / "batch_hbond_saltbridge_changes.csv"
    e2_screen_csv = output / "e2_architecture_screen.csv"
    e2_screen_rows = [{
        "rank_by_e2_screen": row.get("rank_by_e2_screen", ""),
        "rank_by_interface_delta": row.get("rank_by_interface_delta", ""),
        "candidate": row.get("candidate", ""),
        "substitutions": row.get("substitutions", ""),
        "e2_conserved_architecture_score_0_to_4": row.get("e2_conserved_architecture_score_0_to_4", ""),
        "e2_o3_hbond": row.get("e2_o3_hbond", ""),
        "e2_o3_hbond_residues": row.get("e2_o3_hbond_residues", ""),
        "e2_o17_hbond": row.get("e2_o17_hbond", ""),
        "e2_o17_hbond_residues": row.get("e2_o17_hbond_residues", ""),
        "e2_c18_hydrophobic_contact": row.get("e2_c18_hydrophobic_contact", ""),
        "e2_c18_hydrophobic_residues": row.get("e2_c18_hydrophobic_residues", ""),
        "e2_geometric_sandwich": row.get("e2_geometric_sandwich", ""),
        "e2_geometric_sandwich_positive_face": row.get("e2_geometric_sandwich_positive_face", ""),
        "e2_geometric_sandwich_negative_face": row.get("e2_geometric_sandwich_negative_face", ""),
        "e2_a_ring_pi_stack_count": row.get("e2_a_ring_pi_stack_count", ""),
        "e2_a_ring_pi_stack_residues": row.get("e2_a_ring_pi_stack_residues", ""),
        "mutant_arpeggio_pi_stack_detected": row.get("mutant_arpeggio_pi_stack_detected", ""),
        "mutant_arpeggio_pi_stack_sandwich": row.get("mutant_arpeggio_pi_stack_sandwich", ""),
        "mutant_plip_pi_stack_detected": row.get("mutant_plip_pi_stack_detected", ""),
        "mutant_arpeggio_plip_pi_stack_consensus": row.get("mutant_arpeggio_plip_pi_stack_consensus", ""),
        "e2_strict_screen_pass": row.get("e2_strict_screen_pass", ""),
        "interface_delta_REU": row.get("interface_delta_REU", ""),
        "bound_total_delta_REU": row.get("bound_total_delta_REU", ""),
        "comparison_quality_warning": row.get("comparison_quality_warning", ""),
    } for row in screen_order]
    suffix = ""
    try:
        write_csv(summary_csv, ordered)
        write_csv(mutation_csv, mutation_rows)
        write_csv(network_csv, network_rows or [])
        write_csv(bond_csv, bond_rows or [])
        write_csv(e2_screen_csv, e2_screen_rows)
    except PermissionError:
        # Spreadsheet applications on Windows commonly lock an open CSV.
        # Preserve the completed run by writing a timestamped replacement.
        suffix = datetime.now(timezone.utc).strftime("_%Y%m%dT%H%M%SZ")
        summary_csv = output / f"batch_summary{suffix}.csv"
        mutation_csv = output / f"batch_mutation_summary{suffix}.csv"
        network_csv = output / f"batch_interaction_network_summary{suffix}.csv"
        bond_csv = output / f"batch_hbond_saltbridge_changes{suffix}.csv"
        e2_screen_csv = output / f"e2_architecture_screen{suffix}.csv"
        write_csv(summary_csv, ordered)
        write_csv(mutation_csv, mutation_rows)
        write_csv(network_csv, network_rows or [])
        write_csv(bond_csv, bond_rows or [])
        write_csv(e2_screen_csv, e2_screen_rows)

    lines = [
        "# Batch Fv-ligand analysis summary", "",
        f"Completed: **{len(complete)}**. Failed: **{len(failures)}**. "
        "Ranking uses mutant-minus-parent Rosetta interface energy; more negative is favorable.", "",
        "A rank is a pose-level prioritization, not statistical significance. Comparisons marked with "
        "a quality warning require matched relaxation or further inspection.", "",
        "## Ranked complex summary", "",
        "| Rank | Complex | Substitutions | Interface delta | Bound-total delta | Global RMSD | H-CDR3 RMSD | L-CDR3 RMSD | Rosetta H-bonds parent->mutant | Arpeggio H-bonds parent->mutant | PLIP H-bonds parent->mutant | PLIP pi-stacks parent->mutant | Quality warning | Report |",
        "|---:|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---|---|",
    ]
    for row in complete:
        result_dir = Path(str(row["result_folder"]))
        report = (result_dir / "automatic_summary.md").relative_to(output).as_posix()
        lines.append(
            f"| {row['rank_by_interface_delta']} | {row['candidate']} | {row.get('substitutions', '') or 'none'} | "
            f"{_number(row.get('interface_delta_REU')):+.3f} | {_number(row.get('bound_total_delta_REU')):+.3f} | "
            f"{_number(row.get('global_fv_rmsd_A')):.3f} | {_number(row.get('h_cdr3_rmsd_A')):.3f} | "
            f"{_number(row.get('l_cdr3_rmsd_A')):.3f} | "
            f"{_number(row.get('parent_rosetta_ligand_hbonds')):.0f}->{_number(row.get('mutant_rosetta_ligand_hbonds')):.0f} | "
            f"{_number(row.get('parent_arpeggio_ligand_hbonds')):.0f}->{_number(row.get('mutant_arpeggio_ligand_hbonds')):.0f} | "
            f"{_number(row.get('parent_plip_ligand_hbonds')):.0f}->{_number(row.get('mutant_plip_ligand_hbonds')):.0f} | "
            f"{_number(row.get('parent_plip_pi_stacks')):.0f}->{_number(row.get('mutant_plip_pi_stacks')):.0f} | "
            f"{'YES' if row.get('comparison_quality_warning') else 'no'} | [individual report]({report}) |"
        )

    lines += [
        "", "## Conserved E2-binding architecture screen", "",
        "DS supplies O3/O17 conventional H-bonds, calibrated C18 packing and A-ring pi-mediated contacts. The steric sandwich is a separate calibrated coordinate proxy. All-four pass requires O3, O17, C18 and steric features in the same pose; pi-mediated contacts are not a mandatory fifth feature. No Arpeggio or PLIP calls contribute to this screen.", "",
        "| E2 rank | Complex | Score / 4 | DS O3 H-bond | DS O17 H-bond | DS C18 packing | Steric sandwich | DS pi-mediated residues | All four features | Interface delta |",
        "|---:|---|---:|---|---|---|---|---|---|---:|",
    ]
    for row in screen_order:
        lines.append(
            f"| {row.get('rank_by_e2_screen', '')} | {row['candidate']} | "
            f"{row.get('e2_conserved_architecture_score_0_to_4', 'unavailable')} | "
            f"{row.get('e2_o3_hbond', 'unavailable')} | {row.get('e2_o17_hbond', 'unavailable')} | "
            f"{row.get('e2_c18_hydrophobic_contact', 'unavailable')} | {row.get('e2_geometric_sandwich', 'unavailable')} | "
            f"{row.get('e2_a_ring_pi_contact_residues', 'unavailable')} | "
            f"{row.get('e2_strict_screen_pass', 'unavailable')} | {_number(row.get('interface_delta_REU')):+.3f} |"
        )

    lines += [
        "", "## Interaction-network changes", "",
        "Counts are unique residue-pair/interaction-class events around mutated residues; shared mutation pairs are counted once. Positive energy delta is unfavorable and negative is favorable.", "",
        "| Complex | Added events | Removed events | Added classes | Removed classes | Added-pair energy delta | Removed-pair energy delta | All changed-pair energy delta |",
        "|---|---:|---:|---|---|---:|---:|---:|",
    ]
    for row in complete:
        lines.append(
            f"| {row['candidate']} | {int(_number(row.get('network_added_pair_type_events')))} | "
            f"{int(_number(row.get('network_removed_pair_type_events')))} | "
            f"{row.get('network_added_interaction_classes', 'none')} | "
            f"{row.get('network_removed_interaction_classes', 'none')} | "
            f"{_number(row.get('network_added_pairs_energy_delta_REU')):+.3f} | "
            f"{_number(row.get('network_removed_pairs_energy_delta_REU')):+.3f} | "
            f"{_number(row.get('network_changed_unique_pairs_energy_delta_REU')):+.3f} |"
        )

    lines += [
        "", "## Formed and broken hydrogen bonds and salt bridges", "",
        "Counts are unique residue pairs around the mutations. H-bond pairs are deduplicated when both Rosetta and Arpeggio detect the same pair. Salt bridges use Arpeggio ionic annotations and remain protonation-dependent. Pair identities and energies are in `batch_hbond_saltbridge_changes.csv`.", "",
        "| Complex | H-bonds formed intra/inter/ligand | H-bonds broken intra/inter/ligand | Salt bridges formed intra/inter/ligand | Salt bridges broken intra/inter/ligand | H-bond formed energy delta | H-bond broken energy delta | Salt-bridge formed energy delta | Salt-bridge broken energy delta |",
        "|---|---:|---:|---:|---:|---:|---:|---:|---:|",
    ]
    for row in complete:
        def counts(event: str, family: str) -> str:
            return "/".join(str(int(_number(row.get(f"{event}_{scope}_{family}_pair_count")))) for scope in ("intra_chain", "inter_chain", "protein_ligand"))

        def energy(event: str, family: str) -> float:
            return sum(_number(row.get(f"{event}_{scope}_{family}_energy_delta_REU")) for scope in ("intra_chain", "inter_chain", "protein_ligand"))

        lines.append(
            f"| {row['candidate']} | {counts('formed', 'hbond')} | {counts('broken', 'hbond')} | "
            f"{counts('formed', 'salt_bridge')} | {counts('broken', 'salt_bridge')} | "
            f"{energy('formed', 'hbond'):+.3f} | {energy('broken', 'hbond'):+.3f} | "
            f"{energy('formed', 'salt_bridge'):+.3f} | {energy('broken', 'salt_bridge'):+.3f} |"
        )

    lines += ["", "## Mutation-level interpretation", "",
              "| Complex | Mutation | Region | Assessment | Direct ligand evidence | Local-network delta |",
              "|---|---|---|---|---|---:|"]
    for row in mutation_rows:
        lines.append(
            f"| {row.get('candidate', '')} | {row.get('mutation', '')} ({row.get('antibody_number', '')}) | "
            f"{row.get('region', '')} | {row.get('assessment', '')} | {row.get('direct_ligand_evidence', '')} | "
            f"{_number(row.get('local_network_delta_REU')):+.3f} |"
        )
    if failures:
        lines += ["", "## Failed analyses", "", "| Complex | Error | Result folder |", "|---|---|---|"]
        for row in failures:
            lines.append(f"| {row['candidate']} | {row.get('error', '')} | {row['result_folder']} |")
    summary_markdown = output / f"batch_summary{suffix}.md"
    try:
        summary_markdown.write_text("\n".join(lines) + "\n", encoding="utf-8")
    except PermissionError:
        suffix = datetime.now(timezone.utc).strftime("_%Y%m%dT%H%M%SZ")
        summary_markdown = output / f"batch_summary{suffix}.md"
        summary_markdown.write_text("\n".join(lines) + "\n", encoding="utf-8")
    definitions_source = resource_path("docs/batch_output_columns.md")
    definitions_output = output / "COLUMN_DEFINITIONS.md"
    shutil.copyfile(definitions_source, definitions_output)
    return summary_csv, mutation_csv, network_csv, bond_csv, e2_screen_csv, summary_markdown, definitions_output


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--parent", required=True)
    parser.add_argument("--input-folder", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--glob", action="append", default=[], help="Repeatable input glob; defaults to *_relaxed.pdb/cif/mmcif")
    parser.add_argument("--recursive", action="store_true")
    parser.add_argument("--force", action="store_true", help="Rerun analyses already marked complete")
    parser.add_argument(
        "--reuse-mutants-from",
        help=(
            "Existing completed batch whose mutant PyRosetta, Arpeggio, PLIP, and "
            "prepared-structure artifacts should be reused; candidates not present there run normally"
        ),
    )
    parser.add_argument(
        "--reuse-parent-cache",
        help="Reuse an existing complete parent cache after verifying parent, settings, code and artifacts",
    )
    parser.add_argument("--stop-on-error", action="store_true")
    parser.add_argument("--fv-chains", type=_chains, default="H,L")
    parser.add_argument("--ligand-chain", default="d")
    parser.add_argument("--parent-ligand-chain")
    parser.add_argument("--mutant-ligand-chain")
    parser.add_argument("--parent-chain-map", type=_chain_map, default={})
    parser.add_argument("--mutant-chain-map", type=_chain_map, default={})
    parser.add_argument("--rosetta-ligand-chain", default="X")
    parser.add_argument("--ligand-resname", default="E2L")
    parser.add_argument("--ligand-smiles")
    parser.add_argument("--ligand-source")
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
    parser.add_argument("--arpeggio-scope", choices=["neighbourhood", "full"], default="neighbourhood",
                        help="Mutation-neighbourhood queries (default) or full-complex geometric annotations")
    parser.add_argument("--arpeggio-neighbourhood-radius", type=float, default=6.0,
                        help="Heavy-atom mutation-neighbour radius in A, unioned across both states")
    parser.add_argument("--params")
    parser.add_argument("--ligand-template-pdb")
    parser.add_argument("--molfile-to-params", help="Path to your licensed Rosetta molfile_to_params.py (or set FV_MOLFILE_TO_PARAMS)")
    parser.add_argument(
        "--jobs", type=int, default=2,
        help="Persistent parallel PyRosetta workers (default: 2; raise only when RAM permits)",
    )
    return parser


def _map_arg(mapping: dict[str, str]) -> str:
    return ",".join(f"{source}:{target}" for source, target in mapping.items())


def _parent_cache_key(args: argparse.Namespace, parent: Path) -> str:
    payload = {
        "cache_schema": 3,
        "parent": str(parent), "parent_mtime_ns": parent.stat().st_mtime_ns,
        "parent_sha256": sha256(parent.read_bytes()).hexdigest(),
        "fv_chains": args.fv_chains, "ligand_chain": args.ligand_chain,
        "parent_ligand_chain": args.parent_ligand_chain,
        "parent_chain_map": args.parent_chain_map, "weights": args.weights,
        "params": str(Path(args.params).resolve()) if args.params else "",
        "template": str(Path(args.ligand_template_pdb).resolve()) if args.ligand_template_pdb else "",
        "scheme": args.scheme, "cdr_definition": args.cdr_definition,
        "rosetta_ligand_chain": getattr(args, "rosetta_ligand_chain", "X"),
        "ligand_resname": getattr(args, "ligand_resname", "E2L"),
        "ligand_smiles": getattr(args, "ligand_smiles", None),
        "ligand_source_sha256": sha256(Path(args.ligand_source).read_bytes()).hexdigest() if getattr(args, 'ligand_source', None) else '',
    }
    if any(value == "not_run" for value in analysis_status(args).values()):
        payload["analysis_status"] = analysis_status(args)
    if getattr(args, "legacy_pre_talaris", None) is not None or Path(args.weights).is_file() or args.weights != "ligand":
        from .pyrosetta_score import scoring_weights_identity
        from importlib.util import find_spec
        package = find_spec("pyrosetta")
        database = Path(package.origin).parent / "database" if package and package.origin else None
        payload["scoring_weights_identity"] = scoring_weights_identity(args.weights, getattr(args, "legacy_pre_talaris", None), database=database)
    for key, value in (("params_sha256", args.params), ("template_sha256", args.ligand_template_pdb)):
        payload[key] = sha256(Path(value).read_bytes()).hexdigest() if value else ""
    payload["scoring_code_sha256"] = sha256(b"".join(
        (Path(__file__).parent / name).read_bytes() for name in ("pyrosetta_score.py", "hbonds.py", "validation.py", "arpeggio.py")
    )).hexdigest()
    return sha1(json.dumps(payload, sort_keys=True).encode("utf-8")).hexdigest()[:12]


def _ensure_parent_cache(
    args: argparse.Namespace, parent: Path, output: Path, batch_log
) -> Path:
    requested = getattr(args, "reuse_parent_cache", None)
    if requested:
        from .paths import native_path
        from .fast_workflow import _cache_paths
        from .architecture_batch import atom_snapshot

        shared = Path(requested).expanduser().resolve()
        if not _completed(shared):
            raise ValueError(f"Requested shared parent cache is incomplete: {shared}")
        if not args.params or not args.ligand_template_pdb:
            raise ValueError("Shared parent reuse requires explicit matching ligand params and template")
        manifest_path = shared / "workflow_manifest.json"
        manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
        old_parent = native_path(manifest["parent"]).resolve()
        if not old_parent.is_file():
            raise FileNotFoundError(f"Shared cache's original parent is missing: {old_parent}")
        parent_hash = sha256(parent.read_bytes()).hexdigest()
        if sha256(old_parent.read_bytes()).hexdigest() != parent_hash:
            raise ValueError("Requested shared cache has different parent coordinates or atom content")
        expected_name = "_parent_cache_" + _parent_cache_key(args, old_parent)
        if shared.name != expected_name:
            raise ValueError(
                f"Requested shared parent cache settings/code differ: expected {expected_name}, found {shared.name}"
            )
        validate_scoring_reuse(shared, Path(args.params).resolve(), args.weights, getattr(args, "legacy_pre_talaris", None))
        validate_analysis_reuse(shared, args)
        artifacts = _cache_paths(shared)
        missing = [str(path) for path in artifacts.values() if not path.is_file()]
        if missing:
            raise FileNotFoundError("Requested shared parent artifacts are missing: " + "; ".join(missing))
        if atom_snapshot(artifacts["prepared"]) != atom_snapshot(parent):
            raise ValueError("Shared prepared parent differs from requested parent named atoms/coordinates")
        output.mkdir(parents=True, exist_ok=True)
        provenance = {
            "status": "verified", "source_cache": str(shared),
            "original_parent": str(old_parent), "requested_parent": str(parent),
            "parent_sha256": parent_hash, "settings_and_code_cache_key": expected_name,
            "workflow_manifest_sha256": sha256(manifest_path.read_bytes()).hexdigest(),
            "artifact_sha256": {name: sha256(path.read_bytes()).hexdigest() for name, path in artifacts.items()},
            "parent_scoring_repeated": False,
        }
        (output / "shared_parent_cache_reuse.json").write_text(
            json.dumps(provenance, indent=2) + "\n", encoding="utf-8"
        )
        batch_log.write(f"[{datetime.now(timezone.utc).isoformat()}] REUSE VERIFIED SHARED PARENT CACHE {shared}\n")
        batch_log.flush()
        return shared
    cache = output / f"_parent_cache_{_parent_cache_key(args, parent)}"
    if _completed(cache):
        validate_analysis_reuse(cache, args)
        if args.params:
            validate_scoring_reuse(cache, Path(args.params).resolve(), args.weights, getattr(args, "legacy_pre_talaris", None))
        batch_log.write(f"[{datetime.now(timezone.utc).isoformat()}] REUSE PARENT CACHE {cache}\n")
        batch_log.flush()
        return cache
    command = [
        sys.executable, "-m", "fv_energy_analysis.workflow",
        "--parent", str(parent), "--mutant", str(parent), "--output", str(cache),
        "--fv-chains", args.fv_chains, "--ligand-chain", args.ligand_chain,
        "--rosetta-ligand-chain", args.rosetta_ligand_chain,
        "--ligand-resname", args.ligand_resname, "--scheme", args.scheme,
        "--cdr-definition", args.cdr_definition, "--fit-region", args.fit_region,
        "--rmsd-atoms", args.rmsd_atoms,
        "--rotamer-rmsd-cutoff", str(args.rotamer_rmsd_cutoff),
        "--hcdr3-tip-sequence", getattr(args, "hcdr3_tip_sequence", "YYGSS"),
        "--network-energy-presence-cutoff", str(args.network_energy_presence_cutoff),
        "--network-energy-change-cutoff", str(args.network_energy_change_cutoff),
        "--weights", args.weights,
        "--arpeggio-scope", "full",
    ]
    command.extend(analysis_flags(args))
    if getattr(args, "legacy_pre_talaris", None) is not None:
        command.append("--legacy-pre-talaris" if args.legacy_pre_talaris else "--no-legacy-pre-talaris")
    optional = {
        "--parent-ligand-chain": args.parent_ligand_chain,
        "--mutant-ligand-chain": args.parent_ligand_chain,
        "--parent-chain-map": _map_arg(args.parent_chain_map),
        "--mutant-chain-map": _map_arg(args.parent_chain_map),
        "--params": args.params,
        "--ligand-template-pdb": args.ligand_template_pdb,
        "--molfile-to-params": args.molfile_to_params,
        "--ligand-smiles": getattr(args, "ligand_smiles", None),
        "--ligand-source": getattr(args, "ligand_source", None),
    }
    for flag, value in optional.items():
        if value:
            command.extend([
                flag,
                str(Path(value).expanduser().resolve())
                if flag in {"--params", "--ligand-template-pdb", "--molfile-to-params", "--ligand-source"}
                else str(value),
            ])
    batch_log.write(f"[{datetime.now(timezone.utc).isoformat()}] BUILD PARENT CACHE {cache}\n")
    batch_log.flush()
    completed = subprocess.run(command, stdout=batch_log, stderr=subprocess.STDOUT, text=True, check=False)
    if completed.returncode or not _completed(cache):
        raise RuntimeError(f"Parent cache generation failed; see {cache / 'workflow.log'}")
    return cache


def run_batch(args: argparse.Namespace) -> Path:
    started = time.time()
    parent = Path(args.parent).expanduser().resolve()
    input_folder = Path(args.input_folder).expanduser().resolve()
    output = Path(args.output).expanduser().resolve()
    output.mkdir(parents=True, exist_ok=True)
    if (output / 'STOPPED_FOR_VALIDATION.md').exists():
        raise ValueError('This legacy batch was stopped for parameter validation. Preserve it and choose a new output directory for corrected scoring.')
    if not parent.is_file():
        raise FileNotFoundError(f"Parent structure not found: {parent}")
    if not input_folder.is_dir():
        raise NotADirectoryError(f"Input folder not found: {input_folder}")
    if bool(args.params) != bool(args.ligand_template_pdb):
        raise ValueError("Supply both --params and --ligand-template-pdb, or neither")
    if args.jobs < 1:
        raise ValueError("--jobs must be at least 1")
    if not math.isfinite(args.arpeggio_neighbourhood_radius) or args.arpeggio_neighbourhood_radius <= 0:
        raise ValueError("Arpeggio neighbourhood radius must be finite and positive")
    reuse_batch = (
        Path(args.reuse_mutants_from).expanduser().resolve()
        if args.reuse_mutants_from else None
    )
    if reuse_batch is not None and not reuse_batch.is_dir():
        raise NotADirectoryError(f"Saved mutant batch not found: {reuse_batch}")

    candidates = [path for path in _discover(input_folder, args.glob, args.recursive) if path != parent]
    if not candidates:
        raise FileNotFoundError("No candidate complexes matched the requested input globs")
    names = _result_names(candidates, input_folder)
    manifest_path = output / "batch_manifest.json"
    log_path = output / "batch.log"
    manifest: dict[str, object] = {
        "status": "running", "started_utc": datetime.now(timezone.utc).isoformat(),
        "parent": str(parent), "input_folder": str(input_folder), "output": str(output),
        "candidate_count": len(candidates), "candidates": [],
        "analysis_status": analysis_status(args),
    }
    manifest_path.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    rows: list[dict[str, object]] = []
    mutation_rows: list[dict[str, object]] = []
    network_rows: list[dict[str, object]] = []
    bond_rows: list[dict[str, object]] = []

    results_by_candidate: dict[Path, dict[str, object]] = {}
    with log_path.open("a", encoding="utf-8") as batch_log:
        parent_cache = _ensure_parent_cache(args, parent, output, batch_log)
        if args.params:
            batch_params = Path(args.params).resolve()
            batch_template = Path(args.ligand_template_pdb).resolve()
        else:
            batch_params = parent_cache / "prepared" / "ligand_definition" / f"{args.ligand_resname}.params"
            batch_template = parent_cache / "prepared" / "ligand_definition" / f"{args.ligand_resname}_0001.pdb"
        if not batch_params.is_file() or not batch_template.is_file():
            raise FileNotFoundError(
                f"Parent cache is missing ligand params/template: {batch_params}; {batch_template}"
            )
        base_config: dict[str, object] = {
            "parent": str(parent), "parent_cache": str(parent_cache),
            "fv_chains": args.fv_chains, "ligand_chain": args.ligand_chain,
            "mutant_ligand_chain": args.mutant_ligand_chain,
            "mutant_chain_map": args.mutant_chain_map,
            "rosetta_ligand_chain": args.rosetta_ligand_chain,
            "ligand_resname": args.ligand_resname, "scheme": args.scheme,
            "cdr_definition": args.cdr_definition, "fit_region": args.fit_region,
            "rmsd_atoms": args.rmsd_atoms,
            "rotamer_rmsd_cutoff": args.rotamer_rmsd_cutoff,
            "hcdr3_tip_sequence": getattr(args, "hcdr3_tip_sequence", "YYGSS"),
            "network_energy_presence_cutoff": args.network_energy_presence_cutoff,
            "network_energy_change_cutoff": args.network_energy_change_cutoff,
            "weights": args.weights, "params": str(batch_params),
            "arpeggio_scope": args.arpeggio_scope,
            "arpeggio_neighbourhood_radius": args.arpeggio_neighbourhood_radius,
            "ligand_template_pdb": str(batch_template),
            "legacy_pre_talaris": getattr(args, "legacy_pre_talaris", None),
            **{name: getattr(args, name, False) for name in ("skip_arpeggio", "skip_plip", "skip_network_analysis", "skip_rotamer_analysis", "skip_e2_architecture")},
        }
        from .ligand_chemistry import reference_smiles_from_source
        base_config["ligand_smiles"] = getattr(args, "ligand_smiles", None) or reference_smiles_from_source(
            getattr(args, 'ligand_source', None), batch_params, batch_template
        )
        pending: list[tuple[Path, dict[str, object]]] = []
        for candidate in candidates:
            result_dir = output / names[candidate]
            if _completed(result_dir) and not args.force:
                validate_scoring_reuse(result_dir, batch_params, args.weights, getattr(args, "legacy_pre_talaris", None))
                validate_analysis_reuse(result_dir, args)
                validate_comparison_reuse(result_dir, args, candidate)
                saved_selection = json.loads((result_dir / "workflow_manifest.json").read_text()).get("arpeggio_selection", {})
                if not getattr(args, "skip_arpeggio", False) and saved_selection.get("scope", "full") != "full" and (
                    args.arpeggio_scope == "full" or
                    args.arpeggio_neighbourhood_radius > saved_selection.get("radius_A", 0)
                ):
                    raise ValueError("Saved Arpeggio query is narrower than requested; regenerate contacts without rescoring Rosetta")
                results_by_candidate[candidate] = {
                    "candidate": str(candidate), "result_folder": str(result_dir),
                    "status": "complete", "error": "", "reused": True,
                }
            else:
                candidate_config = {
                    **base_config, "candidate": str(candidate), "output": str(result_dir),
                }
                if reuse_batch is not None:
                    source_result = reuse_batch / names[candidate]
                    if _completed(source_result):
                        validate_scoring_reuse(source_result, batch_params, args.weights, getattr(args, "legacy_pre_talaris", None))
                        validate_analysis_reuse(source_result, args)
                        # A basename alone is not sufficient evidence for reusing
                        # scientific results from a differently located input.
                        source_manifest = json.loads((source_result / "workflow_manifest.json").read_text(encoding="utf-8"))
                        source_input = Path(source_manifest["mutant"])
                        if not source_input.is_file() or sha256(source_input.read_bytes()).digest() != sha256(candidate.read_bytes()).digest():
                            raise ValueError(f"Saved mutant input does not match {candidate.name}; refusing reuse")
                        candidate_config["source_result"] = str(source_result)
                pending.append((candidate, candidate_config))
        manifest["optimization"] = {
            "parent_cache": str(parent_cache), "parent_scored_once": True,
            "persistent_pyrosetta_workers": args.jobs,
            "parallel_jobs": args.jobs, "anarcii_sequence_cache_per_worker": True,
            "arpeggio_scope": args.arpeggio_scope,
            "arpeggio_neighbourhood_radius_A": args.arpeggio_neighbourhood_radius,
            "mutant_artifacts_reused_from": str(reuse_batch) if reuse_batch else None,
        }
        batch_log.write(
            f"[{datetime.now(timezone.utc).isoformat()}] FAST BATCH jobs={args.jobs} "
            f"pending={len(pending)} reused={len(candidates) - len(pending)}\n"
        )
        batch_log.flush()
        if pending:
            context = multiprocessing.get_context("spawn")
            with concurrent.futures.ProcessPoolExecutor(
                max_workers=args.jobs, mp_context=context
            ) as executor:
                future_map = {
                    executor.submit(run_rebased_candidate if config.get("source_result") else run_cached_candidate, config): candidate
                    for candidate, config in pending
                }
                for future in concurrent.futures.as_completed(future_map):
                    candidate = future_map[future]
                    try:
                        result = future.result()
                    except Exception as exc:
                        result = {
                            "candidate": str(candidate),
                            "result_folder": str(output / names[candidate]),
                            "status": "failed", "error": f"Worker crashed: {exc}",
                        }
                    result["reused"] = False
                    results_by_candidate[candidate] = result
                    if result["status"] == "failed":
                        record_issue(output, result)
                    batch_log.write(
                        f"[{datetime.now(timezone.utc).isoformat()}] "
                        f"{str(result['status']).upper()} {candidate.name}"
                        + (f" :: {result['error']}" if result.get("error") else "") + "\n"
                    )
                    batch_log.flush()
                    manifest["candidates"] = [
                        results_by_candidate[path] for path in candidates
                        if path in results_by_candidate
                    ]
                    manifest_path.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
                    if result["status"] == "failed" and args.stop_on_error:
                        for queued in future_map:
                            queued.cancel()
                        break

    for candidate in candidates:
        if candidate not in results_by_candidate:
            continue
        result = results_by_candidate[candidate]
        result_dir = Path(str(result["result_folder"]))
        status, error = str(result["status"]), str(result.get("error", ""))
        try:
            row, interpretations, interaction_types, bond_changes = _summary_row(
                candidate, result_dir, status, error
            )
        except Exception as exc:
            status, error = "failed", f"Could not aggregate result: {exc}"
            record_issue(output, {**result, "error": error})
            row, interpretations, interaction_types, bond_changes = _summary_row(
                candidate, result_dir, status, error
            )
        row["reused_existing_result"] = bool(result.get("reused"))
        rows.append(row)
        mutation_rows.extend(interpretations)
        network_rows.extend(interaction_types)
        bond_rows.extend(bond_changes)

    summary_csv, mutation_csv, network_csv, bond_csv, e2_screen_csv, summary_markdown, column_definitions = _write_master_reports(
        output, rows, mutation_rows, network_rows, bond_rows
    )
    failure_count = sum(row["status"] != "complete" for row in rows)
    manifest.update({
        "status": "complete_with_failures" if failure_count else "complete",
        "completed_utc": datetime.now(timezone.utc).isoformat(),
        "elapsed_seconds": round(time.time() - started, 3),
        "completed_count": len(rows) - failure_count, "failure_count": failure_count,
        "issues_csv": str(output / "batch_issues.csv") if (output / "batch_issues.csv").exists() else None,
        "summary_csv": str(summary_csv),
        "mutation_summary_csv": str(mutation_csv),
        "interaction_network_summary_csv": str(network_csv),
        "hbond_saltbridge_changes_csv": str(bond_csv),
        "e2_architecture_screen_csv": str(e2_screen_csv),
        "summary_markdown": str(summary_markdown),
        "column_definitions": str(column_definitions),
    })
    manifest_path.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    return output


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    try:
        output = run_batch(args)
    except Exception as exc:
        path = Path(args.output) / 'batch_manifest.json'
        if path.is_file():
            manifest = json.loads(path.read_text())
            if manifest.get('status') == 'running' and not (path.parent/'STOPPED_FOR_VALIDATION.md').exists():
                manifest.update(status='failed', error=str(exc), completed_utc=datetime.now(timezone.utc).isoformat())
                path.write_text(json.dumps(manifest,indent=2)+'\n')
        raise
    print(f"Batch analysis complete: {output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
