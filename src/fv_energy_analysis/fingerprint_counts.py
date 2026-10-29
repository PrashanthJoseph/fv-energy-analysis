"""Read and normalize the four tools' saved interaction counts."""

from __future__ import annotations

import csv
import json
from collections import Counter
from pathlib import Path
from typing import Any
from .paths import native_path

CountValue = int | str | None  # integer, unsupported "NA", or missing data

COUNT_COLUMNS = (
    "total", "hbond_total", "conventional_hbond", "carbon_hbond", "weak_hbond",
    "water_bridge", "vdw", "vdw_clash", "hydrophobic", "salt_bridge",
    "attractive_charge", "pi_stack", "pi_t_stack", "pi_cation", "pi_anion",
    "pi_alkyl", "pi_sigma", "donor_pi", "lone_pair_pi", "pi_sulfur", "alkyl",
    "halogen_bond", "metal", "unfavorable", "aromatic_contact", "carbon_pi",
    # Native Arpeggio labels and its native plane-plane record type.
    "hbond", "ionic", "aromatic", "CATIONPI", "CARBONPI", "DONORPI", "METSULPHURPI", "plane-plane",
)


def _read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8-sig") as handle:
        return list(csv.DictReader(handle))


def _windows_path(value: str | Path) -> Path:
    return native_path(value)


def _empty() -> dict[str, CountValue]:
    return {name: "NA" for name in COUNT_COLUMNS}


def analysis_status(result_dir: Path, summary: dict[str, Any] | None = None) -> dict[str, str]:
    """Read explicit stage selection; old completed runs retain their defaults."""
    status = dict((summary or {}).get("analysis_status") or {})
    path = result_dir / "workflow_manifest.json"
    if path.is_file():
        status.update(json.loads(path.read_text(encoding="utf-8")).get("analysis_status") or {})
    status.update((summary or {}).get("report_analysis_status") or {})
    return status


def _not_run(status, name):
    return status.get(name) in {"not_run", "skipped", "unavailable"}


def _unavailable_counts():
    return {name: None for name in COUNT_COLUMNS}


def _delta(parent: CountValue, mutant: CountValue) -> CountValue:
    if parent is None or mutant is None:
        return None
    if isinstance(parent, int) and isinstance(mutant, int):
        return mutant - parent
    return "NA"


def _row(candidate: str, result_folder: Path, tool: str, unit: str,
         parent: dict[str, CountValue], mutant: dict[str, CountValue],
         energies: dict[str, Any] | None = None) -> dict[str, object]:
    output: dict[str, object] = {
        "candidate": candidate,
        "result_folder": str(result_folder),
        "tool": tool,
        "count_unit": unit,
    }
    for name in COUNT_COLUMNS:
        output[f"parent_{name}"] = parent[name]
        output[f"mutant_{name}"] = mutant[name]
        output[f"delta_{name}"] = _delta(parent[name], mutant[name])
    output.update(energies or {})
    return output


def _rosetta_counts(summary: dict[str, Any], label: str) -> dict[str, CountValue]:
    values = _empty()
    hbonds = (summary.get("ligand_hbonds") or {}).get(label)
    values["total"] = len(hbonds) if hbonds is not None else None
    values["hbond_total"] = values["total"]
    # Rosetta's scored H-bond records are not classified as conventional vs C-H.
    return values


def _rosetta_energies(summary: dict[str, Any]) -> dict[str, object]:
    terms = {row.get("metric"): row for row in summary.get("interface_terms", [])}
    output: dict[str, object] = {}
    chain = summary.get("ligand_chain", "X")
    for label, metric in (
        ("fa_atr_REU", "if_X_fa_atr"),
        ("fa_rep_REU", "if_X_fa_rep"),
        ("fa_elec_REU", "if_X_fa_elec"),
        ("hbond_sc_REU", "if_X_hbond_sc"),
        ("fa_pair_REU", "if_X_fa_pair"),
        ("hbond_bb_sc_REU", "if_X_hbond_bb_sc"),
        ("solvation_REU", "if_X_fa_sol"),
    ):
        row = terms.get(metric.replace("if_X_", f"if_{chain}_"), {})
        output[f"parent_{label}"] = row.get("parent")
        output[f"mutant_{label}"] = row.get("mutant")
        output[f"delta_{label}"] = row.get("delta_mutant_minus_parent")
    return output


def _ligand_resname(summary: dict[str, Any], result_dir: Path) -> str:
    """Read recorded ligand identity; never guess a chemical identity from a name."""
    if summary.get("ligand_resname"):
        return str(summary["ligand_resname"])
    for data in (summary.get("plip_ligand_interactions") or {}).values():
        if data and data.get("ligand_resname"):
            return str(data["ligand_resname"])
    manifest_path = result_dir / "workflow_manifest.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8")) if manifest_path.is_file() else {}
    if manifest.get("ligand_resname"):
        return str(manifest["ligand_resname"])
    if manifest.get("params"):
        params = native_path(manifest["params"])
        if not params.is_file():
            params = result_dir / "prepared" / "ligand_definition" / params.name
        if params.is_file():
            lines = [line.split() for line in params.read_text().splitlines()]
            for field in ("IO_STRING", "NAME"):
                name = next((parts[1] for parts in lines if len(parts) > 1 and parts[0] == field), None)
                if name:
                    return name
    raise ValueError(f"No recorded ligand residue name in {result_dir}; supply saved ligand metadata or params")


def _is_ligand_protein(row: dict[str, Any], ligand: str) -> bool:
    endpoints = (row.get("bgn", {}), row.get("end", {}))
    has_ligand = any(str(x.get("label_comp_id", "")).upper() == ligand.upper() for x in endpoints)
    has_protein = any(str(x.get("label_comp_type", "")).upper() == "P" for x in endpoints)
    return has_ligand and has_protein


def _arpeggio_counts(path: Path, ligand_resname: str) -> dict[str, CountValue]:
    values = _empty()
    labels: Counter[str] = Counter()
    ring_contacts = 0
    for row in json.loads(path.read_text(encoding="utf-8")):
        if _is_ligand_protein(row, ligand_resname):
            labels.update(str(x) for x in row.get("contact", []))
            if row.get("type") == "plane-plane":
                ring_contacts += 1
    for name in ("hbond", "weak_hbond", "vdw", "vdw_clash", "hydrophobic", "ionic",
                 "aromatic", "CATIONPI", "CARBONPI", "DONORPI", "METSULPHURPI"):
        values[name] = labels[name]
    values["plane-plane"] = ring_contacts
    # Arpeggio reports weak H-bonds but does not identify all of them as C-H bonds.
    values["total"] = sum(labels.values())
    return values

TOOL_PREFIXES = {"Rosetta": "rosetta", "Arpeggio": "arpeggio", "PLIP": "plip", "Discovery Studio": "ds"}

def _arpeggio_raw(path: Path, ligand_resname: str) -> Counter[str]:
    labels: Counter[str] = Counter()
    for row in json.loads(path.read_text(encoding="utf-8")):
        if _is_ligand_protein(row, ligand_resname):
            labels.update(str(x) for x in row.get("contact", []))
    return labels


def _plip_counts(summary: dict[str, Any], label: str) -> dict[str, CountValue]:
    values = _empty()
    data = (summary.get("plip_ligand_interactions") or {}).get(label) or {}
    counts = data.get("interaction_counts")
    if data.get("binding_site_count") == 0:
        counts = None  # No matching ligand is not a verified zero-interaction site.
    mapping = {
        "hbond_total": "hydrogen_bonds", "water_bridge": "water_bridges",
        "hydrophobic": "hydrophobic_interactions", "salt_bridge": "salt_bridges",
        "pi_stack": "pi_stacks", "pi_cation": "pi_cation_interactions",
        "halogen_bond": "halogen_bonds", "metal": "metal_complexes",
    }
    for output, source in mapping.items():
        value = counts.get(source) if counts is not None else None
        values[output] = int(value) if value is not None else None
    values["total"] = (int(data["total_interaction_count"])
                       if counts is not None and data.get("total_interaction_count") is not None else None)
    return values


def _ds_counts(path: Path) -> dict[str, CountValue]:
    if not path.is_file():
        return _unavailable_counts()
    values = _empty()
    rows = _read_csv(path)
    counts = Counter(row["interaction_type"] for row in rows)
    direct = {
        "conventional_hbond": "conventional_hbond", "carbon_hbond": "carbon_hbond",
        "water_bridge": "water_bridge_hbond", "vdw_clash": "unfavorable_vdw_bump",
        "salt_bridge": "salt_bridge", "attractive_charge": "attractive_charge",
        "pi_stack": "pi_pi_stacked", "pi_t_stack": "pi_pi_t_shaped",
        "pi_cation": "pi_cation", "pi_anion": "pi_anion", "pi_alkyl": "pi_alkyl",
        "pi_sigma": "pi_sigma", "donor_pi": "donor_pi", "lone_pair_pi": "lone_pair_pi",
        "pi_sulfur": "pi_sulfur", "alkyl": "alkyl", "halogen_bond": "halogen_cl_br_i",
        "metal": "metal_acceptor",
    }
    for output, source in direct.items():
        values[output] = counts[source]
    # A Trp can contribute two fused-ring records to the same ligand pair.
    # Summary pi-stacks count unique protein residue partners, not rings.
    stack_partners = set()
    for row in rows:
        if row["interaction_type"] != "pi_pi_stacked":
            continue
        partner = row.get("protein_residue") or row.get("residue_pair_id")
        if not partner:
            raise ValueError(f"DS pi-stack record lacks residue identity: {path}")
        stack_partners.add(partner)
    values["pi_stack"] = len(stack_partners)
    values["hbond_total"] = sum(counts[name] for name in (
        "conventional_hbond", "carbon_hbond", "water_bridge_hbond", "water_hbond"
    ))
    values["hydrophobic"] = sum(counts[name] for name in ("alkyl", "pi_alkyl"))
    values["vdw"] = "NA"  # DS monitor does not emit favorable generic VdW as a separate row.
    values["unfavorable"] = sum(value for name, value in counts.items() if name.startswith("unfavorable_"))
    values["total"] = sum(counts.values())
    return values


def _raw_count_rows(candidate: str, result_dir: Path, summary: dict[str, Any],
                    parent_ds_path: Path) -> list[dict[str, object]]:
    from .plip import INTERACTION_TAGS
    ligand_resname = _ligand_resname(summary, result_dir)
    status = analysis_status(result_dir, summary)
    sources = []
    sources.append((
        "Rosetta", "scored ligand H-bond records",
        {"scored_hbond": _rosetta_counts(summary, "parent")["hbond_total"]},
        {"scored_hbond": _rosetta_counts(summary, "mutant")["hbond_total"]},
    ))
    sources.append((
        "Arpeggio", "atom/group/plane contact annotations",
        {} if _not_run(status, "arpeggio") else _arpeggio_raw(result_dir / "arpeggio_parent" / "parent_arpeggio.json", ligand_resname),
        {} if _not_run(status, "arpeggio") else _arpeggio_raw(result_dir / "arpeggio_mutant" / "mutant_arpeggio.json", ligand_resname),
    ))
    plip = summary.get("plip_ligand_interactions") or {}
    def plip_native(label):
        if _not_run(status, "plip"):
            return {name: None for name in INTERACTION_TAGS}
        data = plip.get(label) or {}
        counts = data.get("interaction_counts")
        if counts is None or data.get("binding_site_count") == 0:
            return {name: None for name in INTERACTION_TAGS}
        return {str(k): int(v) if v is not None else None for k, v in counts.items()}
    sources.append((
        "PLIP", "reported interaction entries",
        plip_native("parent"), plip_native("mutant"),
    ))
    sources.append((
        "Discovery Studio", "nonbond-monitor interaction annotations",
        Counter(row["interaction_type"] for row in _read_csv(parent_ds_path))
        if not _not_run(status, "discovery_studio") and parent_ds_path.is_file() else {},
        Counter(row["interaction_type"] for row in _read_csv(result_dir / "discovery_studio_interactions.csv"))
        if not _not_run(status, "discovery_studio") and (result_dir / "discovery_studio_interactions.csv").is_file() else {},
    ))
    rows: list[dict[str, object]] = []
    for tool, unit, parent, mutant in sources:
        key = {"Arpeggio": "arpeggio", "PLIP": "plip", "Discovery Studio": "discovery_studio"}.get(tool)
        if key and _not_run(status, key):
            unit = "not_run (disabled)"
        for interaction_type in sorted(set(parent) | set(mutant)):
            # An absent label in an available detector file is zero. Missing
            # summary entries remain unavailable instead of inventing zeros.
            p = parent.get(interaction_type, 0) if isinstance(parent, Counter) else parent.get(interaction_type)
            m = mutant.get(interaction_type, 0) if isinstance(mutant, Counter) else mutant.get(interaction_type)
            rows.append({
                "candidate": candidate, "result_folder": str(result_dir), "tool": tool,
                "count_unit": unit, "detector_interaction_type": interaction_type,
                "parent_count": p, "mutant_count": m,
                "delta_count": _delta(p, m),
            })
    return rows

def candidate_fingerprints(candidate: str, result_dir: Path,
                           parent_ds_path: Path, summary: dict[str, Any] | None = None) -> list[dict[str, object]]:
    if summary is None:
        summary = json.loads((result_dir / "automatic_summary.json").read_text(encoding="utf-8"))
    status = analysis_status(result_dir, summary)
    ligand_resname = _ligand_resname(summary, result_dir)
    rows = [_row(
        candidate, result_dir, "Rosetta", "scored ligand H-bond records",
        _rosetta_counts(summary, "parent"), _rosetta_counts(summary, "mutant"),
        _rosetta_energies(summary),
    )]
    rows.append(_row(
        candidate, result_dir, "Arpeggio", "atom/group/plane contact annotations",
        _unavailable_counts() if _not_run(status, "arpeggio") else _arpeggio_counts(result_dir / "arpeggio_parent" / "parent_arpeggio.json", ligand_resname),
        _unavailable_counts() if _not_run(status, "arpeggio") else _arpeggio_counts(result_dir / "arpeggio_mutant" / "mutant_arpeggio.json", ligand_resname),
    ))
    rows.append(_row(
        candidate, result_dir, "PLIP", "reported interaction entries",
        _unavailable_counts() if _not_run(status, "plip") else _plip_counts(summary, "parent"),
        _unavailable_counts() if _not_run(status, "plip") else _plip_counts(summary, "mutant"),
    ))
    rows.append(_row(
        candidate, result_dir, "Discovery Studio", "nonbond-monitor interaction annotations",
        _unavailable_counts() if _not_run(status, "discovery_studio") else _ds_counts(parent_ds_path),
        _unavailable_counts() if _not_run(status, "discovery_studio") else _ds_counts(result_dir / "discovery_studio_interactions.csv"),
    ))
    for row in rows:
        key = {"Arpeggio": "arpeggio", "PLIP": "plip", "Discovery Studio": "discovery_studio"}.get(row["tool"])
        if key and _not_run(status, key):
            row["count_unit"] = "not_run (disabled)"
    return rows
