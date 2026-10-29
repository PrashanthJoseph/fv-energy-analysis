"""Write per-mutant reports and the canonical grouped batch summary."""

from __future__ import annotations

import json
from hashlib import sha256
import math
import re
import shutil
from pathlib import Path

from .fingerprint_counts import COUNT_COLUMNS, TOOL_PREFIXES, _read_csv, _windows_path, _raw_count_rows, candidate_fingerprints, _ligand_resname, analysis_status
from .report import write_csv, write_json
from .resources import resource_path

# Essential, separately addressable columns in the four fingerprint sections.
# Rosetta geometric classes that it does not report are deliberately absent.
ESSENTIAL_COUNTS = {
    "Rosetta": (("Scored H-bonds", "hbond_total"),),
    "Arpeggio": (
        ("hbond", "hbond"), ("weak_hbond", "weak_hbond"),
        ("vdw", "vdw"), ("vdw_clash", "vdw_clash"),
        ("hydrophobic", "hydrophobic"), ("ionic", "ionic"),
        ("aromatic", "aromatic"), ("plane-plane", "plane-plane"),
        ("CATIONPI", "CATIONPI"), ("CARBONPI", "CARBONPI"), ("DONORPI", "DONORPI"),
    ),
    "PLIP": (
        ("H-bonds", "hbond_total"), ("Hydrophobic", "hydrophobic"),
        ("Salt bridges", "salt_bridge"), ("Pi-stacks", "pi_stack"),
        ("Pi-cation", "pi_cation"), ("Water bridges", "water_bridge"),
    ),
    "Discovery Studio": (
        ("Conventional H-bonds", "conventional_hbond"), ("Carbon H-bonds", "carbon_hbond"),
        ("Pi-stacks (unique residues)", "pi_stack"), ("Pi T-shaped", "pi_t_stack"),
        ("Pi-cation", "pi_cation"), ("Pi-anion", "pi_anion"),
        ("Pi-alkyl", "pi_alkyl"), ("Alkyl", "alkyl"), ("Pi-sigma", "pi_sigma"),
        ("Attractive charge", "attractive_charge"), ("Salt bridges", "salt_bridge"),
        ("Unfavorable", "unfavorable"),
    ),
}

ARCHITECTURE_COLUMNS = (
    ("e2_architecture_status", "Architecture status", "text"),
    ("e2_o3_hbond", "DS O3 H-bond", "boolean"),
    ("e2_o17_hbond", "DS O17 H-bond", "boolean"),
    ("e2_c18_hydrophobic_contact", "DS C18 packing", "boolean"),
    ("e2_c18_hydrophobic_partner_count", "DS C18 packing residues", "integer"),
    ("e2_geometric_sandwich", "Steric sandwich (geometry)", "boolean"),
    ("e2_steric_alpha_residue_count", "Alpha-face contact residues", "integer"),
    ("e2_steric_beta_residue_count", "Beta-face contact residues", "integer"),
    ("e2_steric_opposing_residue_pair_count", "Opposing residue pairs", "integer"),
    ("e2_steric_overlapping_atom_pair_count", "A-ring atom-pair overlaps", "integer"),
    ("e2_a_ring_pi_contact_count", "DS A-ring pi-mediated residues", "integer"),
    ("e2_a_ring_pi_contact_residues", "DS A-ring pi-mediated contacts", "text"),
    ("e2_conserved_architecture_score_0_to_4", "Recognition features / 4", "integer"),
    ("e2_strict_screen_pass", "All four recognition features", "boolean"),
)


def _saved_architecture(result_dir: Path, role: str, ds_path: Path) -> dict[str, object]:
    """Recompute only cheap geometry from recorded inputs, never run detectors."""
    from .estrogen_architecture import analyze_estrogen_architecture
    from .atom_roles import load_atom_roles
    summary_path = result_dir / "final" / "summary.json"
    manifest_path = result_dir / "workflow_manifest.json"
    if not summary_path.is_file() or not manifest_path.is_file():
        return {"e2_architecture_status": "unavailable: recorded inputs missing"}
    summary = json.loads(summary_path.read_text(encoding="utf-8"))
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    if manifest.get("analysis_status", {}).get("e2_architecture") == "not_run":
        return {"e2_architecture_status": "not_run"}
    if not manifest.get("params"):
        return {"e2_architecture_status": "unavailable: ligand params missing"}
    params = _windows_path(manifest["params"])
    pdb = _windows_path(summary[role])
    if not pdb.is_file():
        pdb = result_dir / "prepared" / pdb.name
    if not params.is_file():
        params = result_dir / "prepared" / "ligand_definition" / params.name
    if not params.is_file() or not pdb.is_file():
        return {"e2_architecture_status": "unavailable: recorded params or PDB missing"}
    template = (_windows_path(manifest["ligand_template_pdb"]) if manifest.get("ligand_template_pdb")
                else params.with_name(params.stem + "_0001.pdb"))
    if not template.is_file():
        template = params.parent / template.name
    if not template.is_file():
        return {"e2_architecture_status": "unavailable: ligand template missing"}
    roles = load_atom_roles(params, template)
    recorded_hash = manifest.get("scoring_protocol", {}).get("params_sha256")
    if recorded_hash and recorded_hash != sha256(params.read_bytes()).hexdigest():
        raise ValueError("Saved scores and current ligand parameters differ; refusing mixed reporting")
    result = analyze_estrogen_architecture(
        pdb,
        ligand_resname=_ligand_resname(summary, result_dir),
        fv_chains=set(summary["fv_chains"]), **roles,
        discovery_studio_csv=ds_path,
    )
    result["e2_atom_role_source"] = "validated parameter/template hash-bound cache"
    result.setdefault("e2_architecture_status", "complete" if ds_path.is_file() else "unavailable: Discovery Studio output missing")
    return result


def _typed(value: object, kind: str = "number") -> object:
    if value in (None, "", "NA", "NaN", "nan"):
        return None
    if kind == "text":
        return str(value)
    if kind == "boolean":
        if str(value).strip().lower() not in ("true", "false"):
            raise ValueError(f"Invalid boolean value: {value!r}")
        return str(value).strip().lower() == "true"
    if isinstance(value, bool):
        raise ValueError(f"Invalid numeric value: {value!r}")
    number = float(value)
    if math.isnan(number):
        return None
    if not math.isfinite(number) or kind == "integer" and not number.is_integer():
        raise ValueError(f"Invalid {kind} value: {value!r}")
    return int(number) if kind == "integer" else number


def _grouped_columns(weights: str = "ref2015") -> list[dict[str, object]]:
    columns: list[dict[str, object]] = [
        {"key": "rank", "label": "Rank", "section": "Identity", "kind": "integer", "width": 7},
        {"key": "mutations", "label": "Mutations", "section": "Identity", "kind": "text", "width": 29},
        {"key": "candidate", "label": "Complex", "section": "Identity", "kind": "text", "width": 62},
    ]
    prefix = TOOL_PREFIXES
    for tool, items in ESSENTIAL_COUNTS.items():
        if tool == "Rosetta":
            columns += [
                {"key": "rosetta_interface_REU", "label": "Interface (REU)", "section": tool, "kind": "energy"},
                {"key": "rosetta_interface_delta_REU", "label": "Interface delta (REU)", "section": tool, "kind": "delta_energy"},
            ]
        for label, name in items:
            columns.append({"key": f"{prefix[tool]}_{name}", "label": label, "section": tool, "kind": "integer"})
        if tool == "Rosetta":
            for label, name in (
                ("VdW attraction (REU)", "fa_atr_REU"),
                ("VdW repulsion (REU)", "fa_rep_REU"),
                ("Electrostatics (REU)", "fa_elec_REU"),
                ("Solvation (REU)", "solvation_REU"),
                ("hbond_sc only (REU)", "hbond_sc_REU"),
            ):
                columns.append({"key": f"rosetta_{name}", "label": label, "section": tool, "kind": "energy"})
    for key, label, kind in (
        ("global_fv_rmsd_A", "Fv RMSD (A)", "number"),
        ("framework_rmsd_A", "Framework RMSD (A)", "number"),
        ("h_cdr3_rmsd_A", "HCDR3 overall RMSD (A)", "number"),
        ("h_cdr3_anchors_rmsd_A", "HCDR3 anchors RMSD (A)", "number"),
        ("h_cdr3_tip_rmsd_A", "HCDR3 YYGSS tip RMSD (A)", "number"),
        ("l_cdr3_rmsd_A", "LCDR3 RMSD (A)", "number"),
        ("network_added_pair_type_events", "Network events added", "integer"),
        ("network_removed_pair_type_events", "Network events removed", "integer"),
        ("network_changed_unique_pairs_energy_delta_REU", "Network energy delta (REU)", "delta_energy"),
        ("comparison_quality_warning", "Quality warning", "boolean"),
    ):
        columns.append({"key": key, "label": label, "section": "Structure & network", "kind": kind})
    for key, label, kind in ARCHITECTURE_COLUMNS:
        columns.append({"key": key, "label": label, "section": "E2 architecture", "kind": kind})
    from .pyrosetta_score import _legacy_weights
    if _legacy_weights(weights):
        from .rescore_report import report_columns
        return report_columns(columns)
    return columns


def _grouped_row(base: dict[str, str], rows: list[dict[str, object]], parent: bool = False,
                 weights: str = "ref2015") -> dict[str, object]:
    label = "parent" if parent else "mutant"
    output: dict[str, object] = {
        "rank": None if parent else _typed(base.get("rank_by_interface_delta"), "integer"),
        "mutations": "Parent" if parent else base.get("substitutions", ""),
        "candidate": "Parent reference" if parent else base["candidate"],
        "rosetta_interface_REU": _typed(base.get(f"interface_{label}_REU")),
        "rosetta_interface_delta_REU": 0.0 if parent else _typed(base.get("interface_delta_REU")),
    }
    prefix = TOOL_PREFIXES
    for row in rows:
        tool = str(row["tool"])
        for _, name in ESSENTIAL_COUNTS[tool]:
            output[f"{prefix[tool]}_{name}"] = _typed(row[f"{label}_{name}"], "integer")
        if tool == "Rosetta":
            names = ["fa_atr_REU", "fa_rep_REU", "fa_elec_REU", "solvation_REU", "hbond_sc_REU"]
            from .pyrosetta_score import _legacy_weights
            if _legacy_weights(weights):
                names += ["fa_pair_REU", "hbond_bb_sc_REU"]
            for name in names:
                output[f"rosetta_{name}"] = _typed(row.get(f"{label}_{name}"))
    for column in _grouped_columns():
        if column["section"] == "Structure & network":
            key = str(column["key"])
            output[key] = None if parent else _typed(base.get(key), str(column["kind"]))
    return output

def _flatten(rows: list[dict[str, object]]) -> dict[str, object]:
    """Flatten four detector rows into one numeric row for one candidate."""
    if not rows:
        return {}
    output: dict[str, object] = {
        "candidate": rows[0]["candidate"],
        "result_folder": rows[0]["result_folder"],
    }
    tool_names = TOOL_PREFIXES
    for row in rows:
        tool = tool_names[str(row["tool"])]
        for name in COUNT_COLUMNS:
            parent, mutant = row[f"parent_{name}"], row[f"mutant_{name}"]
            if parent == "NA" and mutant == "NA":
                continue
            prefix = f"{tool}_{name}"
            output[f"parent_{prefix}"] = parent
            output[f"mutant_{prefix}"] = mutant
            output[f"delta_{prefix}"] = row[f"delta_{name}"]
    return output

def _markdown(rows: list[dict[str, object]]) -> str:
    shown = COUNT_COLUMNS
    lines = [
        "## Tool-by-tool ligand interaction fingerprints", "",
        "Counts are shown as `parent -> mutant (delta)`. `NaN` means missing data; `NA` means that the tool does not "
        "report that interaction class; it does not mean zero. Counts from different tools use "
        "different definitions and must not be added together. DS pi_stack counts unique protein "
        "residue partners; the complete detector-specific table retains ring-level records.", "",
        "| Tool | Count unit | Interaction class | Parent -> mutant (delta) |", "|---|---|---|---:|",
    ]
    for row in rows:
        for name in shown:
            p, m, d = row[f"parent_{name}"], row[f"mutant_{name}"], row[f"delta_{name}"]
            if p == "NA" and m == "NA":
                continue
            delta = f"{d:+d}" if isinstance(d, int) else d
            p, m, delta = ("NaN" if value is None else value for value in (p, m, delta))
            lines.append(f"| {row['tool']} | {row['count_unit']} | {name} | {p} -> {m} ({delta}) |")
    lines += [
        "", "Rosetta is an energy/scoring method, not a general geometric fingerprint detector. "
        "Its unsupported geometric classes are therefore `NA`; its interface energy components "
        "remain in the energy-decomposition tables.", "",
    ]
    return "\n".join(lines)


def _raw_markdown(rows: list[dict[str, object]]) -> str:
    lines = [
        "### Complete detector-specific counts", "",
        "These are the tools' native labels, retained without forcing unlike definitions together.", "",
        "| Tool | Native interaction type | Parent | Mutant | Delta |", "|---|---|---:|---:|---:|",
    ]
    for row in rows:
        p, m, d = (row[key] for key in ("parent_count", "mutant_count", "delta_count"))
        p, m = ("NaN" if value is None else value for value in (p, m))
        delta = f"{d:+d}" if isinstance(d, int) else "NaN" if d is None else str(d)
        lines.append(
            f"| {row['tool']} | {row['detector_interaction_type']} | {p} | {m} | {delta} |"
        )
    return "\n".join(lines) + "\n"


def _replace_section(path: Path, section: str, heading: str | None = None) -> None:
    name = "e2-architecture" if heading else "interaction-fingerprints"
    begin, end = f"<!-- {name}:begin -->", f"<!-- {name}:end -->"
    text = path.read_text(encoding="utf-8")
    block = f"{begin}\n{section.rstrip()}\n{end}"
    if begin in text and end in text:
        prefix, remainder = text.split(begin, 1)
        _, suffix = remainder.split(end, 1)
        text = prefix.rstrip() + "\n\n" + block + suffix
    elif heading and heading + "\n" in text:
        prefix, remainder = text.split(heading + "\n", 1)
        boundary = re.search(r"(?m)^(?:## |<!-- [^\n]+:begin -->)", remainder)
        suffix = remainder[boundary.start():] if boundary else ""
        text = prefix.rstrip() + "\n\n" + block + "\n\n" + suffix
    else:
        text = text.rstrip() + "\n\n" + block + "\n"
    path.write_text(text, encoding="utf-8")

def generate(batch: Path, cif_folder: Path | None = None, parent_cif: Path | None = None,
             scores_csv: Path | None = None, *, include_empty_metadata_sections: bool = False,
             include_e3_discrimination: bool = False, analysis_options: dict | None = None) -> Path:
    """Rebuild reports from saved artifacts; no scoring or detector runs."""
    batch = Path(batch).resolve()
    options = analysis_options or {}
    manifest_path = batch / "batch_manifest.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8")) if manifest_path.is_file() else {}
    saved_confidence = manifest.get("hcdr3_confidence", {})
    saved_scores = manifest.get("model_scores", {})
    # Late-arriving metadata is picked up by the report subprocess, without
    # stopping an already running scientific controller or editing its manifest.
    source_config_path = batch / "report_sources.json"
    source_config = json.loads(source_config_path.read_text(encoding="utf-8")) if source_config_path.is_file() else {}
    scores_enabled = options.get("boltz2_scores") is not False or options.get("gnina_scores") is not False
    if source_config.get("scores_csv") and not scores_csv and scores_enabled:
        configured_scores = _windows_path(source_config["scores_csv"])
        expected = source_config.get("scores_csv_sha256")
        if expected and sha256(configured_scores.read_bytes()).hexdigest() != expected:
            raise ValueError("Configured model-score source changed since validation")
        scores_csv = configured_scores
    cif_folder = cif_folder or (_windows_path(source_config["cif_folder"]) if source_config.get("cif_folder") else None)
    parent_cif = parent_cif or (_windows_path(source_config["parent_cif"]) if source_config.get("parent_cif") else None)
    scores_csv = scores_csv or (Path(saved_scores["source_csv"]) if saved_scores.get("source_csv") else None)
    cif_folder = cif_folder or (Path(saved_confidence["cif_folder"]) if saved_confidence.get("cif_folder") else None)
    parent_cif = parent_cif or (Path(saved_confidence["parent_cif"]) if saved_confidence.get("parent_cif") else None)
    if options.get("hcdr3_plddt") is False:
        cif_folder = parent_cif = None
    if options.get("boltz2_scores") is False and options.get("gnina_scores") is False:
        scores_csv = None
    if parent_cif and not cif_folder:
        raise ValueError("--parent-cif requires --cif-folder or a previously configured CIF folder")
    parent_ds = batch / "discovery_studio_parent_interactions.csv"
    grouped_rows, numeric_rows, native_rows = [], [], []
    parent_reference = None
    parent_architecture = None
    scoring_weights = None
    e3_saved = {}
    availability = {}
    scoring_legacy = None
    for base in _read_csv(batch / "batch_summary.csv"):
        if base.get("status") != "complete":
            continue
        result_dir = _windows_path(base["result_folder"])
        workflow_path = result_dir / "workflow_manifest.json"
        workflow = json.loads(workflow_path.read_text(encoding="utf-8")) if workflow_path.is_file() else {}
        candidate_weights = workflow.get("scoring_protocol", {}).get("weights", "ref2015")
        from .pyrosetta_score import _legacy_weights
        candidate_legacy = workflow.get("scoring_protocol", {}).get("legacy_pre_talaris", _legacy_weights(candidate_weights))
        if scoring_weights is not None and candidate_weights != scoring_weights:
            raise ValueError("Candidates use mixed scoring weights; refusing a combined report")
        if scoring_legacy is not None and candidate_legacy != scoring_legacy:
            raise ValueError("Candidates use mixed scoring options; refusing a combined report")
        scoring_weights = candidate_weights
        scoring_legacy = candidate_legacy
        report_weights = "ligand" if candidate_legacy else "ref2015"
        automatic_path = result_dir / "automatic_summary.json"
        automatic = json.loads(automatic_path.read_text(encoding="utf-8"))
        automatic["report_analysis_status"] = {key: "not_run" for key in
            ("discovery_studio", "e2_architecture") if options.get(key) is False}
        availability[base["candidate"]] = analysis_status(result_dir, automatic)
        rows = candidate_fingerprints(base["candidate"], result_dir, parent_ds, automatic)
        native = _raw_count_rows(base["candidate"], result_dir, automatic, parent_ds)
        numeric = _flatten(rows)
        reference = _grouped_row(base, rows, parent=True, weights=report_weights)
        if parent_architecture is None:
            parent_architecture = ({"e2_architecture_status": "not_run"} if options.get("e2_architecture") is False
                                   else _saved_architecture(result_dir, "parent", parent_ds))
        reference.update({key: parent_architecture.get(key) for key, _, _ in ARCHITECTURE_COLUMNS})
        if parent_reference is not None and reference != parent_reference:
            raise ValueError("Candidates do not share the same parent fingerprint")
        parent_reference = reference
        grouped = _grouped_row(base, rows, weights=report_weights)
        if availability[base["candidate"]].get("network") == "not_run":
            for key in grouped:
                if key.startswith("network_"):
                    grouped[key] = None
        architecture = ({"e2_architecture_status": "not_run"} if options.get("e2_architecture") is False else
                        _saved_architecture(result_dir, "mutant", result_dir / "discovery_studio_interactions.csv"))
        grouped.update({key: architecture.get(key) for key, _, _ in ARCHITECTURE_COLUMNS})
        grouped_rows.append(grouped)
        if include_e3_discrimination:
            from .e3_reports import append_e3_report, load_saved_e3
            saved_e3 = load_saved_e3(result_dir / "e3_discrimination.json", candidate=base["candidate"])
            if saved_e3 is not None:
                e3_saved[base["candidate"]] = saved_e3
                automatic["e3_discrimination"] = saved_e3
                append_e3_report(result_dir / "automatic_summary.md", saved_e3)
        automatic["estrogen_binding_architecture"] = {
            "parent": parent_architecture,
            "mutant": architecture,
        }
        write_json(result_dir / "estrogen_binding_architecture.json", automatic["estrogen_binding_architecture"])
        numeric_rows.append(numeric)
        native_rows.extend(native)
        write_csv(result_dir / "interaction_fingerprint_counts.csv", [numeric])
        write_json(result_dir / "interaction_fingerprint_counts.json", rows)
        write_csv(result_dir / "interaction_fingerprint_all_detector_counts.csv", native)
        markdown = _markdown(rows) + _raw_markdown(native)
        (result_dir / "interaction_fingerprint_summary.md").write_text(markdown, encoding="utf-8")
        architecture_lines = ["## Conserved E2-binding architecture", "",
                              "The four features are DS O3/O17 conventional H-bonds, calibrated DS C18 "
                              "packing and calibrated coordinate steric geometry. Pi-mediated contacts "
                              "are a separate DS unique-residue descriptor, not an additional pass requirement.", "",
                              "| Metric | Parent | Mutant |", "|---|---|---|"]
        for key, label, _ in ARCHITECTURE_COLUMNS:
            p, m = (item.get(key) for item in (parent_architecture, architecture))
            architecture_lines.append(f"| {label} | {p if p is not None else 'NaN'} | {m if m is not None else 'NaN'} |")
        _replace_section(result_dir / "automatic_summary.md", "\n".join(architecture_lines),
                         heading="## Conserved E2-binding architecture")
        _replace_section(result_dir / "automatic_summary.md", markdown)
        automatic["tool_interaction_fingerprints"] = rows
        automatic["tool_interaction_fingerprint_native_counts"] = native
        write_json(automatic_path, automatic)
    if not grouped_rows:
        raise ValueError("No completed candidates were found in batch_summary.csv")
    grouped_path = batch / "batch_summary_grouped_fingerprints.csv"
    numeric_path = batch / "batch_interaction_fingerprint_numeric.csv"
    payload_path = batch / "batch_summary_grouped_fingerprints.json"
    write_csv(grouped_path, grouped_rows)
    write_csv(batch / "e2_architecture_screen.csv", [
        {"candidate": row["candidate"], **{key: row.get(key) for key, _, _ in ARCHITECTURE_COLUMNS}}
        for row in grouped_rows
    ])
    write_csv(numeric_path, numeric_rows)
    write_csv(batch / "batch_interaction_fingerprint_all_detector_counts.csv", native_rows)
    from .pyrosetta_score import _legacy_weights
    legacy = bool(scoring_legacy)
    payload = {
        "title": "Fv-ligand RosettaLigand analysis" if legacy else "Fv-ligand interaction fingerprints",
        "note": f"{len(grouped_rows)} completed analyses only; failed/incomplete inputs remain in batch_summary.csv. The pinned row gives the parent baseline.",
        "columns": _grouped_columns("ligand" if legacy else "ref2015"),
        "rows": grouped_rows,
        "parent": parent_reference,
        "source_note": "Source: saved Rosetta, Arpeggio, PLIP and Discovery Studio analyses. "
                       "Counts follow each tool's native definitions. Arpeggio ring-ring contacts "
                       "are plane-plane records, distinct from aromatic atom contacts. "
                       "Energy deltas are mutant minus parent. E2 architecture: O3/O17 H-bonds and calibrated C18 contacts "
                       "use DS only; steric sandwich uses calibrated coordinate geometry. Pi-mediated contacts include DS stacked, T-shaped, pi-cation "
                       "and pi-anion contacts; each partner residue is counted once across classes. "
                       "The all-four screen requires O3, O17, C18 and steric features, not pi-mediated contacts. "
                       "Numerical proxies are project-calibrated, not published affinity thresholds. Missing DS data leave interaction metrics unavailable.",
        "scoring": {"weights": scoring_weights, "legacy_pre_talaris": legacy},
        "analysis_availability": availability,
    }
    from .ds_architecture import VERSION as architecture_version
    from .steric_sandwich import METHOD_VERSION as steric_version
    payload["e2_recognition_protocol"] = {
        "definition_version": architecture_version,
        "interaction_source": "Discovery Studio",
        "steric_method": steric_version,
        "maximum_surface_gap_A": 1.0,
        "c19_obstruction_included": False,
    }
    if legacy:
        payload["scoring_profile"] = "rosettaligand_legacy"
        payload["source_note"] = ("Rosetta scores and mutation-network energies use the recorded weights "
                                  "with matching legacy pre-Talaris options. " + payload["source_note"])
    disabled = sorted({name for status in availability.values() for name, state in status.items()
                       if state == "not_run"} | {name for name, enabled in options.items() if enabled is False})
    if disabled:
        payload["note"] += " Not run: " + ", ".join(disabled) + ". Blank cells do not mean zero."
    if include_e3_discrimination:
        from .e3_reports import append_e3_to_payload, load_saved_e3
        parent_e3 = load_saved_e3(batch / "parent_e3_discrimination.json")
        payload = append_e3_to_payload(payload, e3_saved, parent_e3)
        write_csv(grouped_path, payload["rows"])
    write_json(payload_path, payload)
    shutil.copyfile(resource_path("docs/interaction_fingerprint_columns.md"), batch / "INTERACTION_FINGERPRINT_COLUMNS.md")
    shutil.copyfile(resource_path("docs/e2_architecture.md"), batch / "E2_DEFINITIONS.md")
    if manifest_path.is_file():
        manifest["interaction_fingerprints"] = {
            "grouped_summary_csv": str(grouped_path),
            "grouped_summary_json": str(payload_path),
            "numeric_flat_csv": str(numeric_path),
        }
        write_json(manifest_path, manifest)
    if cif_folder:
        from .confidence import add_confidence

        grouped_path = add_confidence(batch, cif_folder, parent_cif)
    if scores_csv:
        from .model_scores import add_model_scores

        grouped_path = add_model_scores(batch, scores_csv)
    if include_empty_metadata_sections:
        # Keep the established workbook sections without inventing absent scores.
        from .confidence import confidence_columns
        from .model_scores import score_columns
        payload = json.loads(payload_path.read_text(encoding="utf-8"))
        extra = confidence_columns() + score_columns()
        # Metadata may arrive independently: keep the established section order
        # even when model scores exist but original CIF confidence is missing.
        payload["columns"] = [column for column in payload["columns"]
                              if column["section"] not in {"HCDR3 pLDDT", "Boltz-2", "GNINA"}] + extra
        for row in [payload["parent"], *payload["rows"]]:
            for column in extra:
                row.setdefault(column["key"], None)
        available_score_fields = set(payload.get("model_scores", {}).get("available_fields", []))
        boltz_fields = {column["source_column"] for column in score_columns() if column["section"] == "Boltz-2"}
        gnina_fields = {column["source_column"] for column in score_columns() if column["section"] == "GNINA"}
        def availability(fields):
            if not scores_csv:
                return "unavailable: matching score table not supplied"
            present = len(available_score_fields & fields)
            if present == len(fields):
                return "available"
            return f"partial: {present}/{len(fields)} source fields supplied" if present else "unavailable: these fields are absent from the supplied score table"
        payload["metadata_availability"] = {
            "HCDR3 pLDDT": "available" if cif_folder else "unavailable: matching original CIFs not supplied",
            "Boltz-2": availability(boltz_fields),
            "GNINA": availability(gnina_fields),
        }
        if not cif_folder:
            payload["confidence_note"] = "HCDR3 pLDDT is unavailable: matching original CIF files were not supplied. Cells are blank, not zero."
        if not scores_csv:
            payload["model_scores_note"] = "Boltz-2 and GNINA scores are unavailable: a matching score table was not supplied. Cells are blank, not zero."
        else:
            payload["model_scores_note"] = (
                f"Model scores imported unchanged from {Path(scores_csv).name}. "
                f"Boltz-2: {payload['metadata_availability']['Boltz-2']}. "
                f"GNINA: {payload['metadata_availability']['GNINA']}. "
                "Absent values and unmatched parent-reference scores remain blank, not zero.")
        write_json(payload_path, payload)
        write_csv(grouped_path, payload["rows"])
    disabled_metadata = {section for key, section in (("hcdr3_plddt", "HCDR3 pLDDT"),
                        ("boltz2_scores", "Boltz-2"), ("gnina_scores", "GNINA"))
                        if options.get(key) is False}
    if disabled_metadata:
        payload = json.loads(payload_path.read_text(encoding="utf-8"))
        for column in payload["columns"]:
            if column["section"] in disabled_metadata:
                for row in [payload["parent"], *payload["rows"]]:
                    row[column["key"]] = None
        payload.setdefault("metadata_availability", {}).update({section: "not_run" for section in disabled_metadata})
        write_json(payload_path, payload)
        write_csv(grouped_path, payload["rows"])
    return grouped_path
