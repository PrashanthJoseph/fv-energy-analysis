"""Reporting-only views of saved E3 fixed-pocket steric assessments.

Nothing here imports Rosetta, constructs a ligand or recalculates a contact.
Unassessed and failed measurements stay missing rather than becoming zero.
"""
from __future__ import annotations

import copy
import json
import math
from numbers import Real
from pathlib import Path


SECTION = "E3 discrimination"
E3_COLUMNS = (
    ("e3_status", "E3 assessment status", "text"),
    ("e3_added_o_closest_distance_A", "16alpha-O nearest distance (A)", "number"),
    ("e3_added_o_closest_residue", "16alpha-O nearest residue", "text"),
    ("e3_added_o_deepest_overlap_A", "Maximum vdW overlap (A)", "number"),
    ("e3_added_o_overlap_residue_count", "Large-overlap residues", "integer"),
    ("e3_added_o_non_hbond_overlap_residue_count", "Non-Hbond overlap residues", "integer"),
    ("e3_added_o_hbond_supported_overlap_residue_count", "Hbond-supported overlap residues", "integer"),
    ("e3_added_o_unassessed_overlap_residue_count", "Unassessed overlap residues", "integer"),
    ("e3_added_o_overlap_residues", "Large-overlap residue identities", "text"),
    ("e3_steric_assessment", "Fixed-pocket assessment", "text"),
    ("e3_steric_warning", "Heuristic steric warning", "boolean"),
    ("e3_oh_h_best_sampled_protein_overlap_A", "Best OH-H nonpolar protein overlap (A)", "number"),
)


def e3_columns() -> list[dict]:
    return [{"key": key, "label": label, "kind": kind, "section": SECTION}
            for key, label, kind in E3_COLUMNS]


def _typed(value, kind, key):
    if value is None:
        return None
    if kind == "text":
        if not isinstance(value, str):
            raise ValueError(f"Expected text for {key}")
        return value
    if kind == "boolean":
        if not isinstance(value, bool):
            raise ValueError(f"Expected a boolean E3 descriptor: {key}")
        return value
    if isinstance(value, bool) or not isinstance(value, Real) or not math.isfinite(value):
        raise ValueError(f"Expected a finite numeric E3 measurement: {key}")
    if kind == "integer":
        if not float(value).is_integer() or value < 0:
            raise ValueError(f"Expected a nonnegative integer count: {key}")
        return int(value)
    return float(value)


def flatten_e3_result(result: dict | None, *, missing_status: str = "not_assessed") -> dict:
    """Expose saved steric descriptors without assigning an affinity verdict."""
    output = {key: None for key, _, _ in E3_COLUMNS}
    if result is None:
        output["e3_status"] = missing_status
        return output
    if result.get("schema_version", 1) not in {1, 2}:
        raise ValueError("Unsupported saved E3-discrimination schema version")
    status = result.get("status")
    if status not in {"complete", "failed", "not_assessed", "excluded_input_rejection"}:
        raise ValueError(f"Unrecognized saved E3 assessment status: {status!r}")
    output["e3_status"] = status
    if status != "complete":
        return output
    summary = result.get("summary")
    if not isinstance(summary, dict):
        raise ValueError("A completed E3 assessment must contain a summary object")
    for key, _, kind in E3_COLUMNS:
        if key != "e3_status":
            output[key] = _typed(summary.get(key), kind, key)
    for key in ("e3_added_o_closest_distance_A", "e3_added_o_deepest_overlap_A"):
        if output[key] is None:
            raise ValueError(f"A completed E3 assessment lacks its fixed-pocket measurement: {key}")
    if summary.get("e3_analysis_scope") != "fixed_pocket_steric_only":
        raise ValueError("This report requires an explicitly fixed-pocket steric-only E3 assessment")
    assessment = output["e3_steric_assessment"]
    if assessment not in {"steric_overlap_detected", "hbond_supported_overlaps_only",
                           "no_large_overlap_in_fixed_pose", "unassessed_polar_geometry"}:
        raise ValueError("A completed E3 assessment requires an explicit steric assessment category")
    counts = [output[k] for k in ("e3_added_o_overlap_residue_count",
                                  "e3_added_o_non_hbond_overlap_residue_count",
                                  "e3_added_o_hbond_supported_overlap_residue_count",
                                  "e3_added_o_unassessed_overlap_residue_count")]
    if any(v is None for v in counts):
        raise ValueError("Current E3 reporting requires geometry-resolved overlap counts")
    total, non_hbond, supported, unassessed = counts
    if any(v > total for v in (non_hbond, supported, unassessed)):
        raise ValueError("Each E3 overlap-residue category must be a subset of all overlapping residues")
    expected_assessment, expected_warning = (
        ("steric_overlap_detected", True) if non_hbond else
        ("unassessed_polar_geometry", None) if unassessed else
        ("hbond_supported_overlaps_only", False) if total and supported else
        ("no_large_overlap_in_fixed_pose", False) if not total else (None, None))
    if expected_assessment is None or (assessment, output["e3_steric_warning"]) != (expected_assessment, expected_warning):
        raise ValueError("E3 assessment and three-state warning disagree with geometry-resolved overlap counts")
    distance = output["e3_added_o_closest_distance_A"]
    if distance is not None and distance < 0:
        raise ValueError("Nearest-atom distance cannot be negative")
    hydrogen_overlap = output["e3_oh_h_best_sampled_protein_overlap_A"]
    if hydrogen_overlap is not None and hydrogen_overlap < 0:
        raise ValueError("Sampled OH-H overlap cannot be negative")
    return output


def load_saved_e3(path: str | Path, *, candidate: str | None = None) -> dict | None:
    """Read and validate an existing artifact; absent files are unassessed."""
    path = Path(path)
    if not path.is_file():
        return None
    result = json.loads(path.read_text(encoding="utf-8-sig"))
    if candidate is not None and result.get("candidate") not in (None, candidate):
        raise ValueError(f"E3 result belongs to a different candidate: {path}")
    flatten_e3_result(result)
    return result


def append_e3_to_payload(payload: dict, results_by_candidate: dict[str, dict],
                         parent_result: dict | None = None) -> dict:
    """Return a new grouped payload, preserving order, rank and prior metrics."""
    output = copy.deepcopy(payload)
    names = [row["candidate"] for row in output["rows"]]
    if len(set(names)) != len(names):
        raise ValueError("Duplicate candidates in the source grouped payload")
    unknown = set(results_by_candidate) - set(names)
    if unknown:
        raise ValueError(f"E3 results are not present in the source dataset: {sorted(unknown)}")
    previous_e3_keys = {c["key"] for c in output["columns"] if c["section"] == SECTION}
    current_e3_keys = {key for key, _, _ in E3_COLUMNS}
    obsolete_e3_keys = previous_e3_keys - current_e3_keys
    columns = [c for c in output["columns"] if c["section"] != SECTION]
    # Keep the new screen beside the existing recognition section; the four
    # detector sections and confidence/model columns retain their order.
    insert = next((i + 1 for i in range(len(columns)-1, -1, -1)
                   if columns[i]["section"] == "E2 architecture"), len(columns))
    output["columns"] = columns[:insert] + e3_columns() + columns[insert:]
    for row in output["rows"]:
        for key in obsolete_e3_keys:
            row.pop(key, None)
        saved = results_by_candidate.get(row["candidate"])
        if saved is not None and saved.get("candidate") not in (None, row["candidate"]):
            raise ValueError(f"Mismatched saved E3 candidate: {row['candidate']}")
        missing = "excluded_input_rejection" if row.get("status") == "failed" else "not_assessed"
        row.update(flatten_e3_result(saved, missing_status=missing))
    if output.get("parent") is not None:
        for key in obsolete_e3_keys:
            output["parent"].pop(key, None)
        output["parent"].update(flatten_e3_result(parent_result))
    output["e3_discrimination_protocol"] = {
        "schema_version": 2,
        "input": "saved E3-discrimination artifacts only; no computation during reporting",
        "substitution": "estriol adds 16alpha-OH; C18 is unchanged",
        "analysis_scope": "fixed_pocket_steric_only",
        "coordinate_policy": "shared E2/E3 heavy atoms and all source protein hydrogens are fixed; neutral full-hydrogen E3 probes are constructed",
        "added_oh_hydrogen_sampling": "one consistent saved OH-H orientation for the complete overlap network; 12-torsion sampling is not a continuous minimum or PyRosetta OptH",
        "energies_assessed": False,
        "local_accommodation_assessed": False,
        "interpretation": "mechanistic descriptors, not a calibrated selectivity probability or affinity prediction",
        "warning_interpretation": "true: definite non-Hbond-supported overlap; false: fully assessed without such overlap; null: unresolved geometry without a definite clash",
        "polar_overlap_policy": "N/O/S overlaps require chemically typed explicit-H donor/acceptor geometry; element identity alone never exempts an overlap",
        "overlap_classification": "geometry_resolved_same_oh_orientation",
        "hydrogen_bond_geometry_assessed": True,
        "hydrogen_bond_energies_assessed": False,
        "numerical_specificity_threshold": None,
        "source_references": ["https://pmc.ncbi.nlm.nih.gov/articles/PMC305588/"],
    }
    return output


def _geometry_markdown(contacts: list[dict]) -> list[str]:
    """Show saved trials, including failures; never infer absent geometry."""
    overlapping = [c for c in contacts if c.get("overlap_at_or_above_warning")]
    if not overlapping:
        return []
    lines = ["### Geometry of overlapping contacts", "",
             "These rows are donor-H-acceptor geometry trials for the saved overlap contacts, "
             "not counts of independent hydrogen bonds. Every trial uses the same chosen OH-H "
             "orientation. A dash means that a measurement was not applicable or was not recorded.", "",
             "| Overlap atom | Geometry status | Donor | H | Acceptor | D-A (A) | H-A (A) | D-H-A (deg) | Geometry passed | Reason |",
             "|---|---|---|---|---|---:|---:|---:|---|---|"]
    def cell(value, decimals=None):
        if value is None:
            return "—"
        if decimals is not None:
            return f"{float(value):.{decimals}f}"
        return str(value).replace("\n", " ").replace("|", r"\|")
    for contact in overlapping:
        geometry = contact.get("hbond_geometry") or {}
        status = geometry.get("status", geometry.get("state", "unassessed"))
        reason = geometry.get("reason", "No donor-H-acceptor trial was recorded")
        atom = contact.get("protein_atom", contact.get("protein_residue"))
        for trial in geometry.get("trials", []) or [{}]:
            values = [cell(atom), cell(str(status).replace("_", " ")),
                      cell(trial.get("donor")), cell(trial.get("hydrogen")), cell(trial.get("acceptor")),
                      cell(trial.get("donor_acceptor_distance_A"), 3),
                      cell(trial.get("hydrogen_acceptor_distance_A"), 3),
                      cell(trial.get("angle_DHA_deg"), 1), cell(trial.get("passes_geometry")),
                      cell(trial.get("reason", reason))]
            lines.append("| " + " | ".join(values) + " |")
    return lines + [""]


def render_e3_markdown(result: dict | None) -> str:
    values = flatten_e3_result(result)
    lines = ["## E3 discrimination", "",
             "Estriol adds a 16alpha-OH group; its C18 methyl is unchanged. This section tests "
             "direction-specific steric accommodation of the added oxygen in a fixed pocket. "
             "It is separate from the four E2 recognition features and does not change their score.", "",
             "Neutral estriol probes contain all 24 ligand hydrogens (45 atoms including the 21 heavy atoms). "
             "All source protein hydrogens are retained. Protein and shared steroid heavy-atom coordinates "
             "remain fixed. Oxygen distances and overlap sizes are independent of hydroxyl-H orientation; "
             "polar overlap classification additionally uses explicit-H donor/acceptor geometry. "
             "This is a geometric assessment: "
             "no Rosetta energies, local repacking or minimization are "
             "computed. A fixed-pocket obstruction may be relieved by ligand or protein movement.", "",
             "| Metric | Value |", "|---|---:|"]
    for key, label, kind in E3_COLUMNS:
        value = values[key]
        shown = "Not assessed" if value is None else f"{value:.4f}" if kind == "number" else str(value)
        if key == "e3_steric_assessment":
            shown = shown.replace("_", " ")
        lines.append(f"| {label} | {shown.replace('|', chr(92) + '|')} |")
    if result and result.get("error"):
        lines += ["", "Assessment error: " + str(result["error"]).replace("\n", " ")]
    if result:
        lines += [""] + _geometry_markdown(result.get("contacts", []))
    lines += ["", "Steric evidence comes from explicit 16alpha-OH placement and protein atom identities. "
              "The complete contact list retains distances and van der Waals surface gaps. A short distance or negative "
              "van der Waals surface gap alone is not classified as a clash: legitimate H-bonds can also "
              "have such gaps. The nearest-residue field is a neighbour identity, not a clash verdict.", "",
              "The steric warning is a project-defined geometric heuristic, not a literature-validated "
              "E3 rejection threshold. Large-overlap counts use the configured surface-overlap "
              "threshold (default 0.4 A). Element identity alone does not exempt N/O/S overlaps. "
              "Chemically valid donor/acceptor typing and actual donor-H-acceptor geometry must support "
              "an H-bond explanation; a tested polar overlap failing those criteria contributes to the "
              "non-Hbond overlap category. Missing typing or hydrogen geometry is explicitly unassessed. "
              "The warning is TRUE when a definite non-Hbond-supported overlap is present, FALSE only "
              "when all relevant contacts have been assessed without such an overlap, and blank when "
              "geometry remains unresolved without a definite clash. This section imposes no binary "
              "specificity or affinity verdict.", "",
              "The geometric H-bond convention requires a typed donor/acceptor pair, D-A <4.1 A "
              "and D-H-A >100 degrees ([PLIP thresholds](https://github.com/pharmai/plip/blob/master/plip/basic/config.py), "
              "[detection code](https://github.com/pharmai/plip/blob/master/plip/structure/detection.py)). "
              "The additional H-A >=1.2 A safeguard is a project choice to prevent very short "
              "contacts from being excused as hydrogen bonds; it is not a PLIP threshold. These "
              "custom geometry results are separate from native DS, PLIP and Rosetta counts.", "",
              "When available, the separate OH-H descriptor samples 12 hydroxyl torsions. For each "
              "orientation it measures the largest non-H-bond-compatible protein overlap, then reports the "
              "smallest of those 12 maxima. This sweep descriptor can differ from the single "
              "orientation selected for the complete contact network. It uses a project idealized "
              "polar-H radius of 1.0 A. "
              "The detailed result retains all orientations and contacts. This is a finite torsion "
              "screen, not PyRosetta OptH, an energy minimum or a MolProbity clashscore; an obstruction "
              "in one arbitrary OH-H orientation does not determine the oxygen-based warning.", "",
              "All polar-overlap classifications in one report use the same chosen OH-H orientation. "
              "The choice minimizes certain clashes, then unresolved contacts, and retains the persisted "
              "orientation on ties. H-bond support is never pooled across mutually incompatible OH "
              "torsions. The exported probe records the selected saved H25 position so its visualized "
              "geometry agrees with the assessment.", "",
              "Literature precedent: the [SHBG structural study](https://pmc.ncbi.nlm.nih.gov/articles/PMC305588/) "
              "proposes steric hindrance of estriol's 16alpha-OH near the Asp65-Trp66 peptide bond. "
              "That is evidence for a possible mechanism in SHBG, not proof of the same mechanism in m7. "
              "Energetic discrimination, including the cost of an unsatisfied extra OH, is unassessed "
              "in this steric-only stage. The placement and numerical descriptors here are the project's "
              "operational protocol, not E3-specific cutoffs published by the SHBG paper.", ""]
    return "\n".join(lines)


def append_e3_report(path: str | Path, result: dict) -> None:
    """Add/replace only the E3 block, leaving existing report sections intact."""
    path = Path(path)
    text = path.read_text(encoding="utf-8") if path.is_file() else ""
    begin, end = "<!-- e3-discrimination:begin -->", "<!-- e3-discrimination:end -->"
    block = begin + "\n" + render_e3_markdown(result).rstrip() + "\n" + end
    if begin in text and end in text:
        prefix, rest = text.split(begin, 1)
        _, suffix = rest.split(end, 1)
        text = prefix.rstrip() + "\n\n" + block + suffix
    else:
        text = text.rstrip() + "\n\n" + block + "\n"
    path.write_text(text, encoding="utf-8")
