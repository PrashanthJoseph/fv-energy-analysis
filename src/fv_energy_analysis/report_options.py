"""Project saved results into one requested CSV/XLSX view without recomputation."""

from __future__ import annotations

from copy import deepcopy
import csv
from hashlib import sha256
import json
from pathlib import Path

from .report import write_json


IDENTITY = {"rank", "mutations", "candidate"}


def _list(value, label):
    if value is None:
        return None
    if isinstance(value, str):
        value = [item.strip() for item in value.split(",") if item.strip()]
    if not isinstance(value, (list, tuple)) or not all(isinstance(item, str) and item for item in value):
        raise ValueError(f"{label} must be a list of nonempty names")
    if len(set(value)) != len(value):
        raise ValueError(f"Duplicate {label}")
    return list(value)


def project_payload(payload, sections=None, columns=None):
    """Keep canonical column order and all rows; Identity is always retained.

    Sections and individual column keys form a union when both are supplied.
    Omitting both selects the complete saved schema. Empty lists select Identity.
    """
    sections, columns = _list(sections, "sections"), _list(columns, "columns")
    available = list(dict.fromkeys(column["section"] for column in payload["columns"]))
    by_key = {column["key"]: column for column in payload["columns"]}
    normalize = lambda value: "".join(c for c in value.casefold() if c.isalnum())
    aliases = {normalize(name): name for name in available}
    for alias, canonical in (("ds", "Discovery Studio"), ("structure", "Structure & network"),
                             ("network", "Structure & network"), ("e2", "E2 architecture"),
                             ("e3", "E3 discrimination"), ("plddt", "HCDR3 pLDDT"),
                             ("status", "Analysis status")):
        if canonical in available:
            aliases[alias] = canonical
    for name in ("Rosetta", "RosettaLigand"):
        if name in available:
            aliases["rosetta"] = aliases["rosettaligand"] = name
    selected_sections = []
    for name in sections or []:
        if normalize(name) == "all":
            selected_sections = available.copy()
            continue
        actual = aliases.get(normalize(name))
        if actual is None:
            raise ValueError(f"Unknown report section {name!r}; available: {', '.join(available)}")
        if actual not in selected_sections:
            selected_sections.append(actual)
    unknown = set(columns or []) - by_key.keys()
    if unknown:
        raise ValueError(f"Unknown report columns: {', '.join(sorted(unknown))}")
    selected = set(by_key) if sections is None and columns is None else IDENTITY | set(columns or [])
    selected.update(c["key"] for c in payload["columns"] if c["section"] in selected_sections)
    missing_identity = IDENTITY - by_key.keys()
    if missing_identity:
        raise ValueError(f"Source report lacks Identity fields: {sorted(missing_identity)}")
    output = deepcopy(payload)
    output["columns"] = [c for c in payload["columns"] if c["key"] in selected]
    output["rows"] = [{c["key"]: row.get(c["key"]) for c in output["columns"]} for row in payload["rows"]]
    output["parent"] = {c["key"]: payload.get("parent", {}).get(c["key"]) for c in output["columns"]}
    output["report_projection"] = {"sections": sections, "columns": columns,
                                   "identity_always_included": True,
                                   "selection_mode": "union of sections and column keys"}
    return output


def validate_report_selection(config):
    """Validate configured output knobs before running any scientific stage."""
    from .confidence import confidence_columns
    from .fingerprint_reports import _grouped_columns
    from .model_scores import score_columns
    from .pyrosetta_score import _legacy_weights

    options = config.get("analyses") or {}
    recorded = config.get("weights_identity") or {}
    legacy = recorded.get("legacy_pre_talaris")
    if legacy is None:
        legacy = config.get("legacy_pre_talaris")
    if legacy is None:
        legacy = _legacy_weights(config.get("weights", "ligand"))
    columns = _grouped_columns("ligand" if legacy else "ref2015")
    keep_metadata = config.get("include_empty_metadata_sections", True)
    if keep_metadata or config.get("cif_folder") and options.get("hcdr3_plddt") is not False:
        columns += confidence_columns()
    if keep_metadata or config.get("scores_csv"):
        columns += score_columns()
    e3 = config.get("e3_discrimination", options.get("e3_discrimination", False))
    if e3:
        from .e3_reports import e3_columns
        columns += e3_columns()
    if config.get("combine_sources") or e3:
        columns += [{"key": key, "label": label, "section": "Analysis status", "kind": "text"}
                    for key, label in (("status", "Analysis status"), ("error", "Failure reason"))]
    project_payload({"columns": columns, "parent": {}, "rows": []},
                    config.get("report_sections"), config.get("report_columns"))
    formats = _list(config.get("report_formats", ("csv", "xlsx")), "formats")
    if not formats or set(formats) - {"csv", "xlsx"}:
        raise ValueError("Report formats must contain csv, xlsx, or both")


def prepare_selected_report(folder):
    """Validate a derived view against the exact canonical payload it references."""
    from .paths import native_path
    from .workbook_validation import prepare, typed

    folder = Path(folder).resolve()
    manifest = json.loads((folder / "selected_report_manifest.json").read_text(encoding="utf-8"))
    source = native_path(manifest["source_batch"])
    source_file = source / "batch_summary_grouped_fingerprints.json"
    if sha256(source_file.read_bytes()).hexdigest() != manifest["source_payload_sha256"]:
        raise ValueError("The selected report's canonical source data changed")
    canonical, source_manifest, failures = prepare(source)
    for key in ("candidate_count", "completed_count", "failure_count"):
        if manifest[key] != source_manifest.get(key, 0):
            raise ValueError(f"Selected-report {key} differs from the canonical manifest")
    expected = project_payload(canonical, manifest.get("sections"), manifest.get("columns"))
    payload = json.loads((folder / "batch_summary_grouped_fingerprints.json").read_text(encoding="utf-8"))
    if payload != expected:
        raise ValueError("Selected report differs from its canonical values or selection")
    csv_file = folder / "batch_summary_grouped_fingerprints.csv"
    if "csv" in manifest["formats"]:
        if sha256(csv_file.read_bytes()).hexdigest() != manifest["csv_sha256"]:
            raise ValueError("Selected CSV changed after export")
        with csv_file.open(encoding="utf-8-sig", newline="") as stream:
            reader = csv.DictReader(stream)
            records = list(reader)
            if reader.fieldnames != [c["key"] for c in payload["columns"]] or len(records) != len(payload["rows"]):
                raise ValueError("Selected CSV layout differs from its canonical projection")
        for saved, actual in zip(payload["rows"], records):
            if any(typed(saved.get(c["key"]), c) != typed(actual.get(c["key"]), c) for c in payload["columns"]):
                raise ValueError("Selected CSV values differ from the canonical projection")
    return payload, manifest, failures


def export_selected_report(batch, sections=None, columns=None, formats=("csv", "xlsx"), output=None):
    """Write one explicitly requested view; never rewrite canonical results."""
    from .workbook_validation import prepare, verify
    from .spreadsheet import export_workbook

    batch = Path(batch).resolve()
    formats = _list(formats, "formats")
    if not formats or set(formats) - {"csv", "xlsx"}:
        raise ValueError("Report formats must contain csv, xlsx, or both")
    canonical, source_manifest, failures = prepare(batch)
    payload = project_payload(canonical, sections, columns)
    destination = Path(output).resolve() if output is not None else batch / "selected_report"
    if destination == batch:
        raise ValueError("Selected reports must use a separate output folder")
    destination.mkdir(parents=True, exist_ok=True)
    manifest_file = destination / "selected_report_manifest.json"
    if any(destination.iterdir()) and not manifest_file.is_file():
        raise FileExistsError("Selected-report output contains unrelated files")
    manifest = {"status": "complete", "source_batch": str(batch),
                "source_payload_sha256": sha256((batch / "batch_summary_grouped_fingerprints.json").read_bytes()).hexdigest(),
                "sections": _list(sections, "sections"), "columns": _list(columns, "columns"),
                "formats": formats, "candidate_count": source_manifest["candidate_count"],
                "completed_count": source_manifest["completed_count"], "failure_count": source_manifest.get("failure_count", 0),
                "scientific_analysis_rerun": False}
    # A repeated export of the same view is safe; reusing a folder for a
    # different request requires a new destination to preserve its provenance.
    if manifest_file.is_file():
        previous = json.loads(manifest_file.read_text(encoding="utf-8"))
        for field in ("source_batch", "source_payload_sha256", "sections", "columns", "formats"):
            if previous[field] != manifest[field]:
                raise ValueError("Selected-report destination belongs to a different source or selection")
    write_json(destination / "batch_summary_grouped_fingerprints.json", payload)
    if "csv" in formats:
        csv_file = destination / "batch_summary_grouped_fingerprints.csv"
        with csv_file.open("w", encoding="utf-8", newline="") as stream:
            writer = csv.DictWriter(stream, fieldnames=[column["key"] for column in payload["columns"]])
            writer.writeheader()
            writer.writerows(payload["rows"])
        manifest["csv_sha256"] = sha256(csv_file.read_bytes()).hexdigest()
    write_json(manifest_file, manifest)
    if "xlsx" in formats:
        workbook = destination / "batch_summary_grouped_fingerprints.xlsx"
        export_workbook(payload, workbook, failures=failures)
        verify(destination, workbook)
    return destination
