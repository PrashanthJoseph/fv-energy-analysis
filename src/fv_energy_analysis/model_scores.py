"""Join saved Boltz-2 and GNINA scores to existing reports without rescoring."""

from __future__ import annotations

import argparse
import csv
import json
import math
import re
from pathlib import Path

from .fingerprint_counts import _read_csv, _windows_path
from .report import write_csv, write_json


BOLTZ2_FIELDS = (
    "confidence_score", "ptm", "iptm", "ligand_iptm", "protein_iptm",
    "complex_plddt", "complex_iplddt", "complex_pde", "complex_ipde",
    "affinity_pred_value", "affinity_probability_binary", "affinity_pred_value1",
    "affinity_probability_binary1", "affinity_pred_value2", "affinity_probability_binary2",
)
BOLTZ2_AFFINITY_FIELDS = BOLTZ2_FIELDS[9:]
GNINA_FIELDS = tuple(
    f"{pose}_{metric}" for pose in ("rank1", "best")
    for metric in ("minimizedAffinity", "CNNscore", "CNNaffinity", "CNN_VS", "CNNaffinity_variance")
)
SECTIONS = ("Boltz-2", "GNINA")


def score_columns() -> list[dict[str, object]]:
    """Summary schema; source names and score scales are deliberately retained."""
    labels = {
        "confidence_score": "Confidence score", "ptm": "pTM", "iptm": "iPTM",
        "ligand_iptm": "Ligand iPTM", "protein_iptm": "Protein iPTM",
        "complex_plddt": "Complex pLDDT (0-1)", "complex_iplddt": "Complex iPLDDT (0-1)",
        "complex_pde": "Complex PDE", "complex_ipde": "Complex iPDE",
        "affinity_pred_value": "Affinity value", "affinity_probability_binary": "Binder probability",
        "affinity_pred_value1": "Affinity value 1", "affinity_probability_binary1": "Binder probability 1",
        "affinity_pred_value2": "Affinity value 2", "affinity_probability_binary2": "Binder probability 2",
    }
    for pose in ("rank1", "best"):
        for metric, label in (("minimizedAffinity", "minimizedAffinity"), ("CNNscore", "CNN score"),
                              ("CNNaffinity", "CNN affinity"), ("CNN_VS", "CNN VS"),
                              ("CNNaffinity_variance", "CNN affinity variance")):
            labels[f"{pose}_{metric}"] = f"{'Rank 1' if pose == 'rank1' else 'Best'} {label}"
    return [
        {"key": f"{prefix}_{field}",
         "label": labels[field],
         "section": section, "kind": "number", "source_column": field}
        for prefix, section, fields in (("boltz2", "Boltz-2", BOLTZ2_FIELDS), ("gnina", "GNINA", GNINA_FIELDS))
        for field in fields
    ]


def _model_key(value: str) -> str:
    """Exact basename identity apart from extension and a terminal relaxed tag."""
    name = str(value).strip().replace("\\", "/").rsplit("/", 1)[-1]
    lowered = name.casefold()
    for suffix in (".mmcif", ".cif", ".pdb"):
        if lowered.endswith(suffix):
            name = name[:-len(suffix)]
            break
    return name.removesuffix("_relaxed")


def _score(value: str | None, field: str, identity: str) -> float | None:
    if value is None or not value.strip():
        return None
    try:
        score = float(value)
    except ValueError as error:
        raise ValueError(f"Invalid numeric score {field}={value!r} for {identity}") from error
    if not math.isfinite(score):
        raise ValueError(f"Nonfinite score {field}={value!r} for {identity}")
    return score


def _source_index(path: Path) -> dict[str, dict[str, object]]:
    with path.open(newline="", encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle)
        fields = reader.fieldnames or []
        native = "job_name" in fields and not {"cif_file", "input"} & set(fields)
        required = BOLTZ2_AFFINITY_FIELDS if native else BOLTZ2_FIELDS + GNINA_FIELDS
        missing = set(required) - set(fields)
        if missing:
            raise ValueError(f"Score CSV is missing required columns: {', '.join(sorted(missing))}")
        if not native and not {"cif_file", "input"} & set(fields):
            raise ValueError("Score CSV requires cif_file or input for exact model matching")
        available = [column["source_column"] for column in score_columns() if column["source_column"] in fields]
        index: dict[str, dict[str, object]] = {}
        for line, row in enumerate(reader, 2):
            if native:
                if None in row:
                    raise ValueError(f"Malformed score CSV row {line}: more values than header columns")
                if not any(value and value.strip() for value in row.values()):
                    continue
                identity = (row.get("job_name") or "").strip()
                if not identity or _model_key(identity) != identity:
                    raise ValueError(f"Missing or invalid job_name in score CSV row {line}")
                if "result_status" in fields and (row.get("result_status") or "").strip().casefold() != "success":
                    raise ValueError(f"Unsuccessful result_status for {identity} in score CSV row {line}")
            else:
                identities = {_model_key(row.get(field) or "") for field in ("cif_file", "input")}
                identities.discard("")
                if len(identities) != 1:
                    raise ValueError(f"Missing or conflicting model identity in score CSV row {line}")
                identity = identities.pop()
            if identity in index:
                raise ValueError(f"Duplicate score mapping for {identity}")
            index[identity] = {"source_row": line, "raw_source_row": row,
                               "native_boltz_affinity": native, "available_fields": available}
    return index


def _resolve_source_key(value: str, index: dict[str, dict[str, object]]) -> str | None:
    """Resolve exact model names or a verified model-0 native Boltz wrapper only."""
    identity = _model_key(value)
    if identity in index:
        return identity
    wrapped = re.fullmatch(r"(?P<job>[^/\\]+)__(?P=job)_model_0", identity)
    if wrapped and index.get(wrapped["job"], {}).get("native_boltz_affinity"):
        return wrapped["job"]
    return None


def _extract(source: dict[str, object], identity: str, source_csv: Path) -> dict[str, object]:
    raw = source["raw_source_row"]
    values = {column["key"]: _score(raw.get(column["source_column"]), column["source_column"], identity)
              for column in score_columns()}
    numeric: dict[str, float | None] = {}
    for field, value in raw.items():
        if field in BOLTZ2_FIELDS + GNINA_FIELDS or field.startswith(("chain_ptm_", "pair_iptm_")):
            numeric[field] = _score(value, field, identity)
        elif value and value.strip():
            try:
                number = float(value)
            except ValueError:
                continue
            if not math.isfinite(number):
                raise ValueError(f"Nonfinite source value {field}={value!r} for {identity}")
            numeric[field] = number
    return {"source_csv": str(source_csv), "source_model": identity,
            "source_row": source["source_row"], "values": values,
            "numeric_source_fields": numeric, "raw_source_row": raw}


def _write_detail(result_dir: Path, detail: dict[str, object]) -> None:
    if not result_dir.is_dir():
        return
    write_json(result_dir / "model_scores.json", detail)
    json_path = result_dir / "automatic_summary.json"
    if json_path.is_file():
        summary = json.loads(json_path.read_text(encoding="utf-8"))
        summary["model_scores"] = detail
        write_json(json_path, summary)
    md_path = result_dir / "automatic_summary.md"
    if not md_path.is_file():
        return
    lines = ["## Boltz-2 and GNINA scores", "",
             f"Source: `{detail['source_csv']}`; exact model `{detail['source_model']}`; "
             f"CSV row {detail['source_row']} (header is row 1).", "",
             "Original source values and labels are retained. Boltz complex_plddt and complex_iplddt "
             "use the source 0–1 scale; no scores were recomputed or rescaled.", "",
             "| Section | Source score | Value |", "|---|---|---:|"]
    for column in score_columns():
        score = detail["values"][column["key"]]
        lines.append(f"| {column['section']} | {column['source_column']} | {score if score is not None else 'NA'} |")
    begin, end = "<!-- model-scores:begin -->", "<!-- model-scores:end -->"
    block = begin + "\n" + "\n".join(lines) + "\n" + end
    text = md_path.read_text(encoding="utf-8")
    if begin in text and end in text:
        prefix, remainder = text.split(begin, 1)
        _, suffix = remainder.split(end, 1)
        text = prefix.rstrip() + "\n\n" + block + suffix
    else:
        text = text.rstrip() + "\n\n" + block + "\n"
    md_path.write_text(text, encoding="utf-8")


def _safe_csv(path: Path, rows: list[dict[str, object]]) -> Path:
    try:
        write_csv(path, rows)
        return path
    except PermissionError:
        alternative = path.with_stem(path.stem + "_with_model_scores")
        write_csv(alternative, rows)
        return alternative


def add_model_scores(batch: Path, score_csv: Path) -> Path:
    """Append two score sections, validating every model join before any writes."""
    batch, score_csv = Path(batch).resolve(), Path(score_csv).resolve()
    payload_path = batch / "batch_summary_grouped_fingerprints.json"
    payload = json.loads(payload_path.read_text(encoding="utf-8"))
    manifest_path = batch / "batch_manifest.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8")) if manifest_path.is_file() else {}
    index = _source_index(score_csv)
    matched: dict[str, dict[str, object]] = {}
    for row in payload["rows"]:
        candidate = row["candidate"]
        identity = _resolve_source_key(candidate, index)
        if candidate in matched:
            raise ValueError(f"Duplicate candidate in grouped report: {candidate}")
        if identity is None:
            raise ValueError(f"No exact score mapping for {candidate}")
        matched[candidate] = _extract(index[identity], identity, score_csv)

    numeric_path = batch / "batch_interaction_fingerprint_numeric.csv"
    numeric = _read_csv(numeric_path) if numeric_path.is_file() else None
    if numeric is not None:
        for row in numeric:
            if row["candidate"] not in matched:
                raise ValueError(f"Numeric report candidate is absent from grouped report: {row['candidate']}")
            row.update(matched[row["candidate"]]["values"])

    parent_values = {column["key"]: None for column in score_columns()}
    parent_key = _resolve_source_key(manifest.get("parent") or payload.get("parent", {}).get("candidate", ""), index)
    if parent_key is not None:
        parent_values = _extract(index[parent_key], parent_key, score_csv)["values"]
    available = next(iter(index.values()), {}).get("available_fields", [])
    metadata = {"source_csv": str(score_csv), "matched_candidates": len(matched),
                "parent_scores_available": parent_key is not None,
                "available_fields": available,
                "available_sections": [section for section in SECTIONS if any(
                    column["section"] == section and column["source_column"] in available
                    for column in score_columns())],
                "mapping": {candidate: {"source_model": detail["source_model"], "source_row": detail["source_row"]}
                            for candidate, detail in matched.items()}}
    payload["columns"] = [c for c in payload["columns"] if c["section"] not in SECTIONS] + score_columns()
    for row in payload["rows"]:
        row.update(matched[row["candidate"]]["values"])
    payload.setdefault("parent", {}).update(parent_values)
    payload["model_scores"] = metadata
    # Optional detailed reports are updated only where a result folder is known.
    summary_path = batch / "batch_summary.csv"
    if summary_path.is_file():
        for base in _read_csv(summary_path):
            if base.get("candidate") in matched and base.get("result_folder"):
                _write_detail(_windows_path(base["result_folder"]), matched[base["candidate"]])
    write_json(payload_path, payload)
    csv_path = _safe_csv(batch / "batch_summary_grouped_fingerprints.csv", payload["rows"])
    if numeric is not None:
        numeric_path = _safe_csv(numeric_path, numeric)
    manifest["model_scores"] = {**metadata, "summary_csv": str(csv_path),
                                "numeric_csv": str(numeric_path) if numeric is not None else None}
    manifest.setdefault("interaction_fingerprints", {})["grouped_summary_csv"] = str(csv_path)
    if numeric is not None:
        manifest["interaction_fingerprints"]["numeric_flat_csv"] = str(numeric_path)
    write_json(manifest_path, manifest)
    return csv_path


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--batch-results", required=True)
    parser.add_argument("--scores-csv", required=True)
    args = parser.parse_args(argv)
    result = add_model_scores(Path(args.batch_results), Path(args.scores_csv))
    print(f"Boltz-2 and GNINA scores added: {result}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
