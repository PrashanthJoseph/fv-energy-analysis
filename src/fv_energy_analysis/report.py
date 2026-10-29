from __future__ import annotations

import csv
import json
from pathlib import Path


IDENTITY_COLUMNS = {"pair_id", "residue1", "residue1_name", "residue2", "residue2_name"}


def write_csv(path: str | Path, rows: list[dict[str, object]]) -> None:
    path = Path(path)
    keys: list[str] = []
    seen: set[str] = set()
    for row in rows:
        for key in row:
            if key not in seen:
                keys.append(key)
                seen.add(key)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=keys or ["empty"])
        writer.writeheader()
        writer.writerows(rows)


def write_json(path: str | Path, value) -> None:
    Path(path).write_text(json.dumps(value, indent=2, allow_nan=False, default=str) + "\n", encoding="utf-8")


def outer_compare(parent: list[dict[str, object]], mutant: list[dict[str, object]]) -> list[dict[str, object]]:
    p_by = {str(row["pair_id"]): row for row in parent}
    m_by = {str(row["pair_id"]): row for row in mutant}
    rows: list[dict[str, object]] = []
    for key in sorted(set(p_by) | set(m_by)):
        p, m = p_by.get(key, {}), m_by.get(key, {})
        out: dict[str, object] = {"pair_id": key}
        for col in ("residue1", "residue1_name", "residue2", "residue2_name"):
            out[f"parent_{col}"] = p.get(col, "")
            out[f"mutant_{col}"] = m.get(col, "")
        columns = (set(p) | set(m)) - IDENTITY_COLUMNS
        for col in sorted(columns):
            pv, mv = p.get(col, 0), m.get(col, 0)
            out[f"parent_{col}"] = pv
            out[f"mutant_{col}"] = mv
            if isinstance(pv, (int, float)) and isinstance(mv, (int, float)):
                out[f"delta_{col}"] = float(mv) - float(pv)
        rows.append(out)
    return rows


def merge_by_pair(*tables: list[dict[str, object]]) -> list[dict[str, object]]:
    merged: dict[str, dict[str, object]] = {}
    for table in tables:
        for row in table:
            key = str(row["pair_id"])
            merged.setdefault(key, {"pair_id": key}).update(row)
    return [merged[key] for key in sorted(merged)]
