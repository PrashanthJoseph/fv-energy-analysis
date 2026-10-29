from __future__ import annotations

import subprocess
from pathlib import Path

from .models import Structure, canonical_pair, pair_id


NON_ENERGY_COLUMNS = {
    "pose_id", "resi1", "resi2", "pdbid1", "pdbid2", "restype1", "restype2", "description",
}


def _number(value: str):
    try:
        return int(value)
    except ValueError:
        try:
            return float(value)
        except ValueError:
            return value


def read_residue_energy_breakdown(path: str | Path, structure: Structure) -> list[dict[str, object]]:
    header: list[str] | None = None
    rows: list[dict[str, object]] = []
    for raw in Path(path).read_text(encoding="utf-8", errors="replace").splitlines():
        if not raw.startswith("SCORE:"):
            continue
        fields = raw.split()[1:]
        if "resi1" in fields and "restype1" in fields:
            header = fields
            continue
        if header is None or len(fields) < len(header):
            continue
        record = {key: _number(value) for key, value in zip(header, fields)}
        if str(record.get("resi2", "--")) == "--" or str(record.get("restype2", "")) == "onebody":
            continue
        try:
            r1 = structure.pose_residue(int(record["resi1"]))
            r2 = structure.pose_residue(int(record["resi2"]))
        except (ValueError, KeyError, IndexError) as exc:
            raise ValueError(f"Cannot map Rosetta pose numbering in {path}: {record}") from exc
        for suffix, residue in (("1", r1), ("2", r2)):
            pdb_id = f"{residue.key.number}{residue.key.insertion}{residue.key.chain}"
            native_pdb_id = f"{residue.key.number}{residue.key.chain}" + (f":{residue.key.insertion}" if residue.key.insertion else "")
            if f"pdbid{suffix}" in record and str(record[f"pdbid{suffix}"]) not in {pdb_id, native_pdb_id}:
                raise ValueError(f"Rosetta/PDB residue mapping mismatch in {path}: expected {pdb_id}, got {record[f'pdbid{suffix}']}")
            if f"restype{suffix}" in record and str(record[f"restype{suffix}"]) != residue.name:
                raise ValueError(f"Rosetta/PDB residue type mismatch in {path}: {residue.key.label()}")
        left, right = canonical_pair(r1.key, r2.key)
        if left != r1.key:
            r1, r2 = r2, r1
        out: dict[str, object] = {
            "pair_id": pair_id(left, right),
            "residue1": left.label(), "residue1_name": r1.name,
            "residue2": right.label(), "residue2_name": r2.name,
        }
        for key, value in record.items():
            if key not in NON_ENERGY_COLUMNS and isinstance(value, (int, float)):
                out[key] = float(value)
        out["vdw_net"] = float(out.get("fa_atr", 0.0)) + float(out.get("fa_rep", 0.0))
        out["hbond_total"] = sum(float(v) for k, v in out.items() if k.startswith("hbond_") and k != "hbond_total")
        total = float(out.get("total", 0.0))
        out["nonbonded_total"] = total - sum(float(out.get(term, 0.0)) for term in ("rama_prepro", "dslf_fa13", "cart_bonded", "pro_close"))
        out["pair_energy_note"] = "total includes coupled conformational/covalent terms; nonbonded_total excludes rama_prepro, dslf_fa13, cart_bonded, pro_close"
        out["pair_energy_strength"] = energy_strength(total)
        rows.append(out)
    if header is None:
        raise ValueError(f"No residue_energy_breakdown header found in {path}")
    return rows


def read_scorefile_metrics(path: str | Path) -> dict[str, float]:
    """Read the last numeric SCORE row from a standard Rosetta scorefile."""
    header: list[str] | None = None
    last: dict[str, float] | None = None
    for raw in Path(path).read_text(encoding="utf-8", errors="replace").splitlines():
        if not raw.startswith("SCORE:"):
            continue
        fields = raw.split()[1:]
        if not fields:
            continue
        if any(name in fields for name in ("description", "total_score", "score")) and not _looks_numeric(fields[0]):
            header = fields
            continue
        if header is None:
            continue
        record: dict[str, float] = {}
        for key, value in zip(header, fields):
            if key == "description":
                continue
            try:
                record[key] = float(value)
            except ValueError:
                continue
        if record:
            last = record
    if last is None:
        raise ValueError(f"No numeric Rosetta SCORE row found in {path}")
    return last


def _looks_numeric(value: str) -> bool:
    try:
        float(value)
        return True
    except ValueError:
        return False


def energy_strength(value: float, strong: float = -2.0, moderate: float = -1.0) -> str:
    """Heuristic label; the numeric REU value remains authoritative."""
    if value <= strong:
        return "strong_favorable"
    if value <= moderate:
        return "moderate_favorable"
    if value < 0:
        return "weak_favorable"
    if value < 1.0:
        return "neutral_or_weak_unfavorable"
    return "unfavorable"


def run_residue_energy_breakdown(
    executable: str,
    pdb: str | Path,
    output: str | Path,
    params: list[str] | None = None,
    weights: str = "ref2015",
    database: str | None = None,
) -> None:
    command = [executable, "-in:file:s", str(pdb), "-out:file:silent", str(output), "-score:weights", weights]
    if database:
        command += ["-database", database]
    for param in params or []:
        command += ["-extra_res_fa", param]
    completed = subprocess.run(command, check=False, text=True, capture_output=True)
    if completed.returncode:
        raise RuntimeError(
            f"Rosetta failed with exit code {completed.returncode}\nSTDOUT:\n{completed.stdout}\nSTDERR:\n{completed.stderr}"
        )
