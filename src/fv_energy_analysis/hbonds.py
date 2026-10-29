from __future__ import annotations

import csv
from collections import defaultdict
from pathlib import Path

from .models import ResidueKey, pair_id
from .report import write_csv


def _pdb_key(pdb_info, pose_index: int) -> ResidueKey:
    chain = str(pdb_info.chain(pose_index)).strip() or "_"
    number = int(pdb_info.number(pose_index))
    insertion = str(pdb_info.icode(pose_index)).strip()
    return ResidueKey(chain, number, insertion)


def extract_pyrosetta_hbonds(
    pdb: str | Path,
    output: str | Path,
    params: list[str] | None = None,
    weights: str = "ref2015",
    legacy_pre_talaris: bool | None = None,
) -> None:
    """Write one row per scored Rosetta hydrogen bond.

    PyRosetta is imported only here, keeping it an optional licensed dependency.
    """
    from .pyrosetta_score import _initialize, _scorefunction
    from .pyrosetta_score import resolve_legacy_pre_talaris
    pyrosetta = _initialize(params, legacy=resolve_legacy_pre_talaris(weights, legacy_pre_talaris))
    pose = pyrosetta.pose_from_pdb(str(pdb))
    scorefxn = _scorefunction(pyrosetta, weights)
    write_pose_hbonds(pose, output, scorefxn=scorefxn)


def write_pose_hbonds(pose, output: str | Path, scorefxn=None) -> None:
    """Write Rosetta H-bonds from an already loaded pose.

    This is used by persistent batch workers so PyRosetta and the pose are not
    reinitialized merely to extract H-bonds after scoring.
    """
    try:
        from pyrosetta.rosetta.core.scoring.hbonds import HBondSet, fill_hbond_set, hb_eval_type_weight  # type: ignore
    except ImportError as exc:  # pragma: no cover - optional licensed dependency
        raise RuntimeError("PyRosetta is required for H-bond extraction") from exc
    if scorefxn is not None:
        scorefxn(pose)
    hbond_set = HBondSet(scorefxn.energy_method_options().hbond_options()) if scorefxn is not None else HBondSet()
    fill_hbond_set(pose, False, hbond_set)
    pdb_info = pose.pdb_info()
    rows: list[dict[str, object]] = []
    for index in range(1, hbond_set.nhbonds() + 1):
        hb = hbond_set.hbond(index)
        if not hbond_set.allow_hbond(index):
            continue
        donor_i, acceptor_i = int(hb.don_res()), int(hb.acc_res())
        donor_res, acceptor_res = pose.residue(donor_i), pose.residue(acceptor_i)
        donor_h = int(hb.don_hatm())
        donor_heavy = int(donor_res.atom_base(donor_h))
        acceptor_atom = int(hb.acc_atm())
        donor_key, acceptor_key = _pdb_key(pdb_info, donor_i), _pdb_key(pdb_info, acceptor_i)
        rows.append({
            "hbond_index": index,
            "donor_residue": donor_key.label(),
            "donor_residue_name": donor_res.name3(),
            "donor_atom": donor_res.atom_name(donor_heavy).strip(),
            "donor_hydrogen": donor_res.atom_name(donor_h).strip(),
            "acceptor_residue": acceptor_key.label(),
            "acceptor_residue_name": acceptor_res.name3(),
            "acceptor_atom": acceptor_res.atom_name(acceptor_atom).strip(),
            "hbond_energy": float(hb.energy()),  # backward-compatible, explicitly RAW
            "hbond_raw_energy": float(hb.energy()),
            "hbond_environment_weight": float(hb.weight()),
            "hbond_scorefunction_weight": float(hb_eval_type_weight(hb.eval_type(), scorefxn.weights(), donor_i == acceptor_i)) if scorefxn is not None else "",
            "hbond_weighted_energy_REU": (float(hb.energy()) * float(hb.weight()) * float(hb_eval_type_weight(hb.eval_type(), scorefxn.weights(), donor_i == acceptor_i))) if scorefxn is not None else "",
            "hbond_energy_note": "hbond_energy is raw; use hbond_weighted_energy_REU for score-weighted energy",
            "hbond_evaluation_type": str(hb.eval_type()),
            "pair_id": pair_id(donor_key, acceptor_key),
        })
    write_csv(output, rows)


def aggregate_hbond_csv(path: str | Path) -> list[dict[str, object]]:
    grouped: dict[str, list[float]] = defaultdict(list)
    weighted: dict[str, list[float]] = defaultdict(list)
    identity: dict[str, tuple[str, str]] = {}
    with Path(path).open(newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
            key = row.get("pair_id", "")
            if not key:
                donor = row["donor_residue"]
                acceptor = row["acceptor_residue"]
                key = "--".join(sorted((donor, acceptor)))
            grouped[key].append(float(row["hbond_energy"]))
            if row.get('hbond_weighted_energy_REU') not in (None,''):
                weighted[key].append(float(row['hbond_weighted_energy_REU']))
            identity[key] = tuple(key.split("--", 1))  # type: ignore[assignment]
    rows = []
    for key, energies in grouped.items():
        residue1, residue2 = identity[key]
        rows.append({
            "pair_id": key,
            "residue1": residue1,
            "residue2": residue2,
            "count_rosetta_hbond": len(energies),
            "sum_individual_hbond_energy": sum(energies),
            "mean_individual_hbond_energy": sum(energies) / len(energies),
            "strongest_individual_hbond_energy": min(energies),
            "individual_hbond_energy_units": "legacy individual_hbond_energy fields are raw; weighted energies are separate",
            "sum_individual_hbond_weighted_energy_REU": sum(weighted[key]) if len(weighted[key]) == len(energies) else '',
            "mean_individual_hbond_weighted_energy_REU": sum(weighted[key])/len(energies) if len(weighted[key]) == len(energies) else '',
            "strongest_individual_hbond_weighted_energy_REU": min(weighted[key]) if len(weighted[key]) == len(energies) else '',
        })
    return rows
