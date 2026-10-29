from __future__ import annotations

from collections import defaultdict
from dataclasses import dataclass

import numpy as np

from .models import Residue, ResidueKey, Structure
from .structure import residue_sequence


BACKBONE = ("N", "CA", "C", "O")


@dataclass(frozen=True)
class Transform:
    rotation: np.ndarray
    translation: np.ndarray

    def apply(self, xyz: np.ndarray) -> np.ndarray:
        return xyz @ self.rotation + self.translation


def kabsch(mobile: np.ndarray, reference: np.ndarray) -> Transform:
    if mobile.shape != reference.shape or mobile.ndim != 2 or mobile.shape[1] != 3:
        raise ValueError(f"Coordinate arrays must both be N x 3; got {mobile.shape} and {reference.shape}")
    if len(mobile) < 3:
        raise ValueError("At least three atoms are required for a rigid-body fit")
    if not np.isfinite(mobile).all() or not np.isfinite(reference).all():
        raise ValueError("Rigid-body fitting requires finite coordinates")
    mobile_center = mobile.mean(axis=0)
    reference_center = reference.mean(axis=0)
    covariance = (mobile - mobile_center).T @ (reference - reference_center)
    left, _, right_t = np.linalg.svd(covariance)
    correction = np.eye(3)
    correction[-1, -1] = np.sign(np.linalg.det(left @ right_t))
    rotation = left @ correction @ right_t
    translation = reference_center - mobile_center @ rotation
    return Transform(rotation, translation)


def rmsd(mobile: np.ndarray, reference: np.ndarray) -> float:
    if mobile.shape != reference.shape or len(mobile) == 0:
        return float("nan")
    return float(np.sqrt(np.mean(np.sum((mobile - reference) ** 2, axis=1))))


def paired_atoms(
    parent: Structure,
    mutant: Structure,
    parent_keys: list[ResidueKey],
    mutant_keys: list[ResidueKey],
    atom_names: tuple[str, ...] = BACKBONE,
) -> tuple[np.ndarray, np.ndarray, int]:
    p_by = parent.by_key
    m_by = mutant.by_key
    p_coords: list[np.ndarray] = []
    m_coords: list[np.ndarray] = []
    residues_used = 0
    for p_key, m_key in zip(parent_keys, mutant_keys):
        p_res, m_res = p_by[p_key], m_by[m_key]
        before = len(p_coords)
        for name in atom_names:
            p_atom, m_atom = p_res.atom(name), m_res.atom(name)
            if p_atom is not None and m_atom is not None:
                p_coords.append(p_atom.coord)
                m_coords.append(m_atom.coord)
        residues_used += int(len(p_coords) > before)
    return np.asarray(p_coords).reshape(-1, 3), np.asarray(m_coords).reshape(-1, 3), residues_used


def hcdr3_segments(chain_pairs, parent_regions, heavy_chain_ids=("H",), tip_sequence="YYGSS"):
    """Map the unique parent tip motif to corresponding mutant positions.

    Never search the mutant sequence: a substitution in the tip must not change
    which positions are compared. Anchors are the complement *within HCDR3*.
    """
    if not tip_sequence or any(aa not in "ACDEFGHIKLMNPQRSTVWY" for aa in tip_sequence):
        raise ValueError("HCDR3 tip sequence must contain canonical amino-acid letters")
    output = {}
    for chain in heavy_chain_ids:
        if chain not in chain_pairs:
            continue
        pairs = [(p, m) for p, m in zip(*chain_pairs[chain]) if parent_regions.get(p.key) == "CDR3"]
        sequence = residue_sequence([p for p, _ in pairs])
        starts = [i for i in range(len(sequence) - len(tip_sequence) + 1)
                  if sequence[i:i + len(tip_sequence)] == tip_sequence]
        status = "complete" if len(starts) == 1 else "tip_motif_not_found" if not starts else "tip_motif_ambiguous"
        for segment in ("ANCHORS", "TIP"):
            selected = []
            if len(starts) == 1:
                tip_positions = set(range(starts[0], starts[0] + len(tip_sequence)))
                selected = [pair for i, pair in enumerate(pairs) if (i in tip_positions) == (segment == "TIP")]
            output[f"{chain}_CDR3_{segment}"] = {
                "parent_keys": [p.key for p, _ in selected],
                "mutant_keys": [m.key for _, m in selected],
                "segment_status": status if selected or status != "complete" else "empty_segment",
                "tip_reference_sequence": tip_sequence,
                "parent_cdr3_sequence": sequence,
            }
    return output


def structural_rmsds(
    parent: Structure,
    mutant: Structure,
    chain_pairs: dict[str, tuple[list[Residue], list[Residue]]],
    parent_regions: dict[ResidueKey, str],
    fit_region: str = "fv",
    atom_names: tuple[str, ...] = BACKBONE,
    heavy_chain_ids: tuple[str, ...] = ("H",),
    hcdr3_tip_sequence: str = "YYGSS",
) -> tuple[Transform, list[dict[str, object]]]:
    groups: dict[str, tuple[list[ResidueKey], list[ResidueKey]]] = defaultdict(lambda: ([], []))
    all_parent: list[ResidueKey] = []
    all_mutant: list[ResidueKey] = []
    fr_parent: list[ResidueKey] = []
    fr_mutant: list[ResidueKey] = []
    for chain_id, (p_residues, m_residues) in chain_pairs.items():
        for p_res, m_res in zip(p_residues, m_residues):
            region = parent_regions.get(p_res.key, "UNNUMBERED")
            all_parent.append(p_res.key)
            all_mutant.append(m_res.key)
            if region.startswith("FR"):
                fr_parent.append(p_res.key)
                fr_mutant.append(m_res.key)
            if region.startswith("CDR") or region.startswith("FR"):
                label = f"{chain_id}_{region}"
                groups[label][0].append(p_res.key)
                groups[label][1].append(m_res.key)
    if fit_region == "framework":
        fit_parent, fit_mutant = fr_parent, fr_mutant
    else:
        fit_parent, fit_mutant = all_parent, all_mutant
    p_fit, m_fit, fit_residues = paired_atoms(parent, mutant, fit_parent, fit_mutant, atom_names)
    transform = kabsch(m_fit, p_fit)

    report_groups = {
        "GLOBAL_FV": (all_parent, all_mutant),
        "FRAMEWORK": (fr_parent, fr_mutant),
        **dict(groups),
    }
    segments = hcdr3_segments(chain_pairs, parent_regions, heavy_chain_ids, hcdr3_tip_sequence)
    report_groups.update({label: (item["parent_keys"], item["mutant_keys"]) for label, item in segments.items()})
    rows: list[dict[str, object]] = []
    for label, (p_keys, m_keys) in report_groups.items():
        p_xyz, m_xyz, n_residues = paired_atoms(parent, mutant, p_keys, m_keys, atom_names)
        after_global = rmsd(transform.apply(m_xyz), p_xyz)
        local_value = float("nan")
        if len(p_xyz) >= 3:
            local_transform = kabsch(m_xyz, p_xyz)
            local_value = rmsd(local_transform.apply(m_xyz), p_xyz)
        segment = segments.get(label, {})
        rows.append({
            "region": label,
            "fit_basis": fit_region,
            "atom_names": ",".join(atom_names),
            "n_residues": n_residues,
            "n_atoms": len(p_xyz),
            "rmsd_after_global_fit_A": after_global if np.isfinite(after_global) else "",
            "local_fit_rmsd_A": local_value if np.isfinite(local_value) else "",
            "parent_residues": ";".join(key.label() for key in p_keys),
            "mutant_residues": ";".join(key.label() for key in m_keys),
            "missing_atom_pairs": len(p_keys) * len(atom_names) - len(p_xyz),
            "segment_status": segment.get("segment_status", "complete" if len(p_xyz) else "no_coordinates"),
            "tip_reference_sequence": segment.get("tip_reference_sequence", ""),
            "parent_cdr3_sequence": segment.get("parent_cdr3_sequence", ""),
            "fit_n_residues": fit_residues,
            "fit_n_atoms": len(p_fit),
        })
    return transform, rows
