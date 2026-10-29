from __future__ import annotations

import math

import numpy as np

from .models import Residue


CHI_ATOMS: dict[str, tuple[tuple[str, str, str, str], ...]] = {
    "ARG": (("N", "CA", "CB", "CG"), ("CA", "CB", "CG", "CD"), ("CB", "CG", "CD", "NE"), ("CG", "CD", "NE", "CZ")),
    "ASN": (("N", "CA", "CB", "CG"), ("CA", "CB", "CG", "OD1")),
    "ASP": (("N", "CA", "CB", "CG"), ("CA", "CB", "CG", "OD1")),
    "CYS": (("N", "CA", "CB", "SG"),),
    "GLN": (("N", "CA", "CB", "CG"), ("CA", "CB", "CG", "CD"), ("CB", "CG", "CD", "OE1")),
    "GLU": (("N", "CA", "CB", "CG"), ("CA", "CB", "CG", "CD"), ("CB", "CG", "CD", "OE1")),
    "HIS": (("N", "CA", "CB", "CG"), ("CA", "CB", "CG", "ND1")),
    "ILE": (("N", "CA", "CB", "CG1"), ("CA", "CB", "CG1", "CD1")),
    "LEU": (("N", "CA", "CB", "CG"), ("CA", "CB", "CG", "CD1")),
    "LYS": (("N", "CA", "CB", "CG"), ("CA", "CB", "CG", "CD"), ("CB", "CG", "CD", "CE"), ("CG", "CD", "CE", "NZ")),
    "MET": (("N", "CA", "CB", "CG"), ("CA", "CB", "CG", "SD"), ("CB", "CG", "SD", "CE")),
    "PHE": (("N", "CA", "CB", "CG"), ("CA", "CB", "CG", "CD1")),
    "PRO": (("N", "CA", "CB", "CG"), ("CA", "CB", "CG", "CD")),
    "SER": (("N", "CA", "CB", "OG"),),
    "THR": (("N", "CA", "CB", "OG1"),),
    "TRP": (("N", "CA", "CB", "CG"), ("CA", "CB", "CG", "CD1")),
    "TYR": (("N", "CA", "CB", "CG"), ("CA", "CB", "CG", "CD1")),
    "VAL": (("N", "CA", "CB", "CG1"),),
}


def dihedral_degrees(a: np.ndarray, b: np.ndarray, c: np.ndarray, d: np.ndarray) -> float:
    b0 = -(b - a)
    b1 = c - b
    b2 = d - c
    b1 = b1 / np.linalg.norm(b1)
    v = b0 - np.dot(b0, b1) * b1
    w = b2 - np.dot(b2, b1) * b1
    angle = math.degrees(math.atan2(np.dot(np.cross(b1, v), w), np.dot(v, w)))
    return ((angle + 180.0) % 360.0) - 180.0


def chi_angles(residue: Residue) -> list[float]:
    output: list[float] = []
    for names in CHI_ATOMS.get(residue.name, ()):
        atoms = [residue.atom(name) for name in names]
        if any(atom is None for atom in atoms):
            break
        output.append(dihedral_degrees(*(atom.coord for atom in atoms if atom is not None)))
    return output


def angle_bin(angle: float) -> str:
    if not math.isfinite(angle):
        raise ValueError("Cannot assign a torsion bin to a non-finite angle")
    angle = ((angle + 180.0) % 360.0) - 180.0
    if 0.0 <= angle < 120.0:
        return "p"  # gauche+
    if -120.0 <= angle < 0.0:
        return "m"  # gauche-
    return "t"      # trans


def rotamer_bin(residue: Residue) -> tuple[str | None, list[float]]:
    angles = chi_angles(residue)
    if not angles:
        return None, []
    return "".join(angle_bin(angle) for angle in angles), angles


def compare_rotamers(parent: Residue, mutant: Residue) -> dict[str, object]:
    p_bin, p_angles = rotamer_bin(parent)
    m_bin, m_angles = rotamer_bin(mutant)
    comparable = min(len(p_angles), len(m_angles))
    same_chi1 = (p_bin[0] == m_bin[0]) if p_bin and m_bin else None
    same_common = (p_bin[:comparable] == m_bin[:comparable]) if comparable else None
    p_complete = len(p_angles) == len(CHI_ATOMS.get(parent.name, ()))
    m_complete = len(m_angles) == len(CHI_ATOMS.get(mutant.name, ()))
    # Different amino acids do not share a residue-specific Dunbrack rotamer
    # identity merely because their coarse torsion strings happen to match.
    same_full = (p_bin == m_bin) if p_bin and m_bin and parent.name == mutant.name and p_complete and m_complete else None
    return {
        "parent_rotamer_bin": p_bin or "NA",
        "mutant_rotamer_bin": m_bin or "NA",
        "parent_chi_degrees": ";".join(f"{x:.1f}" for x in p_angles),
        "mutant_chi_degrees": ";".join(f"{x:.1f}" for x in m_angles),
        "n_comparable_chi": comparable,
        "same_chi1_bin": same_chi1,
        "same_common_chi_bins": same_common,
        "same_full_rotamer_bin": same_full,
        "rotamer_method": "canonical_chi_coarse_torsion_bins_not_dunbrack",
        "parent_chi_complete": p_complete,
        "mutant_chi_complete": m_complete,
        "rotamer_comparison_note": (
            "Coarse chi-angle bins only; not a Dunbrack rotamer assignment. "
            "Cross-amino-acid comparisons describe torsions, not identical rotamers."
        ),
    }
