from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Iterable

import numpy as np


@dataclass(frozen=True, order=True)
class ResidueKey:
    chain: str
    number: int
    insertion: str = ""

    def label(self) -> str:
        return f"{self.chain}:{self.number}{self.insertion}"


@dataclass
class Atom:
    serial: int
    name: str
    altloc: str
    residue_name: str
    residue_key: ResidueKey
    coord: np.ndarray
    element: str
    occupancy: float
    line_index: int


@dataclass
class Residue:
    key: ResidueKey
    name: str
    atoms: dict[str, Atom] = field(default_factory=dict)
    hetero: bool = False

    def atom(self, name: str) -> Atom | None:
        return self.atoms.get(name)


@dataclass
class Structure:
    path: Path
    lines: list[str]
    residues: list[Residue]

    @property
    def by_key(self) -> dict[ResidueKey, Residue]:
        return {r.key: r for r in self.residues}

    def chain_residues(self, chain: str, protein_only: bool = True) -> list[Residue]:
        values = [r for r in self.residues if r.key.chain == chain]
        return [r for r in values if not r.hetero] if protein_only else values

    def pose_residue(self, one_based_index: int) -> Residue:
        if one_based_index < 1 or one_based_index > len(self.residues):
            raise IndexError(f"Pose residue {one_based_index} is outside 1..{len(self.residues)}")
        return self.residues[one_based_index - 1]


def residue_key_sort(key: ResidueKey) -> tuple[str, int, str]:
    return key.chain, key.number, key.insertion


def canonical_pair(a: ResidueKey, b: ResidueKey) -> tuple[ResidueKey, ResidueKey]:
    return (a, b) if residue_key_sort(a) <= residue_key_sort(b) else (b, a)


def pair_id(a: ResidueKey, b: ResidueKey) -> str:
    left, right = canonical_pair(a, b)
    return f"{left.label()}--{right.label()}"


def coords(atoms: Iterable[Atom]) -> np.ndarray:
    return np.asarray([atom.coord for atom in atoms], dtype=float)
