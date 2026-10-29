from __future__ import annotations

from pathlib import Path

import numpy as np

from .models import Atom, Residue, ResidueKey, Structure


THREE_TO_ONE = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
    "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
    "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
    "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
    "MSE": "M",
}


def _element(atom_name: str, element_field: str) -> str:
    value = element_field.strip().upper()
    if value:
        return value
    letters = "".join(c for c in atom_name if c.isalpha()).upper()
    return (letters[:2] if letters[:2] in {"CL", "BR"} else letters[:1]) or "?"


def read_pdb(path: str | Path) -> Structure:
    path = Path(path)
    lines = path.read_text(encoding="utf-8", errors="replace").splitlines(keepends=True)
    residues: list[Residue] = []
    lookup: dict[ResidueKey, Residue] = {}
    seen_model = False
    for line_index, line in enumerate(lines):
        record = line[:6].strip()
        if record == "MODEL":
            if seen_model:
                break
            seen_model = True
            continue
        if record == "ENDMDL" and seen_model:
            break
        if record not in {"ATOM", "HETATM"}:
            continue
        padded = line.rstrip("\r\n").ljust(80)
        altloc = padded[16].strip()
        if altloc not in {"", "A"}:
            continue
        chain = padded[21].strip() or "_"
        try:
            number = int(padded[22:26])
            serial = int(padded[6:11])
            xyz = np.array([float(padded[30:38]), float(padded[38:46]), float(padded[46:54])])
        except ValueError as exc:
            raise ValueError(f"Cannot parse PDB coordinate line {line_index + 1} in {path}") from exc
        key = ResidueKey(chain, number, padded[26].strip())
        resname = padded[17:20].strip().upper()
        residue = lookup.get(key)
        if residue is None:
            residue = Residue(key=key, name=resname, hetero=(record == "HETATM" and resname != "MSE"))
            lookup[key] = residue
            residues.append(residue)
        atom_name = padded[12:16].strip()
        occupancy_text = padded[54:60].strip()
        atom = Atom(
            serial=serial,
            name=atom_name,
            altloc=altloc,
            residue_name=resname,
            residue_key=key,
            coord=xyz,
            element=_element(atom_name, padded[76:78]),
            occupancy=float(occupancy_text) if occupancy_text else 1.0,
            line_index=line_index,
        )
        existing = residue.atoms.get(atom_name)
        if existing is None or atom.occupancy > existing.occupancy:
            residue.atoms[atom_name] = atom
    if not residues:
        raise ValueError(f"No ATOM/HETATM records found in {path}")
    return Structure(path=path, lines=lines, residues=residues)


def residue_sequence(residues: list[Residue]) -> str:
    unknown = [(r.key.label(), r.name) for r in residues if r.name not in THREE_TO_ONE]
    if unknown:
        details = ", ".join(f"{key}={name}" for key, name in unknown[:10])
        raise ValueError(f"Noncanonical polymer residues cannot be sequenced: {details}")
    return "".join(THREE_TO_ONE[r.name] for r in residues)


def write_transformed_pdb(structure: Structure, rotation: np.ndarray, translation: np.ndarray, path: str | Path) -> None:
    lines = list(structure.lines)
    for residue in structure.residues:
        for atom in residue.atoms.values():
            xyz = atom.coord @ rotation + translation
            raw = lines[atom.line_index].rstrip("\r\n").ljust(80)
            newline = f"{raw[:30]}{xyz[0]:8.3f}{xyz[1]:8.3f}{xyz[2]:8.3f}{raw[54:]}"
            lines[atom.line_index] = newline.rstrip() + "\n"
    Path(path).write_text("".join(lines), encoding="utf-8")
