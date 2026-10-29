"""Input and ligand-topology validation helpers."""

from __future__ import annotations

import math
from pathlib import Path


def validate_ligand_against_params(
    pdb: str | Path, params: str | Path, chain: str, template: str | Path | None = None,
) -> dict[str, object]:
    pdb, params = Path(pdb), Path(params)
    atoms: dict[str, tuple[str, tuple[float, float, float]]] = {}
    duplicate_names: list[str] = []
    for line in pdb.read_text(encoding="utf-8").splitlines():
        if not line.startswith("HETATM") or line[21:22].strip() != chain:
            continue
        name = line[12:16].strip()
        element = line[76:78].strip() or name[0]
        xyz = tuple(float(value) for value in (line[30:38], line[38:46], line[46:54]))
        if name in atoms:
            duplicate_names.append(name)
        atoms[name] = (element, xyz)

    parameter_atoms: list[str] = []
    bonds: list[tuple[str, str]] = []
    for line in params.read_text(encoding="utf-8").splitlines():
        fields = line.split()
        if fields and fields[0] == "ATOM":
            parameter_atoms.append(fields[1])
        elif fields and fields[0] in {"BOND_TYPE", "BOND"}:
            bonds.append((fields[1], fields[2]))

    parameter_set, observed_set = set(parameter_atoms), set(atoms)
    missing = sorted(parameter_set - observed_set)
    missing_hydrogens = [name for name in missing if name.startswith("H")]
    missing_heavy = [name for name in missing if not name.startswith("H")]
    extra = sorted(observed_set - parameter_set)
    measured_bonds: list[dict[str, object]] = []
    outliers: list[dict[str, object]] = []
    for atom1, atom2 in bonds:
        if atom1 not in atoms or atom2 not in atoms:
            continue
        distance = math.dist(atoms[atom1][1], atoms[atom2][1])
        record = {"atom1": atom1, "atom2": atom2, "distance_A": distance}
        measured_bonds.append(record)
        involves_hydrogen = atom1.startswith("H") or atom2.startswith("H")
        lower, upper = (0.65, 1.25) if involves_hydrogen else (0.9, 1.9)
        if not lower <= distance <= upper:
            outliers.append(record)

    explicit_hydrogen_count = sum(element == "H" for element, _ in atoms.values())
    expected_hydrogen_count = sum(name.startswith("H") for name in parameter_atoms)
    valid = bool(parameter_atoms) and bool(atoms) and not missing and not extra and not duplicate_names and not outliers
    chemistry = {}
    if valid:
        from .ligand_chemistry import inspect_params_chemistry
        chemistry = inspect_params_chemistry(params, template=template)
        valid = valid and chemistry['valid']
        if valid:
            from .atom_roles import save_atom_roles
            save_atom_roles(params, template, chemistry)
    status = (
        "complete_explicit_hydrogens"
        if valid and not missing_hydrogens else
        "missing_hydrogens_rosetta_rebuild_required"
        if not missing_heavy and not extra and not duplicate_names and not outliers and missing_hydrogens else "invalid"
    )
    return {
        "pdb": str(pdb), "params": str(params), "chain": chain, "status": status,
        "valid": valid, "observed_atom_count": len(atoms),
        "chemical_topology": chemistry,
        "expected_atom_count": len(parameter_atoms),
        "explicit_hydrogen_count": explicit_hydrogen_count,
        "expected_hydrogen_count": expected_hydrogen_count,
        "missing_hydrogens": missing_hydrogens, "missing_heavy_atoms": missing_heavy,
        "extra_atoms": extra, "duplicate_atom_names": duplicate_names,
        "expected_bond_count": len(bonds), "measured_bond_count": len(measured_bonds),
        "bond_length_min_A": min((item["distance_A"] for item in measured_bonds), default=None),
        "bond_length_max_A": max((item["distance_A"] for item in measured_bonds), default=None),
        "bond_length_outliers": outliers,
    }
