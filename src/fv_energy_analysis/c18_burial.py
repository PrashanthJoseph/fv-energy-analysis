"""Protein-induced SASA loss of a chemically mapped estradiol C18 carbon.

The complete ligand is present in both SASA states. Protein atoms farther than
the sum of their and C18's probe-expanded radii cannot shield a C18 surface
point, so discarding only those atoms is exact for this atom's SASA. This
diagnostic has no packing acceptance threshold and does not change coordinates.
"""
from __future__ import annotations

from numbers import Integral
from pathlib import Path
from typing import Sequence
import math

import numpy as np

from .steric_sandwich import BONDI_RADII_A, StericAtom, read_steric_atoms


PROBE_RADIUS_A = 1.4
METHOD = "biopython_shrake_rupley_c18_exact_local_occluders_v1"


def _bio_structure(atoms: Sequence[StericAtom]):
    """Copy selected atoms into a public Biopython SASA input structure."""
    from Bio.PDB.Atom import Atom
    from Bio.PDB.Chain import Chain
    from Bio.PDB.Model import Model
    from Bio.PDB.Residue import Residue
    from Bio.PDB.Structure import Structure

    structure = Structure("c18_sasa")
    model = Model(0)
    structure.add(model)
    chains, residues, converted = {}, {}, []
    for index, original in enumerate(atoms):
        if original.chain not in chains:
            chains[original.chain] = Chain(original.chain)
            model.add(chains[original.chain])
        residue_id = (
            " " if original.record == "ATOM" else "H_" + original.residue,
            original.number,
            original.insertion or " ",
        )
        key = original.chain, residue_id
        if key not in residues:
            residues[key] = Residue(residue_id, original.residue, " ")
            chains[original.chain].add(residues[key])
        atom = Atom(
            original.name, np.array(original.xyz, dtype=float, copy=True),
            0.0, original.occupancy, " ", original.name, index + 1,
            element=original.element,
        )
        residues[key].add(atom)
        converted.append(atom)
    return structure, converted


def _validate_atoms(atoms: Sequence[StericAtom]) -> None:
    for atom in atoms:
        if atom.element not in BONDI_RADII_A:
            raise ValueError(f"Unsupported heavy-atom SASA radius: {atom.label} ({atom.element})")
        if np.shape(atom.xyz) != (3,) or not np.isfinite(atom.xyz).all():
            raise ValueError(f"Nonfinite or invalid atom coordinates: {atom.label}")


def _burial_from_atoms(
    ligand: Sequence[StericAtom],
    protein: Sequence[StericAtom],
    *,
    methyl_atom: str = "C10",
    n_points: int = 1920,
) -> dict:
    """Calculate C18 SASA; also supports small synthetic validation systems."""
    import Bio
    from Bio.PDB.SASA import ShrakeRupley

    if isinstance(n_points, bool) or not isinstance(n_points, Integral) or n_points < 1:
        raise ValueError("n_points must be a positive integer")
    ligand, protein = list(ligand), list(protein)
    _validate_atoms(ligand + protein)
    if len({atom.name for atom in ligand}) != len(ligand):
        raise ValueError("Selected ligand atom names must be unique")
    targets = [(index, atom) for index, atom in enumerate(ligand) if atom.name == methyl_atom]
    if len(targets) != 1 or targets[0][1].element != "C":
        raise ValueError(f"Expected one carbon named {methyl_atom} in the selected ligand")
    target_index, target = targets[0]

    target_radius = BONDI_RADII_A[target.element] + PROBE_RADIUS_A
    # The same expanded radii are used by ShrakeRupley. The small tolerance only
    # retains additional boundary atoms; it cannot remove a possible occluder.
    occluders = [
        atom for atom in protein
        if np.linalg.norm(atom.xyz - target.xyz)
        <= target_radius + BONDI_RADII_A[atom.element] + PROBE_RADIUS_A + 1e-12
    ]
    free_structure, free_atoms = _bio_structure(ligand)
    bound_structure, bound_atoms = _bio_structure(ligand + occluders)
    algorithm = ShrakeRupley(
        probe_radius=PROBE_RADIUS_A, n_points=int(n_points),
        radii_dict=dict(BONDI_RADII_A),
    )
    algorithm.compute(free_structure, level="A")
    algorithm.compute(bound_structure, level="A")
    isolated = float(free_atoms[target_index].sasa)
    bound = float(bound_atoms[target_index].sasa)
    if not math.isfinite(isolated) or isolated <= 0:
        raise ValueError("C18 has zero or invalid isolated-full-ligand SASA; burial fraction is undefined")
    if not math.isfinite(bound) or bound < 0 or bound > isolated + 1e-8:
        raise ValueError("Invalid C18 complex SASA subtraction")
    buried = isolated - bound
    fraction = buried / isolated
    return {
        "ligand_atom": methyl_atom,
        "ligand_residue_id": target.residue_id,
        "isolated_full_ligand_SASA_A2": isolated,
        "complex_SASA_A2": bound,
        "buried_SASA_A2": buried,
        "fraction_buried": fraction,
        "percent_buried": 100.0 * fraction,
        "ligand_heavy_atom_count": len(ligand),
        "protein_heavy_atom_count": len(protein),
        "protein_occluder_atom_count": len(occluders),
        "full_complex_heavy_atom_count": len(ligand) + len(protein),
        "measured_complex_heavy_atom_count": len(ligand) + len(occluders),
        "n_points_per_atom": int(n_points),
        "probe_radius_A": PROBE_RADIUS_A,
        "radii_A": dict(BONDI_RADII_A),
        "method": METHOD,
        "biopython_version": Bio.__version__,
        "scope": "C18 carbon only; ligand-side protein-induced SASA loss with full ligand self-occlusion",
        "occluder_selection": "distance <= sum of C18 and protein-atom probe-expanded Bondi radii",
        "full_ligand_present_in_both_states": True,
        "protein_backbone_included": True,
        "waters_and_hydrogens_included": False,
        "coordinates_changed": False,
        "acceptance_threshold": None,
    }


def analyze_c18_burial(
    pdb_path: str | Path,
    methyl_atom: str = "C10",
    ligand_chain: str = "X",
    ligand_resname: str = "E2L",
    protein_chains: Sequence[str] = ("H", "L"),
    n_points: int = 1920,
) -> dict:
    """Measure the mapped C18 carbon in a prepared twenty-heavy-atom E2 ligand.

    C10 is C18 only in this project's validated E2L atom map; callers using
    other atom naming must provide their chemically established methyl name.
    First-model/highest-occupancy-altloc selection follows the shared PDB reader.
    Protein ATOM records in the requested chains include backbone and sidechain.
    HETATM waters, other ligands, hydrogens and deuterium do not contribute.
    """
    atoms = read_steric_atoms(pdb_path)
    ligand = [
        atom for atom in atoms
        if atom.record == "HETATM" and atom.chain == ligand_chain
        and atom.residue == ligand_resname
    ]
    if len({atom.residue_id for atom in ligand}) != 1:
        raise ValueError("Expected exactly one selected E2 ligand residue")
    if len(ligand) != 20:
        raise ValueError(f"Expected exactly twenty E2 ligand heavy atoms; found {len(ligand)}")
    chains = set(protein_chains)
    protein = [atom for atom in atoms if atom.record == "ATOM" and atom.chain in chains]
    return _burial_from_atoms(ligand, protein, methyl_atom=methyl_atom, n_points=n_points)
