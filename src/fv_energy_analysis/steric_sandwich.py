"""Coordinate-only two-sided enclosure of a chemically mapped A-ring.

This project-defined screen was calibrated against the opposing walls described
by Nahoum et al. (2003), not supplied as a numerical definition by that paper.
It measures opposing side-chain contacts, not SASA, affinity, or specificity.
It does not call Discovery Studio or any interaction-fingerprint tool.
"""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Sequence
import math

import numpy as np


BONDI_RADII_A = {"C": 1.70, "N": 1.55, "O": 1.52, "S": 1.80}
BACKBONE_ATOMS = frozenset({"N", "CA", "C", "O", "OXT"})
METHOD_VERSION = "a_ring_two_sided_enclosure_v2"


@dataclass(frozen=True)
class StericAtom:
    name: str
    residue: str
    chain: str
    number: int
    insertion: str
    element: str
    xyz: np.ndarray
    occupancy: float = 1.0
    record: str = "ATOM"

    @property
    def residue_id(self) -> str:
        return f"{self.chain}:{self.number}{self.insertion}:{self.residue}"

    @property
    def label(self) -> str:
        return f"{self.residue_id}:{self.name}"


def read_steric_atoms(pdb_path: str | Path) -> list[StericAtom]:
    """Read the first PDB model; highest occupancy altloc, first wins ties.

    Hydrogens, deuterium and zero-occupancy records are omitted. Blank element
    fields are inferred from the first alphabetic atom-name character for the
    protein/steroid atoms supported here (not a general metal-element parser).
    """
    atoms: dict[tuple, StericAtom] = {}
    for line in Path(pdb_path).read_text().splitlines():
        if line.startswith("ENDMDL"):
            break
        if not line.startswith(("ATOM  ", "HETATM")):
            continue
        occupancy = float(line[54:60].strip() or "1")
        if occupancy <= 0:
            continue
        name = line[12:16].strip()
        element = line[76:78].strip().upper()
        if not element:
            element = next((c.upper() for c in name if c.isalpha()), "")
        if element in {"H", "D"}:
            continue
        atom = StericAtom(
            name, line[17:20].strip(), line[21:22].strip(),
            int(line[22:26]), line[26:27].strip(), element,
            np.array([float(line[i:i + 8]) for i in (30, 38, 46)]),
            occupancy, line[:6].strip(),
        )
        if not math.isfinite(occupancy) or not np.isfinite(atom.xyz).all():
            raise ValueError(f"Nonfinite PDB atom: {atom.label}")
        key = atom.chain, atom.number, atom.insertion, atom.residue, atom.name
        if key not in atoms or occupancy > atoms[key].occupancy:
            atoms[key] = atom
    return list(atoms.values())


def legacy_steric_wall_rays(
    ligand_atoms: Iterable[StericAtom],
    protein_atoms: Iterable[StericAtom],
    *,
    ring_atom_names: Sequence[str],
    methyl_atom: str | None = None,
    methyl_parent_atom: str | None = None,
    maximum_surface_gap_A: float = 1.5,
) -> dict:
    """Archived twelve-ray rule, retained only for historical calibration.

    Production uses :func:`steric_wall_coverage` below.

    Supply the six *chemical* A-ring carbon names explicitly. Standard steroid
    names are C1,C2,C3,C4,C5,C10; the validated project E2L map uses C1-C6.
    Optional C13-to-C18 methyl direction gives reproducible alpha/beta labels;
    without it the labels have arbitrary sign but the binary result is invariant.
    All intersected atoms and nearest wall identities are retained for audit.
    """
    names = tuple(ring_atom_names)
    if len(names) != 6 or len(set(names)) != 6:
        raise ValueError("Exactly six distinct chemically mapped A-ring atoms are required")
    if not math.isfinite(maximum_surface_gap_A) or maximum_surface_gap_A < 0:
        raise ValueError("maximum_surface_gap_A must be finite and nonnegative")
    if (methyl_atom is None) != (methyl_parent_atom is None):
        raise ValueError("Supply both methyl_atom and methyl_parent_atom, or neither")
    ligand = list(ligand_atoms)
    lookup = {a.name: a for a in ligand}
    if len(lookup) != len(ligand):
        raise ValueError("Ligand atom names must be unique after alternate-location selection")
    missing = set(names) - lookup.keys()
    if missing:
        raise ValueError(f"Missing mapped A-ring atoms: {sorted(missing)}")
    if any(lookup[name].element != "C" for name in names):
        raise ValueError("All six mapped A-ring atoms must be carbon")
    ring = np.array([lookup[name].xyz for name in names], dtype=float)
    if ring.shape != (6, 3) or not np.isfinite(ring).all():
        raise ValueError("Mapped A-ring coordinates must be finite 3-vectors")
    centroid = ring.mean(axis=0)
    _, singular, vh = np.linalg.svd(ring - centroid, full_matrices=False)
    if singular[1] < 1e-6:
        raise ValueError("Degenerate A-ring plane")
    normal = vh[-1]
    orientation = "arbitrary normal sign; alpha/beta labels are not stereochemical assignments"
    if methyl_atom is not None:
        if methyl_atom not in lookup or methyl_parent_atom not in lookup:
            raise ValueError("Mapped methyl orientation atoms are missing from ligand")
        projection = float(np.dot(lookup[methyl_atom].xyz - lookup[methyl_parent_atom].xyz, normal))
        if abs(projection) < 1e-8:
            raise ValueError("Mapped methyl direction cannot orient the A-ring plane normal")
        if projection < 0:
            normal = -normal
        orientation = f"beta normal follows {methyl_parent_atom}-to-{methyl_atom} methyl direction"
    walls = [a for a in protein_atoms
             if a.record == "ATOM" and a.name not in BACKBONE_ATOMS
             and a.element not in {"H", "D"} and a.occupancy > 0]
    unknown = sorted({a.element for a in walls} - BONDI_RADII_A.keys())
    if unknown:
        raise ValueError(f"No Bondi radius specified for protein side-chain elements: {unknown}")
    coordinates = np.array([a.xyz for a in walls], dtype=float).reshape((-1, 3))
    if not np.isfinite(coordinates).all():
        raise ValueError("Protein wall coordinates must be finite")
    radii = np.array([BONDI_RADII_A[a.element] for a in walls])
    rays = []
    minimum_gap = None
    for name, point in zip(names, ring):
        delta = coordinates - point
        distance2 = np.einsum("ij,ij->i", delta, delta)
        if walls:
            pair_gap = np.sqrt(distance2) - radii - BONDI_RADII_A["C"]
            this_minimum = float(pair_gap.min())
            minimum_gap = this_minimum if minimum_gap is None else min(minimum_gap, this_minimum)
        for face, direction in (("beta", normal), ("alpha", -normal)):
            height = delta @ direction
            lateral2 = np.maximum(0.0, distance2 - height ** 2)
            intersect = np.flatnonzero((height > 0) & (lateral2 <= radii ** 2))
            hits = []
            for index in intersect:
                entry = float(height[index] - math.sqrt(max(0.0, radii[index] ** 2 - lateral2[index])))
                atom = walls[index]
                hits.append({"wall_atom": atom.label, "wall_residue": atom.residue_id,
                             "wall_element": atom.element, "entry_A": entry,
                             "surface_gap_A": entry - BONDI_RADII_A["C"]})
            hits.sort(key=lambda item: (item["entry_A"], item["wall_atom"]))
            nearest = hits[0] if hits else {"wall_atom": None, "wall_residue": None,
                                           "wall_element": None, "entry_A": None, "surface_gap_A": None}
            gap = nearest["surface_gap_A"]
            rays.append({"ring_atom": name, "face": face, **nearest,
                         "near_wall": gap is not None and 0 <= gap <= maximum_surface_gap_A,
                         "overlap": gap is not None and gap < 0,
                         "intersections": hits})
    by_ray = {(r["ring_atom"], r["face"]): r for r in rays}
    by_face = {face: [name for name in names if by_ray[name, face]["near_wall"]]
               for face in ("alpha", "beta")}
    paired = [name for name in names if all(by_ray[name, face]["near_wall"] for face in by_face)
              and by_ray[name, "alpha"]["wall_residue"] != by_ray[name, "beta"]["wall_residue"]]
    overlaps = sum(ray["overlap"] for ray in rays)
    passed = min(map(len, by_face.values())) >= 2 and len(paired) >= 1 and overlaps == 0
    return {
        "method": "a_ring_sidechain_wall_rays_v1",
        "steric_sandwich": passed,
        "ring_atom_names": list(names),
        "ring_centroid_A": centroid.tolist(),
        "beta_normal": normal.tolist(),
        "face_orientation": orientation,
        "ring_plane_rms_deviation_A": float(np.sqrt(np.mean(((ring - centroid) @ normal) ** 2))),
        "coverage": {"alpha_carbon_rays": len(by_face["alpha"]),
                     "beta_carbon_rays": len(by_face["beta"]), "paired_carbon_rays": len(paired),
                     "overlapping_carbon_rays": overlaps,
                     "alpha_ring_atoms": by_face["alpha"], "beta_ring_atoms": by_face["beta"],
                     "paired_ring_atoms": paired},
        "thresholds": {"maximum_surface_gap_A": maximum_surface_gap_A,
                       "minimum_rays_each_face": 2, "minimum_paired_rays_distinct_residues": 1,
                       "maximum_overlapping_rays": 0},
        "rays": rays,
        "minimum_ring_sidechain_vdw_gap_A": minimum_gap,
        "largest_ring_sidechain_vdw_overlap_A": max(0.0, -minimum_gap) if minimum_gap is not None else 0.0,
        "protein_sidechain_heavy_atom_count": len(walls),
        "radii_A": dict(BONDI_RADII_A),
        "provenance": {
            "definition": "Project-defined exploratory paper-reference-calibrated geometric proxy, not a published numerical rule",
            "source_concept_doi": "10.1096/fj.02-0524fje",
            "radius_source_doi": "10.1021/j100785a001",
            "calibration": "Seven systems / sixteen crystallographic bound copies; largest required gap 1.335 A, rounded to 1.5 A",
            "limitations": "Discrete projected side-chain wall coverage, not SASA, percent burial, favorable energy, or validated affinity/specificity prediction; excludes backbone and waters",
            "implementation": "Coordinate geometry, not a native Discovery Studio interaction classification",
        },
    }


def steric_wall_coverage(
    ligand_atoms: Iterable[StericAtom],
    protein_atoms: Iterable[StericAtom],
    *,
    ring_atom_names: Sequence[str],
    methyl_atom: str | None = None,
    methyl_parent_atom: str | None = None,
    maximum_surface_gap_A: float = 1.0,
) -> dict:
    """Find distinct side chains enclosing opposite A-ring faces.

    A contact has a nearest ring-carbon surface gap in [0, maximum], and
    its atom sphere lies wholly on one side of the least-squares ring plane.
    An opposing pair qualifies when its connecting segment crosses the union
    of the six carbon van der Waals disks projected onto that plane. Every
    qualifying atom is searched, rather than preselecting one atom per residue.
    Overlapping atom pairs are reported independently and receive no credit.
    The optional methyl direction only orients alpha/beta labels.
    """
    names = tuple(ring_atom_names)
    if len(names) != 6 or len(set(names)) != 6:
        raise ValueError("Exactly six distinct chemically mapped A-ring atoms are required")
    if not math.isfinite(maximum_surface_gap_A) or maximum_surface_gap_A < 0:
        raise ValueError("maximum_surface_gap_A must be finite and nonnegative")
    if (methyl_atom is None) != (methyl_parent_atom is None):
        raise ValueError("Supply both methyl_atom and methyl_parent_atom, or neither")
    ligand = list(ligand_atoms)
    lookup = {a.name: a for a in ligand}
    if len(lookup) != len(ligand):
        raise ValueError("Ligand atom names must be unique after alternate-location selection")
    missing = set(names) - lookup.keys()
    if missing:
        raise ValueError(f"Missing mapped A-ring atoms: {sorted(missing)}")
    if any(lookup[name].element != "C" for name in names):
        raise ValueError("All six mapped A-ring atoms must be carbon")
    ring = np.array([lookup[name].xyz for name in names], dtype=float)
    if ring.shape != (6, 3) or not np.isfinite(ring).all():
        raise ValueError("Mapped A-ring coordinates must be finite 3-vectors")
    centroid = ring.mean(axis=0)
    _, singular, vh = np.linalg.svd(ring - centroid, full_matrices=False)
    if singular[1] < 1e-6:
        raise ValueError("Degenerate A-ring plane")
    normal = vh[-1]
    orientation = "arbitrary normal sign; alpha/beta labels are not stereochemical assignments"
    if methyl_atom is not None:
        if methyl_atom not in lookup or methyl_parent_atom not in lookup:
            raise ValueError("Mapped methyl orientation atoms are missing from ligand")
        direction = np.asarray(lookup[methyl_atom].xyz) - np.asarray(lookup[methyl_parent_atom].xyz)
        projection = float(np.dot(direction, normal))
        if not math.isfinite(projection) or abs(projection) < 1e-8:
            raise ValueError("Mapped methyl direction cannot orient the A-ring plane normal")
        if projection < 0:
            normal = -normal
        orientation = f"beta normal follows {methyl_parent_atom}-to-{methyl_atom} methyl direction"
    walls = sorted((a for a in protein_atoms
                    if a.record == "ATOM" and a.name not in BACKBONE_ATOMS
                    and a.element not in {"H", "D"} and a.occupancy > 0),
                   key=lambda a: a.label)
    unknown = sorted({a.element for a in walls} - BONDI_RADII_A.keys())
    if unknown:
        raise ValueError(f"No Bondi radius specified for protein side-chain elements: {unknown}")
    coordinates = np.array([a.xyz for a in walls], dtype=float).reshape((-1, 3))
    if not np.isfinite(coordinates).all() or any(not math.isfinite(a.occupancy) for a in walls):
        raise ValueError("Protein wall coordinates and occupancies must be finite")
    radii = np.array([BONDI_RADII_A[a.element] for a in walls])
    carbon_radius = BONDI_RADII_A["C"]
    tolerance = 1e-9  # Roundoff only, not an additional physical allowance.
    distances = np.linalg.norm(coordinates[:, None, :] - ring[None, :, :], axis=2)
    gaps = distances - radii[:, None] - carbon_radius
    heights = (coordinates - centroid) @ normal
    contacts, overlaps = [], []
    for index, atom in enumerate(walls):
        for ring_index in np.flatnonzero(gaps[index] < -tolerance):
            overlaps.append({"wall_atom": atom.label, "wall_residue": atom.residue_id,
                             "wall_element": atom.element, "ring_atom": names[ring_index],
                             "distance_A": float(distances[index, ring_index]),
                             "surface_gap_A": float(gaps[index, ring_index])})
        nearest = int(np.argmin(distances[index]))
        gap = float(gaps[index, nearest])
        if not (-tolerance <= gap <= maximum_surface_gap_A + tolerance
                and abs(heights[index]) + tolerance >= radii[index]):
            continue
        relative = coordinates[index] - ring[nearest]
        axial = float(np.dot(relative, normal))
        lateral = float(np.linalg.norm(relative - axial * normal))
        contacts.append({"wall_atom": atom.label, "wall_residue": atom.residue_id,
                         "wall_element": atom.element, "xyz_A": atom.xyz.tolist(),
                         "nearest_ring_atom": names[nearest],
                         "distance_A": float(distances[index, nearest]), "surface_gap_A": gap,
                         "signed_height_A": float(heights[index]),
                         "face": "alpha" if heights[index] < 0 else "beta",
                         "lateral_offset_A": lateral,
                         "angle_to_nearest_carbon_normal_degrees": math.degrees(math.atan2(lateral, abs(axial)))})
    by_face = {face: [c for c in contacts if c["face"] == face] for face in ("alpha", "beta")}
    opposing_pairs = []
    residue_pairs = set()
    for alpha in by_face["alpha"]:
        for beta in by_face["beta"]:
            if alpha["wall_residue"] == beta["wall_residue"]:
                continue
            a, b = np.array(alpha["xyz_A"]), np.array(beta["xyz_A"])
            fraction = -alpha["signed_height_A"] / (beta["signed_height_A"] - alpha["signed_height_A"])
            crossing = a + fraction * (b - a)
            delta = ring - crossing
            projected_delta = delta - (delta @ normal)[:, None] * normal
            nearest = float(np.linalg.norm(projected_delta, axis=1).min())
            if nearest > carbon_radius + tolerance:
                continue
            va, vb = a - centroid, b - centroid
            cosine = float(np.dot(va, vb) / (np.linalg.norm(va) * np.linalg.norm(vb)))
            opposing_pairs.append({"alpha": alpha, "beta": beta,
                                   "plane_crossing_xyz_A": crossing.tolist(),
                                   "plane_crossing_nearest_projected_ring_carbon_A": nearest,
                                   "plane_crossing_inside_projected_carbon_vdw_footprint": True,
                                   "centroid_subtended_angle_degrees": math.degrees(math.acos(float(np.clip(cosine, -1, 1)))),
                                   "protein_atom_separation_A": float(np.linalg.norm(b - a))})
            residue_pairs.add(tuple(sorted((alpha["wall_residue"], beta["wall_residue"]))))
    residues = {face: sorted({c["wall_residue"] for c in by_face[face]}) for face in by_face}
    minimum_gap = float(gaps.min()) if walls else None
    return {
        "method": METHOD_VERSION, "steric_sandwich": bool(opposing_pairs),
        "ring_atom_names": list(names), "ring_centroid_A": centroid.tolist(),
        "beta_normal": normal.tolist(), "face_orientation": orientation,
        "ring_plane_rms_deviation_A": float(np.sqrt(np.mean(((ring - centroid) @ normal) ** 2))),
        "coverage": {"alpha_residue_count": len(residues["alpha"]),
                     "beta_residue_count": len(residues["beta"]),
                     "opposing_residue_pair_count": len(residue_pairs),
                     "overlapping_atom_pair_count": len(overlaps),
                     "alpha_residues": residues["alpha"], "beta_residues": residues["beta"],
                     "opposing_residue_pairs": [list(pair) for pair in sorted(residue_pairs)]},
        "thresholds": {"minimum_surface_gap_A": 0.0, "maximum_surface_gap_A": maximum_surface_gap_A,
                       "minimum_face_height": "absolute signed height >= protein atom Bondi radius",
                       "minimum_opposing_distinct_residue_pairs": 1,
                       "plane_crossing": "inside union of six projected carbon van der Waals disks",
                       "numerical_tolerance_A": tolerance},
        "contacts": contacts, "opposing_pairs": opposing_pairs, "overlaps": overlaps,
        "minimum_ring_sidechain_vdw_gap_A": minimum_gap,
        "largest_ring_sidechain_vdw_overlap_A": max(0.0, -minimum_gap) if minimum_gap is not None else 0.0,
        "protein_sidechain_heavy_atom_count": len(walls), "radii_A": dict(BONDI_RADII_A),
        "provenance": {
            "definition": "Project-defined paper-reference-compatible two-sided enclosure, not a published numerical rule",
            "source_concept_doi": "10.1096/fj.02-0524fje", "radius_source_doi": "10.1021/j100785a001",
            "calibration": "Seven systems / sixteen bound reference copies; largest named-wall gap 0.844 A, rounded to 1.0 A",
            "limitations": "Two-sided enclosure, not SASA, favorable energy, C19 obstruction, or validated E2 specificity; excludes backbone and waters",
            "implementation": "Coordinate geometry, not a native Discovery Studio interaction classification",
        },
    }


def analyze_steric_sandwich(
    pdb_path: str | Path,
    *,
    ligand_chain: str,
    ligand_resname: str,
    ring_atom_names: Sequence[str],
    protein_chains: Sequence[str],
    ligand_resseq: int | None = None,
    ligand_icode: str | None = None,
    methyl_atom: str | None = None,
    methyl_parent_atom: str | None = None,
    maximum_surface_gap_A: float = 1.0,
) -> dict:
    """PDB entry point, with explicit ligand identity and protein-chain scope."""
    chains = set(protein_chains)
    if not chains:
        raise ValueError("At least one explicit protein chain is required")
    atoms = read_steric_atoms(pdb_path)
    ligand = [a for a in atoms if a.chain == ligand_chain and a.residue == ligand_resname
              and (ligand_resseq is None or a.number == ligand_resseq)
              and (ligand_icode is None or a.insertion == ligand_icode)]
    identities = {(a.chain, a.number, a.insertion, a.residue) for a in ligand}
    if len(identities) != 1:
        raise ValueError(f"Expected one ligand residue for {ligand_chain}:{ligand_resname}; found {len(identities)}")
    protein = [a for a in atoms if a.chain in chains and a.record == "ATOM"
               and (a.chain, a.number, a.insertion, a.residue) not in identities]
    result = steric_wall_coverage(ligand, protein, ring_atom_names=ring_atom_names,
                                 methyl_atom=methyl_atom, methyl_parent_atom=methyl_parent_atom,
                                 maximum_surface_gap_A=maximum_surface_gap_A)
    result["ligand_residue"] = ligand[0].residue_id
    result["protein_chains"] = sorted(chains)
    result["input_pdb"] = str(Path(pdb_path))
    return result
