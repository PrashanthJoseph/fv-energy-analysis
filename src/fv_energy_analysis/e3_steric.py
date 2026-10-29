"""Fixed-pocket heavy-atom screen for a correctly placed estriol 16alpha-OH.

The caller supplies an E3 probe from e3_chemistry, not a literal C18 hydroxyl.
This module neither changes coordinates nor calculates energies. It reports
continuous Bondi gaps and a conservative overlap warning, not E3 specificity.
"""
from __future__ import annotations

from hashlib import sha256
import math
from pathlib import Path
from typing import Sequence

import numpy as np

from .steric_sandwich import BACKBONE_ATOMS, BONDI_RADII_A, StericAtom, read_steric_atoms


METHOD_VERSION = "e3_16alpha_added_oxygen_fixed_pocket_hbond_geometry_v2"
SCOPE = "fixed_pocket_steric_only"
MOLPROBITY_REFERENCE = "https://pmc.ncbi.nlm.nih.gov/articles/PMC2803126/"
SHBG_REFERENCE = "https://pmc.ncbi.nlm.nih.gov/articles/PMC305588/"


def _may_form_hbond_with_hydroxyl(atom: StericAtom) -> bool:
    """Route N/O/S to explicit chemistry and geometry, without exempting them."""
    return atom.element in {"N", "O", "S"}


def _saved_hydrogen_orientations(hydrogen_positions, current_hydrogen):
    """Normalize saved positions without generating or moving any coordinates."""
    if isinstance(hydrogen_positions, dict):
        hydrogen_positions = hydrogen_positions.get("orientations")
    if hydrogen_positions is None:
        return [{"orientation_index": 0,
                 "hydrogen_xyz_A": None if current_hydrogen is None else current_hydrogen.xyz.tolist(),
                 "steric_eligibility": None}], 0 if current_hydrogen is not None else None
    rows = []
    for index, original in enumerate(hydrogen_positions):
        row = dict(original) if isinstance(original, dict) else {"hydrogen_xyz_A": original}
        xyz = np.asarray(row.get("hydrogen_xyz_A"), dtype=float)
        if xyz.shape != (3,) or not np.isfinite(xyz).all():
            raise ValueError("Saved hydroxyl-H positions must contain finite XYZ triplets")
        # The corrected exported probe uses standard PDB 0.001 A coordinates;
        # evaluate precisely that saved position rather than a different
        # full-precision point that could flip an angular cutoff on export.
        row["hydrogen_xyz_A"] = np.round(xyz, 3).tolist()
        row.setdefault("orientation_index", index)
        has_eligibility = "nonpolar_warning_atom_count" in row and "short_polar_contact_count" in row
        row["steric_eligibility"] = (row["nonpolar_warning_atom_count"] == 0
                                      and row["short_polar_contact_count"] == 0) if has_eligibility else None
        rows.append(row)
    if not rows or len({row["orientation_index"] for row in rows}) != len(rows):
        raise ValueError("Saved hydroxyl-H orientations must be nonempty and have unique indices")
    persisted = None
    if current_hydrogen is not None:
        closest = min(rows, key=lambda r: float(np.linalg.norm(np.asarray(r["hydrogen_xyz_A"])-current_hydrogen.xyz)))
        if float(np.linalg.norm(np.asarray(closest["hydrogen_xyz_A"])-current_hydrogen.xyz)) <= .002:
            persisted = closest["orientation_index"]
    return rows, persisted


def analyze_added_oxygen_contacts(
    protein: Sequence[StericAtom],
    added_o_xyz_A: Sequence[float],
    *,
    overlap_warning_A: float = 0.4,
    contact_cutoff_A: float = 6.0,
    polar_evaluations: dict[str, dict] | None = None,
) -> dict:
    """Measure nonbonded contacts of one added oxygen against protein atoms.

    Radius-overlap >= 0.4 A is inspired by the published non-donor/acceptor
    MolProbity convention. Bondi heavy-atom spheres here are not the PROBE
    algorithm, radii or all-atom clashscore. Polar overlaps are exempt only when
    the supplied chemistry/directional-H geometry supports an H-bond. Missing
    geometry is unassessed, rather than a negative warning or an exemption.
    """
    xyz = np.asarray(added_o_xyz_A, dtype=float)
    if xyz.shape != (3,) or not np.isfinite(xyz).all():
        raise ValueError("Added oxygen must have three finite coordinates")
    if not math.isfinite(overlap_warning_A) or overlap_warning_A <= 0:
        raise ValueError("overlap_warning_A must be positive and finite")
    if not math.isfinite(contact_cutoff_A) or contact_cutoff_A <= 0:
        raise ValueError("contact_cutoff_A must be positive and finite")
    protein = list(protein)
    if not protein:
        raise ValueError("No protein heavy atoms selected")
    unsupported = sorted({a.element for a in protein} - BONDI_RADII_A.keys())
    if unsupported:
        raise ValueError(f"Unsupported protein heavy-atom radii: {unsupported}")
    if any(np.shape(a.xyz) != (3,) or not np.isfinite(a.xyz).all() for a in protein):
        raise ValueError("Nonfinite or invalid protein coordinates")
    first_residue = {}
    for atom in protein:
        first_residue.setdefault(atom.chain, atom.residue_id)
    rows = []
    for atom in protein:
        distance = float(np.linalg.norm(atom.xyz - xyz))
        gap = distance - BONDI_RADII_A["O"] - BONDI_RADII_A[atom.element]
        overlap = max(0.0, -gap)
        potential_polar = _may_form_hbond_with_hydroxyl(atom)
        legacy_nonpolar = (not potential_polar or
                           (atom.residue == "PRO" and atom.name == "N"
                            and atom.residue_id != first_residue[atom.chain]))
        above = overlap >= overlap_warning_A - 1e-12
        geometry = (polar_evaluations or {}).get(atom.label) if above and potential_polar else None
        if geometry is not None and geometry.get("state") not in {"supported", "not_supported", "unassessed"}:
            raise ValueError(f"Invalid H-bond assessment for {atom.label}: {geometry.get('state')!r}")
        state = geometry["state"] if geometry is not None else "unassessed" if potential_polar else "not_supported"
        assessment = ("below_overlap_warning" if not above else
                      "hbond_supported_overlap" if state == "supported" else
                      "non_hbond_overlap" if state == "not_supported" else
                      "unassessed_polar_overlap")
        rows.append({
            "protein_residue": atom.residue_id, "protein_atom": atom.label,
            "protein_atom_name": atom.name, "protein_element": atom.element,
            "is_backbone": atom.name in BACKBONE_ATOMS,
            "distance_A": distance, "vdw_surface_gap_A": gap,
            "vdw_overlap_A": overlap, "overlap_at_or_above_warning": above,
            "potential_polar_contact": potential_polar,
            "legacy_nonpolar_diagnostic": legacy_nonpolar,
            "contact_assessment": assessment,
            "hbond_assessed": above and potential_polar and state != "unassessed",
            "hbond_geometry": geometry,
        })
    rows.sort(key=lambda r: (r["distance_A"], r["protein_atom"]))
    closest = rows[0]
    overlaps = [r for r in rows if r["overlap_at_or_above_warning"]]
    nonpolar = [r for r in overlaps if r["legacy_nonpolar_diagnostic"]]
    polar = [r for r in overlaps if not r["legacy_nonpolar_diagnostic"]]
    definite = [r for r in overlaps if r["contact_assessment"] == "non_hbond_overlap"]
    supported = [r for r in overlaps if r["contact_assessment"] == "hbond_supported_overlap"]
    unassessed = [r for r in overlaps if r["contact_assessment"] == "unassessed_polar_overlap"]
    residue_labels = lambda selected: sorted({r["protein_residue"] for r in selected})
    warning = True if definite else None if unassessed else False
    assessment = ("steric_overlap_detected" if definite else
                  "unassessed_polar_geometry" if unassessed else
                  "hbond_supported_overlaps_only" if supported else
                  "no_large_overlap_in_fixed_pose")
    contacts = [r for r in rows if r["distance_A"] <= contact_cutoff_A
                or r["overlap_at_or_above_warning"]]
    return {
        "schema_version": 2, "status": "complete",
        "summary": {
            "e3_steric_warning": warning,
            "e3_steric_assessment": assessment,
            "e3_added_o_closest_distance_A": closest["distance_A"],
            "e3_added_o_closest_residue": closest["protein_residue"],
            "e3_added_o_closest_atom": closest["protein_atom"],
            "e3_added_o_deepest_overlap_A": max(r["vdw_overlap_A"] for r in rows),
            "e3_added_o_overlap_residue_count": len(residue_labels(overlaps)),
            "e3_added_o_overlap_residues": ";".join(residue_labels(overlaps)) or "none",
            "e3_added_o_nonpolar_overlap_residue_count": len(residue_labels(nonpolar)),
            "e3_added_o_nonpolar_overlap_residues": ";".join(residue_labels(nonpolar)) or "none",
            "e3_added_o_polar_close_residue_count": len(residue_labels(polar)),
            "e3_added_o_polar_close_residues": ";".join(residue_labels(polar)) or "none",
            "e3_added_o_non_hbond_overlap_residue_count": len(residue_labels(definite)),
            "e3_added_o_non_hbond_overlap_residues": ";".join(residue_labels(definite)) or "none",
            "e3_added_o_hbond_supported_overlap_residue_count": len(residue_labels(supported)),
            "e3_added_o_hbond_supported_overlap_residues": ";".join(residue_labels(supported)) or "none",
            "e3_added_o_unassessed_overlap_residue_count": len(residue_labels(unassessed)),
            "e3_added_o_unassessed_overlap_residues": ";".join(residue_labels(unassessed)) or "none",
            "e3_analysis_scope": SCOPE,
        },
        "contacts": contacts,
        "closest_contact": closest,
        "protein_heavy_atom_count": len(protein),
        "added_o_xyz_A": xyz.tolist(),
        "protocol": {
            "method": METHOD_VERSION, "added_group": "estriol 16alpha-OH oxygen",
            "scope": SCOPE, "protein_backbone_included": True,
            "hydrogens_included": False, "waters_included": False,
            "coordinates_changed": False, "protein_repacking": False,
            "minimization": False, "energies_assessed": False,
            "hbond_geometry_assessed": polar_evaluations is not None,
            "hydrogens_used_for_hbond_geometry": polar_evaluations is not None,
            "radii_A": dict(BONDI_RADII_A),
            "gap_formula": "distance - radius_added_oxygen - radius_protein_atom",
            "warning_overlap_A": overlap_warning_A,
            "contact_output_cutoff_A": contact_cutoff_A,
            "warning_rule": "TRUE when any >= warning_overlap_A pair is not H-bond compatible; FALSE only when every large overlap is H-bond-supported or no large overlap exists; null when unresolved polar geometry remains with no definite clash",
            "polar_policy": "N/O/S receive no blanket exemption: residue/atom donor-acceptor chemistry and explicit donor-H directional geometry must support an H-bond; missing information remains unassessed",
            "legacy_nonpolar_diagnostic": "Historical C/internal-PRO-N count is retained for reference only; the current warning uses the chemistry-and-geometry-assessed non-Hbond overlap count",
            "geometry_combination_policy": "Per-pair geometric support across available hydroxyl-H positions does not prove all supported interactions coexist in one hydroxyl orientation",
            "cutoff_provenance": "0.4 A is a project heavy-atom warning inspired by MolProbity non-donor/acceptor overlap convention; not an E3-calibrated cutoff or a PROBE clashscore",
            "references": {"steric_overlap_precedent": MOLPROBITY_REFERENCE,
                           "estriol_16alpha_steric_precedent": SHBG_REFERENCE},
            "limitations": [
                "A fixed-pose warning does not establish E2/E3 binding specificity.",
                "A lack of warning does not establish E3 accommodation; energetic costs and alternate binding poses are unassessed.",
                "The added oxygen is tested at the chemically validated graft position; no conformational adaptation is sampled.",
            ],
        },
    }


def analyze_e3_steric(
    pdb_path: str | Path,
    added_o_atom: str = "O3",
    ligand_chain: str = "X",
    ligand_resname: str = "E3L",
    protein_chains: Sequence[str] = ("H", "L"),
    overlap_warning_A: float = 0.4,
    contact_cutoff_A: float = 6.0,
    added_o_xyz_A: Sequence[float] | None = None,
    hydrogen_positions=None,
    added_h_atom: str = "H25",
) -> dict:
    """Analyze an exported E3 probe, or coordinates supplied by its builder.

    Supplying coordinates permits measurement against the unchanged source E2
    complex. The caller remains responsible for verified steroid atom mapping
    and 16alpha stereochemistry; this function does not infer them from names.
    """
    from .e3_hbond import classify_polar_overlap

    path = Path(pdb_path)
    _, records = _read_all_atom_records(path.read_text(encoding="utf-8"))
    all_atoms = [a for a, _ in records]
    atoms = [a for a in all_atoms if a.element != "H"]
    ligand_all = [a for a in all_atoms if a.record == "HETATM"
                  and a.chain == ligand_chain and a.residue == ligand_resname]
    if added_o_xyz_A is None:
        ligand = [a for a in atoms if a.record == "HETATM"
                  and a.chain == ligand_chain and a.residue == ligand_resname]
        if len({a.residue_id for a in ligand}) != 1 or len(ligand) != 21:
            raise ValueError("Expected one twenty-one-heavy-atom E3 probe ligand")
        selected = [a for a in ligand if a.name == added_o_atom]
        if len(selected) != 1 or selected[0].element != "O":
            raise ValueError("Expected one mapped added hydroxyl oxygen")
        added_o_xyz_A = selected[0].xyz
    chains = set(protein_chains)
    protein = [a for a in atoms if a.record == "ATOM" and a.chain in chains]
    protein_all = [a for a in all_atoms if a.record == "ATOM" and a.chain in chains]
    current = [a for a in ligand_all if a.name == added_h_atom and a.element == "H"]
    if len(current) > 1:
        raise ValueError("Expected at most one mapped hydroxyl H")
    orientations, persisted_index = _saved_hydrogen_orientations(
        hydrogen_positions, current[0] if current else None)
    xyz = np.asarray(added_o_xyz_A, dtype=float)
    if xyz.shape != (3,) or not np.isfinite(xyz).all():
        raise ValueError("Added oxygen must have three finite coordinates")
    if not protein:
        raise ValueError("No protein heavy atoms selected")
    unknown = sorted({a.element for a in protein} - BONDI_RADII_A.keys())
    if unknown:
        raise ValueError(f"Unsupported protein heavy-atom radii: {unknown}")
    coordinates = np.array([a.xyz for a in protein])
    distances = np.linalg.norm(coordinates-xyz, axis=1)
    overlaps = (BONDI_RADII_A["O"] + np.array([BONDI_RADII_A[a.element] for a in protein])
                - distances)
    large = [a for a, overlap in zip(protein, overlaps) if overlap >= overlap_warning_A-1e-12]
    fixed_nonpolar = [a for a in large if not _may_form_hbond_with_hydroxyl(a)]
    large_polar = [a for a in large if _may_form_hbond_with_hydroxyl(a)]
    # Protein-donor -> added-O support is unchanged by new-OH torsion, but the
    # reverse direction must be evaluated for one shared H position at a time.
    evaluated = []
    for orientation in orientations:
        hxyz = orientation["hydrogen_xyz_A"]
        geometry = {}
        for atom in large_polar:
            value = dict(classify_polar_overlap(
                atom, protein_all, xyz.tolist(),
                added_h_positions_A=[] if hxyz is None else [hxyz]))
            state = value.get("state", value.get("status"))
            if state not in {"supported", "not_supported", "unassessed"}:
                raise ValueError(f"Invalid H-bond geometry state for {atom.label}: {state!r}")
            value["state"] = state
            geometry[atom.label] = value
        definite = fixed_nonpolar + [a for a in large_polar if geometry[a.label]["state"] == "not_supported"]
        unassessed = [a for a in large_polar if geometry[a.label]["state"] == "unassessed"]
        supported = [a for a in large_polar if geometry[a.label]["state"] == "supported"]
        evaluated.append({
            "orientation_index": orientation["orientation_index"],
            "hydrogen_xyz_A": hxyz,
            "steric_eligibility": orientation["steric_eligibility"],
            "certain_clash_atom_count": len(definite),
            "certain_clash_residue_count": len({a.residue_id for a in definite}),
            "unassessed_atom_count": len(unassessed),
            "unassessed_residue_count": len({a.residue_id for a in unassessed}),
            "hbond_supported_atom_count": len(supported),
            "hbond_supported_residue_count": len({a.residue_id for a in supported}),
            "polar_evaluations": geometry,
        })
    eligible = [row for row in evaluated if row["steric_eligibility"] is True]
    candidates = eligible or evaluated
    selected = min(candidates, key=lambda row: (
        row["certain_clash_atom_count"], row["unassessed_atom_count"],
        row["orientation_index"] != persisted_index, row["orientation_index"]))
    result = analyze_added_oxygen_contacts(
        protein, added_o_xyz_A, overlap_warning_A=overlap_warning_A,
        contact_cutoff_A=contact_cutoff_A, polar_evaluations=selected["polar_evaluations"],
    )
    result["selected_hydrogen_orientation_index"] = selected["orientation_index"]
    result["selected_hydrogen_xyz_A"] = selected["hydrogen_xyz_A"]
    result["hydrogen_orientation_selection"] = {
        "orientation_count": len(evaluated), "persisted_orientation_index": persisted_index,
        "eligible_orientation_count": len(eligible),
        "steric_eligibility_assessed": any(r["steric_eligibility"] is not None for r in evaluated),
        "eligibility_policy": "Use saved orientations with no nonpolar H warnings and no short-polar H contacts when any exist; otherwise retain and report the H-channel limitation separately",
        "selection": "One common H orientation minimizes certain O-overlap clashes, then unassessed O-overlap pairs; persisted orientation breaks ties",
        "per_orientation": [{k: v for k, v in row.items() if k != "polar_evaluations"} for row in evaluated],
    }
    result["protocol"]["geometry_combination_policy"] = (
        "All O-overlap exemptions use ONE common new-OH hydrogen orientation. "
        "H-bonds from fixed protein donors to added O do not depend on this torsion. "
        "No atom is moved by this analysis; the selected saved H position is returned separately.")
    result["protocol"]["hbond_geometry_scope"] = "Chemistry-aware, explicit donor-H geometry using the helper protocol; no DS, PLIP, Rosetta, minimization or repacking executed"
    result["protocol"]["saved_hydrogen_coordinate_precision_A"] = .001
    result["protocol"]["hbond_geometry_criteria"] = {
        "donor_acceptor_distance_A": "< 4.1",
        "angle_DHA_degrees": "> 100",
        "hydrogen_acceptor_distance_A": ">= 1.2",
        "criteria_provenance": "PLIP-style donor-acceptor distance and D-H-A angle screen plus project 1.2 A short-contact guard; geometry postprocessing, not execution of PLIP",
        "short_contact_guard_status": "Project safeguard, not an experimentally validated E3 discrimination cutoff",
    }
    result["added_o_atom"] = added_o_atom
    result["provenance"] = {"pdb": str(path.resolve()),
                            "pdb_sha256": sha256(path.read_bytes()).hexdigest(),
                            "protein_chains": list(protein_chains)}
    return result


def _read_all_atom_records(text: str):
    """First-model, highest-occupancy PDB atoms, retaining explicit H and lines."""
    selected = {}
    lines = text.splitlines(keepends=True)
    for line_index, line in enumerate(lines):
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
        if element == "D":
            element = "H"
        atom = StericAtom(
            name, line[17:20].strip(), line[21:22].strip(), int(line[22:26]),
            line[26:27].strip(), element,
            np.array([float(line[p:p+8]) for p in (30, 38, 46)]),
            occupancy=occupancy, record=line[:6].strip(),
        )
        key = atom.record, atom.chain, atom.number, atom.insertion, atom.residue, atom.name
        if key not in selected or selected[key][0].occupancy < occupancy:
            selected[key] = atom, line_index
    return lines, list(selected.values())


def _potential_acceptor_for_hydroxyl_h(atom: StericAtom, atoms: Sequence[StericAtom]) -> bool:
    """Conservative standard-residue acceptor possibility, not an H-bond call."""
    if atom.element in {"O", "S"}:
        return True
    if atom.element != "N":
        return False
    if atom.residue not in {"HIS", "HID", "HIE", "HIP"} or atom.name not in {"ND1", "NE2"}:
        # Protein backbone/amide, Lys, Arg and Trp nitrogens are not acceptors.
        return False
    return not any(a.element == "H" and a.residue_id == atom.residue_id
                   and float(np.linalg.norm(a.xyz - atom.xyz)) <= 1.25 for a in atoms)


def _hydroxyl_h_contact_score(
    position: np.ndarray,
    targets: Sequence[StericAtom],
    *,
    overlap_warning_A: float,
    short_polar_distance_A: float,
    target_arrays=None,
) -> dict:
    radii = {**BONDI_RADII_A, "H": 1.2}
    if target_arrays is None:
        target_arrays = (
            np.array([atom.xyz for atom in targets]),
            np.array([radii[atom.element] for atom in targets]),
            np.array([_potential_acceptor_for_hydroxyl_h(atom, targets) for atom in targets], dtype=bool),
            np.array([atom.record == "ATOM" for atom in targets], dtype=bool),
        )
    coordinates, target_radii, acceptors, is_protein = target_arrays
    distances = np.linalg.norm(coordinates-position, axis=1)
    overlaps = np.maximum(0.0, 1.0+target_radii-distances)
    nonpolar = ~acceptors
    protein_nonpolar = nonpolar & is_protein
    short = acceptors & (distances < short_polar_distance_A)
    warnings = nonpolar & (overlaps >= overlap_warning_A - 1e-12)
    contacts = []
    for index in np.flatnonzero((distances <= 3.5) | short):
        atom = targets[index]
        distance, overlap, acceptor = float(distances[index]), float(overlaps[index]), bool(acceptors[index])
        contacts.append({
            "atom": atom.label, "element": atom.element,
            "record": atom.record, "distance_A": distance,
            "vdw_overlap_A": overlap,
            "potential_hbond_acceptor": acceptor,
            "short_polar_contact": acceptor and distance < short_polar_distance_A,
            "nonpolar_overlap_warning": not acceptor and overlap >= overlap_warning_A - 1e-12,
        })
    return {
        "maximum_nonpolar_overlap_A": float(np.max(overlaps[nonpolar], initial=0.0)),
        "maximum_protein_nonpolar_overlap_A": float(np.max(overlaps[protein_nonpolar], initial=0.0)),
        "nonpolar_warning_atom_count": int(warnings.sum()),
        "short_polar_contact_count": int(short.sum()),
        "short_polar_deficit_sum_A": float(np.sum(short_polar_distance_A-distances[short])),
        "nonpolar_overlap_square_sum_A2": float(np.sum(overlaps[nonpolar]**2)),
        "contacts": sorted((r for r in contacts if r["distance_A"] <= 3.5
                             or r["short_polar_contact"]), key=lambda r: (r["distance_A"], r["atom"])),
    }


def orient_e3_hydroxyl(
    probe_pdb: str | Path,
    added_o_atom: str = "O3",
    added_h_atom: str = "H25",
    c16_atom: str = "C12",
    n_orientations: int = 12,
    *,
    ligand_chain: str = "X",
    ligand_resname: str = "E3L",
    protein_chains: Sequence[str] = ("H", "L"),
    output_pdb: str | Path | None = None,
    overlap_warning_A: float = 0.4,
    short_polar_distance_A: float = 1.2,
) -> dict:
    """Sample only the new hydroxyl H around its fixed C16–O bond.

    This is deterministic steric sampling, not minimization, repacking, energy
    optimization or a Boltzmann ensemble. Only the added H line is modified.
    Full neutral E3 must already contain 21 heavy atoms and 24 hydrogens.
    """
    if isinstance(n_orientations, bool) or not isinstance(n_orientations, int) or n_orientations < 1:
        raise ValueError("n_orientations must be a positive integer")
    if not math.isfinite(overlap_warning_A) or overlap_warning_A <= 0:
        raise ValueError("overlap_warning_A must be positive and finite")
    if not math.isfinite(short_polar_distance_A) or short_polar_distance_A <= 0:
        raise ValueError("short_polar_distance_A must be positive and finite")
    source = Path(probe_pdb)
    original = source.read_bytes()
    lines, records = _read_all_atom_records(original.decode("utf-8"))
    ligand_records = [(a, i) for a, i in records if a.record == "HETATM"
                      and a.chain == ligand_chain and a.residue == ligand_resname]
    ligand = [a for a, _ in ligand_records]
    if (len({a.residue_id for a in ligand}) != 1 or len(ligand) != 45
            or sum(a.element == "H" for a in ligand) != 24
            or sum(a.element != "H" for a in ligand) != 21):
        raise ValueError("Expected neutral full-H E3: one ligand with 21 heavy atoms and 24 hydrogens")
    def selected(name, element):
        result = [(a, i) for a, i in ligand_records if a.name == name]
        if len(result) != 1 or result[0][0].element != element:
            raise ValueError(f"Expected exactly one mapped {element} atom {name}")
        return result[0]
    oxygen, _ = selected(added_o_atom, "O")
    hydrogen, hydrogen_line = selected(added_h_atom, "H")
    carbon, _ = selected(c16_atom, "C")
    axis = carbon.xyz - oxygen.xyz
    length = float(np.linalg.norm(axis))
    if not np.isfinite(axis).all() or length < 1e-8:
        raise ValueError("Invalid C16–O hydroxyl axis")
    axis /= length
    # Anchor the torsion lattice to the builder's original O–H projection,
    # rather than global XYZ axes, so rigid transforms preserve the sampling.
    projection = hydrogen.xyz - oxygen.xyz
    projection -= float(np.dot(projection, axis)) * axis
    norm = float(np.linalg.norm(projection))
    if not np.isfinite(projection).all() or norm < 1e-8:
        raise ValueError("The builder's O–H vector must not be collinear with C16–O")
    first = projection / norm
    second = np.cross(axis, first)
    chains = set(protein_chains)
    protein = [a for a, _ in records if a.record == "ATOM" and a.chain in chains]
    if not protein:
        raise ValueError("No protein atoms selected for hydroxyl-H orientation sampling")
    excluded = {added_h_atom, added_o_atom, c16_atom}
    # Exclude H–O (1–2) and H–C16 (1–3); retain all other intraligand atoms,
    # including C16's remaining H (a 1–4 contact), as steric diagnostics.
    targets = protein + [a for a in ligand if a.name not in excluded]
    unknown = sorted({a.element for a in targets} - {"C", "N", "O", "S", "H"})
    if unknown:
        raise ValueError(f"Unsupported target radii in hydroxyl sampling: {unknown}")
    if any(not np.isfinite(a.xyz).all() for a in targets):
        raise ValueError("Nonfinite coordinates in hydroxyl-H target atoms")
    target_arrays = (
        np.array([atom.xyz for atom in targets]),
        np.array([1.2 if atom.element == "H" else BONDI_RADII_A[atom.element] for atom in targets]),
        np.array([_potential_acceptor_for_hydroxyl_h(atom, targets) for atom in targets], dtype=bool),
        np.array([atom.record == "ATOM" for atom in targets], dtype=bool),
    )
    angle = math.radians(109.5)
    orientations = []
    for index in range(n_orientations):
        torsion = 2 * math.pi * index / n_orientations
        direction = (math.cos(angle) * axis + math.sin(angle)
                     * (math.cos(torsion) * first + math.sin(torsion) * second))
        xyz = oxygen.xyz + .96 * direction
        score = _hydroxyl_h_contact_score(
            xyz, targets, overlap_warning_A=overlap_warning_A,
            short_polar_distance_A=short_polar_distance_A,
            target_arrays=target_arrays,
        )
        orientations.append({"orientation_index": index, "torsion_offset_degrees": 360*index/n_orientations,
                             "hydrogen_xyz_A": xyz.tolist(), **score})
    def objective(row):
        return (row["short_polar_contact_count"], row["short_polar_deficit_sum_A"],
                row["nonpolar_warning_atom_count"], row["maximum_nonpolar_overlap_A"],
                row["nonpolar_overlap_square_sum_A2"], row["orientation_index"])
    chosen = min(orientations, key=objective)
    persisted_xyz = np.round(np.asarray(chosen["hydrogen_xyz_A"]), 3)
    line = lines[hydrogen_line]
    lines[hydrogen_line] = (line[:30] + "".join(f"{value:8.3f}" for value in persisted_xyz) + line[54:])
    target = Path(output_pdb) if output_pdb is not None else source
    target.parent.mkdir(parents=True, exist_ok=True)
    target.write_bytes("".join(lines).encode("utf-8"))
    persisted_score = _hydroxyl_h_contact_score(
        persisted_xyz, targets, overlap_warning_A=overlap_warning_A,
        short_polar_distance_A=short_polar_distance_A,
        target_arrays=target_arrays,
    )
    return {
        "method": "new_16alpha_hydroxyl_h_torsion_steric_sampling_v1",
        "orientation_count": n_orientations,
        "selected_orientation_index": chosen["orientation_index"],
        "minimum_nonpolar_overlap_A": min(r["maximum_nonpolar_overlap_A"] for r in orientations),
        "minimum_protein_nonpolar_overlap_A": min(r["maximum_protein_nonpolar_overlap_A"] for r in orientations),
        "maximum_nonpolar_overlap_A": max(r["maximum_nonpolar_overlap_A"] for r in orientations),
        "selected_nonpolar_overlap_A": chosen["maximum_nonpolar_overlap_A"],
        "selected_persisted_nonpolar_overlap_A": persisted_score["maximum_nonpolar_overlap_A"],
        "nonpolar_blocked_orientation_count": sum(r["nonpolar_warning_atom_count"] > 0 for r in orientations),
        "short_polar_orientation_count": sum(r["short_polar_contact_count"] > 0 for r in orientations),
        "added_h_atom": added_h_atom, "added_h_xyz_A": persisted_xyz.tolist(),
        "selected_persisted_contacts": persisted_score["contacts"],
        "orientations": orientations,
        "protocol": {
            "scope": "new_hydroxyl_h_sterics_only_not_energy",
            "hydroxyl_O_H_bond_A": .96, "C_O_H_angle_degrees": 109.5,
            "added_polar_hydrogen_radius_A": 1.0, "target_hydrogen_radius_A": 1.2,
            "heavy_atom_radii_A": dict(BONDI_RADII_A),
            "warning_overlap_A": overlap_warning_A,
            "short_polar_distance_A": short_polar_distance_A,
            "short_polar_policy": "Potential acceptor contacts below this project diagnostic distance are prioritized for avoidance; other H–acceptor overlaps are not counted as nonpolar steric overlaps",
            "acceptor_policy": "O/S are possible acceptors; unprotonated HIS ring N is a possible acceptor; ordinary protein amide/Lys/Arg/Trp N are not",
            "excluded_ligand_atoms": sorted(excluded),
            "exclusion_policy": "O and C16 are bonded 1–2/1–3 exclusions; other ligand contacts including 1–4 are retained as geometry diagnostics",
            "selection": "Lexicographic: short-polar count, short-polar deficit, non-Hbond-compatible warning count, max overlap, sum squared overlaps, torsion index",
            "minimum_statistic": "Minimum across all sampled orientations of maximum non-Hbond-compatible overlap; not a continuum minimum or an energy",
            "summary_statistic": "The separately reported best sampled protein overlap excludes intraligand contacts; intraligand contacts only contribute to orientation selection and detailed diagnostics",
            "protein_atoms_included": "all selected protein heavy atoms and explicit hydrogens",
            "coordinates_moved": [added_h_atom], "minimization": False, "repacking": False,
            "energies_assessed": False, "hbond_detection": False,
            "radius_precedent": "Polar H 1.0 A follows MolProbity-era conventions, but the mixed Bondi/target-H screen is project-specific and is not PROBE",
            "reference": "https://doi.org/10.1093/nar/gkh398",
        },
        "provenance": {"input_pdb": str(source.resolve()), "input_sha256": sha256(original).hexdigest(),
                       "output_pdb": str(target.resolve()), "output_sha256": sha256(target.read_bytes()).hexdigest()},
    }
