"""DS-only interaction features plus the separate coordinate steric proxy.

No other interaction classifier is read. Native DS group annotations are not
atom-level energy measurements, and the reference-calibrated flags do not
establish affinity or estrogen-versus-androgen (or E2-versus-E3) selectivity.
"""
from __future__ import annotations

import csv
import hashlib
import json
import re
from collections import defaultdict
from pathlib import Path

import numpy as np

from .steric_sandwich import analyze_steric_sandwich


VERSION = "ds_e2_recognition_v3_enclosure"
PI_TYPES = {"pi_pi_stacked": "pi_stack", "pi_pi_t_shaped": "pi_t_stack",
            "pi_cation": "pi_cation", "pi_anion": "pi_anion"}
RAW_TYPES = {"conventionalHBondType": "conventional_hbond", "carbonHBondType": "carbon_hbond",
             "piPiStackedType": "pi_pi_stacked", "piPiTShapedType": "pi_pi_t_shaped",
             "piCationType": "pi_cation", "piAnionType": "pi_anion",
             "alkylType": "alkyl", "piAlkylType": "pi_alkyl"}
PROTEIN_HYDROGEN_ALIASES = {
    "ASN": {"HD21": "1HD2", "HD22": "2HD2"},
    "GLN": {"HE21": "1HE2", "HE22": "2HE2"},
    "ARG": {"HH11": "1HH1", "HH12": "2HH1", "HH21": "1HH2", "HH22": "2HH2"},
}


def _group(text: str) -> set[tuple]:
    atoms = set()
    for token in text.split(";"):
        match = re.fullmatch(r"([^:]+):([A-Za-z0-9]{3})(-?\d+)([A-Za-z]?):([^:]+)", token.strip())
        if not match:
            raise ValueError(f"Cannot parse Discovery Studio atom: {token!r}")
        chain, residue, number, insertion, name = match.groups()
        atoms.add((chain, int(number), insertion, residue.upper(), name))
    return atoms


def _label(identity: tuple) -> str:
    chain, number, insertion, residue = identity
    return f"{chain}:{number}{insertion}:{residue}"


def _resolve_protein_hydrogens(group: set[tuple], atoms: dict) -> tuple[set[tuple], list[dict]]:
    """Allow only enumerated receptor-H spelling aliases, never ligand/heavy renaming."""
    resolved, aliases = set(), []
    for key in sorted(group):
        actual = key
        if key not in atoms:
            alternate = PROTEIN_HYDROGEN_ALIASES.get(key[3], {}).get(key[4])
            candidate = (*key[:4], alternate) if alternate is not None else None
            if candidate in atoms and atoms[candidate][0] in {"H", "D"}:
                actual = candidate
                aliases.append({"discovery_studio_atom": f"{_label(key[:4])}:{key[4]}",
                                "actual_pdb_atom": f"{_label(actual[:4])}:{actual[4]}"})
        if actual in resolved:
            raise ValueError("Protein hydrogen alias collision: multiple DS atoms resolve to one PDB atom")
        resolved.add(actual)
    return resolved, aliases


def _atoms(path: Path) -> dict:
    """First model; highest-occupancy altloc, first wins ties, matching steric reader."""
    output = {}
    for line in path.read_text(encoding="utf-8").splitlines():
        if line.startswith("ENDMDL"):
            break
        if not line.startswith(("ATOM  ", "HETATM")):
            continue
        occupancy = float(line[54:60].strip() or "1")
        if occupancy <= 0:
            continue
        name = line[12:16].strip()
        key = (line[21:22].strip(), int(line[22:26]), line[26:27].strip(), line[17:20].strip(), name)
        element = line[76:78].strip().upper() or next((c.upper() for c in name if c.isalpha()), "")
        xyz = np.array([float(line[i:i + 8]) for i in (30, 38, 46)])
        if not np.isfinite(xyz).all() or not np.isfinite(occupancy):
            raise ValueError(f"Nonfinite PDB atom: {key}")
        if key not in output or occupancy > output[key][2]:
            output[key] = (element, xyz, occupancy)
    return output


def _hydroxyl_hydrogens(ligand: dict, oxygen: str, explicit) -> tuple[list[str], dict]:
    heavy = {name: atom for name, atom in ligand.items() if atom[0] not in {"H", "D"}}
    assigned = []
    distances = {}
    for name, (element, xyz, _) in ligand.items():
        if element not in {"H", "D"}:
            continue
        ranked = sorted((float(np.linalg.norm(xyz - atom[1])), target) for target, atom in heavy.items())
        distance, nearest = ranked[0]
        if nearest == oxygen and distance <= 1.25:
            if len(ranked) > 1 and abs(distance - ranked[1][0]) < 1e-6:
                raise ValueError(f"Ambiguous ligand hydrogen attachment: {name}")
            assigned.append(name)
            distances[name] = distance
    if explicit is not None and set(explicit) != set(assigned):
        raise ValueError(f"Explicit hydroxyl hydrogen mapping disagrees with coordinates for {oxygen}")
    if len(assigned) != 1:
        raise ValueError(f"Expected exactly one explicit hydroxyl hydrogen attached to {oxygen}; found {assigned}")
    return sorted(assigned), {"method": "unique nearest ligand heavy atom within 1.25 A",
                              "oxygen_atom": oxygen, "hydrogen_distances_A": distances}


def _load_calibration(value: dict | None) -> dict | None:
    if value is None:
        path = Path(__file__).with_name("data") / "c18_packing_calibration.json"
        value = json.loads(path.read_text(encoding="utf-8")) if path.is_file() else None
    if value is not None:
        rule = value.get("production_rule", {})
        if (set(rule.get("interaction_types", [])) != {"alkyl", "pi_alkyl"}
                or rule.get("ligand_endpoint") != "singleton_mapped_C18"
                or rule.get("minimum_unique_residues") != 1
                or rule.get("maximum_atom_distance_A") is not None):
            raise ValueError("Unsupported C18 calibration rule; expected singleton DS alkyl/pi-alkyl presence")
    return value


def analyze_ds_architecture(
    pdb_path: str | Path,
    discovery_studio_csv: str | Path | None,
    *,
    roles: dict,
    ligand_resname: str = "E2L",
    ligand_chain: str = "X",
    fv_chains: set[str] | None = None,
    ligand_resseq: int | None = None,
    c18_calibration: dict | None = None,
) -> dict:
    """Extract one pose's DS recognition calls without mixing classifier sources.

    Required roles: o3_atom, o17_atom, c18_methyl_atom, six a_ring_atoms.
    Optional o3_hydrogens/o17_hydrogens are checked against explicit coordinates.
    Missing DS data yields unavailable interaction flags, never negative calls.
    """
    pdb_path = Path(pdb_path)
    fv_chains = {"H", "L"} if fv_chains is None else set(fv_chains)
    if not fv_chains:
        raise ValueError("At least one protein chain is required")
    if roles.get("validated") is False:
        raise ValueError("Unvalidated E2 atom roles are not accepted")
    role_metadata = roles
    roles = roles.get("roles", roles)
    required = ("o3_atom", "o17_atom", "c18_methyl_atom", "a_ring_atoms")
    if any(not roles.get(key) for key in required):
        raise ValueError("Explicit chemically validated O3/O17/C18/A-ring atom roles are required")
    ring = tuple(roles["a_ring_atoms"])
    if len(ring) != 6 or len(set(ring)) != 6:
        raise ValueError("Exactly six distinct A-ring carbon names are required")
    atoms = _atoms(pdb_path)
    identities = {key[:4] for key in atoms if key[0] == ligand_chain and key[3] == ligand_resname
                  and (ligand_resseq is None or key[1] == ligand_resseq)}
    if len(identities) != 1:
        raise ValueError(f"Expected one selected ligand residue; found {len(identities)}")
    identity = next(iter(identities))
    ligand = {key[4]: atom for key, atom in atoms.items() if key[:4] == identity}
    for name in {roles["o3_atom"], roles["o17_atom"], roles["c18_methyl_atom"], *ring}:
        if name not in ligand:
            raise ValueError(f"Missing mapped ligand atom {name}")
    for role in ("o3_atom", "o17_atom"):
        if ligand[roles[role]][0] != "O":
            raise ValueError(f"Mapped {role} must be oxygen")
    if ligand[roles["c18_methyl_atom"]][0] != "C":
        raise ValueError("Mapped C18 methyl must be carbon")
    methyl_parent = roles.get("c18_parent_atom") or roles.get("c18_methyl_parent_atom")
    # A caller-supplied graph mapping or the validated E2L map orients the faces;
    # otherwise the pass is unchanged but alpha/beta have arbitrary normal sign.
    if methyl_parent is None and ligand_resname == "E2L" and roles["c18_methyl_atom"] == "C10" and "C16" in ligand:
        methyl_parent = "C16"
    steric = analyze_steric_sandwich(
        pdb_path, ligand_chain=ligand_chain, ligand_resname=ligand_resname,
        ring_atom_names=ring, protein_chains=sorted(fv_chains), ligand_resseq=identity[1],
        ligand_icode=identity[2], methyl_atom=roles["c18_methyl_atom"] if methyl_parent else None,
        methyl_parent_atom=methyl_parent, maximum_surface_gap_A=1.0,
    )
    ds_path = Path(discovery_studio_csv) if discovery_studio_csv is not None else None
    available = ds_path is not None and ds_path.is_file()
    calibration = _load_calibration(c18_calibration)
    matches = {key: [] for key in ("o3", "o17", "c18", "c18_group_diagnostics", "pi")}
    partners = defaultdict(set)
    hydrogen_mapping = {}
    protein_aliases = {}
    if available:
        hydroxyl = {}
        for role in ("o3", "o17"):
            hydrogens, audit = _hydroxyl_hydrogens(ligand, roles[f"{role}_atom"], roles.get(f"{role}_hydrogens"))
            hydroxyl[role] = {roles[f"{role}_atom"], *hydrogens}
            hydrogen_mapping[role] = audit
        with ds_path.open(newline="", encoding="utf-8-sig") as handle:
            reader = csv.DictReader(handle)
            if not {"interaction_type", "protein_atoms", "ligand_atoms"}.issubset(reader.fieldnames or []):
                raise ValueError("Normalized DS CSV must have interaction_type, protein_atoms and ligand_atoms")
            for row in reader:
                kind = RAW_TYPES.get(row["interaction_type"], row["interaction_type"])
                if kind not in {"conventional_hbond", "alkyl", "pi_alkyl", *PI_TYPES}:
                    continue
                ligand_group, protein_group = _group(row["ligand_atoms"]), _group(row["protein_atoms"])
                if {key[:4] for key in ligand_group} != {identity}:
                    raise ValueError("Discovery Studio CSV identifies a different ligand residue")
                protein_ids = {key[:4] for key in protein_group}
                if len(protein_ids) != 1:
                    raise ValueError("Discovery Studio endpoint must identify exactly one protein residue")
                protein_id = next(iter(protein_ids))
                if protein_id[0] not in fv_chains:
                    continue
                protein_group, aliases = _resolve_protein_hydrogens(protein_group, atoms)
                missing = (ligand_group | protein_group) - atoms.keys()
                if missing:
                    raise ValueError(f"Discovery Studio atoms missing from exact PDB: {sorted(missing)}")
                label = _label(protein_id)
                names = {key[4] for key in ligand_group}
                record = {**row, "interaction_type": kind, "protein_residue": label}
                if aliases:
                    if kind != "conventional_hbond" or len(protein_group) != 1 or len(ligand_group) != 1:
                        raise ValueError("Protein-H aliases require a single-atom conventional DS H-bond for distance verification")
                    protein_key, ligand_key = next(iter(protein_group)), next(iter(ligand_group))
                    actual_distance = float(np.linalg.norm(atoms[protein_key][1] - atoms[ligand_key][1]))
                    native_distance = float(row["distance_A"])
                    if not np.isfinite(native_distance) or abs(actual_distance - native_distance) > 0.002:
                        raise ValueError("Protein hydrogen alias distance disagrees with native DS H-bond distance")
                    record["protein_hydrogen_atom_aliases"] = aliases
                    record["protein_hydrogen_alias_distance_check"] = {
                        "actual_pdb_distance_A": actual_distance, "native_ds_distance_A": native_distance,
                        "maximum_absolute_difference_A": 0.002,
                    }
                    for alias in aliases:
                        protein_aliases[alias["discovery_studio_atom"]] = alias
                if kind == "conventional_hbond":
                    for role in ("o3", "o17"):
                        if names & hydroxyl[role]:
                            matches[role].append(record)
                            partners[role].add(label)
                if kind in {"alkyl", "pi_alkyl"} and roles["c18_methyl_atom"] in names:
                    methyl_xyz = ligand[roles["c18_methyl_atom"]][1]
                    distances = sorted((float(np.linalg.norm(methyl_xyz - atoms[key][1])), key)
                                       for key in protein_group if atoms[key][0] not in {"H", "D"})
                    if not distances:
                        raise ValueError("DS hydrophobic protein group has no heavy atoms")
                    carbon_sulfur = [item for item in distances if atoms[item[1]][0] in {"C", "S"}]
                    distance, closest = distances[0]
                    diagnostic = {**record, "singleton_mapped_c18_endpoint": names == {roles["c18_methyl_atom"]},
                                  "nearest_protein_heavy_distance_A": distance,
                                  "nearest_protein_heavy_atom": f"{_label(closest[:4])}:{closest[4]}",
                                  "nearest_protein_carbon_or_sulfur_distance_A": carbon_sulfur[0][0] if carbon_sulfur else None,
                                  "distance_policy": "Diagnostic only; DS native call is not recut by atom distance"}
                    matches["c18_group_diagnostics"].append(diagnostic)
                    if diagnostic["singleton_mapped_c18_endpoint"]:
                        matches["c18"].append(diagnostic)
                        partners["c18"].add(label)
                if kind in PI_TYPES and set(ring).issubset(names):
                    matches["pi"].append(record)
                    partners[PI_TYPES[kind]].add(label)
                    partners["pi"].add(label)
    c18_available = available and calibration is not None
    o3, o17 = (bool(partners[key]) if available else None for key in ("o3", "o17"))
    c18 = bool(partners["c18"]) if c18_available else None
    score = sum((o3, o17, c18, steric["steric_sandwich"])) if c18_available else None
    coverage = steric["coverage"]
    def labels(key, ready=available):
        return ";".join(sorted(partners[key])) or ("none" if ready else "unavailable")
    output = {
        "e2_architecture_status": "complete" if c18_available else "unavailable: Discovery Studio output missing" if not available else "partial: C18 calibration missing",
        "e2_architecture_definition_version": VERSION,
        "e2_atom_roles": dict(roles),
        "e2_atom_role_source": "explicit chemically validated mapping supplied by caller",
        "e2_interaction_source": "BIOVIA Discovery Studio native saved interaction classifications" if available else "unavailable",
        "e2_o3_hbond": o3, "e2_o17_hbond": o17,
        "e2_o3_hbond_residues": labels("o3"), "e2_o17_hbond_residues": labels("o17"),
        "e2_o3_hbond_annotation_count": len(matches["o3"]) if available else None,
        "e2_o17_hbond_annotation_count": len(matches["o17"]) if available else None,
        "e2_both_end_hbonds": o3 and o17 if available else None,
        "e2_c18_hydrophobic_contact": c18,
        "e2_c18_hydrophobic_partner_count": len(partners["c18"]) if c18_available else None,
        "e2_c18_hydrophobic_residues": labels("c18", c18_available),
        "e2_c18_definition": "Native DS alkyl/pi-alkyl with ligand endpoint exactly the mapped C18 methyl carbon; >=1 unique protein residue; no additional atom-distance cutoff",
        "e2_c18_calibration": calibration,
        "e2_a_ring_pi_contact_count": len(partners["pi"]) if available else None,
        "e2_a_ring_pi_contact_residues": labels("pi"),
        "e2_pi_contact_source": "DS A-ring pi-mediated contacts: stacked, T-shaped, pi-cation, pi-anion; unique residues",
        "e2_a_ring_pi_stack_count": len(partners["pi_stack"]) if available else None,
        "e2_a_ring_pi_stack_ring_count": sum(row["interaction_type"] == "pi_pi_stacked" for row in matches["pi"]) if available else None,
        "e2_a_ring_pi_stack_residues": labels("pi_stack"),
        "e2_geometric_sandwich": steric["steric_sandwich"],
        "e2_steric_method": steric["method"],
        "e2_steric_alpha_residue_count": coverage["alpha_residue_count"],
        "e2_steric_beta_residue_count": coverage["beta_residue_count"],
        "e2_steric_opposing_residue_pair_count": coverage["opposing_residue_pair_count"],
        "e2_steric_overlapping_atom_pair_count": coverage["overlapping_atom_pair_count"],
        "e2_steric_sandwich": steric,
        "e2_conserved_architecture_score_0_to_4": score,
        "e2_strict_screen_pass": score == 4 if score is not None else None,
        "e2_feature_sources": {"o3_hbond": "DS conventional_hbond only", "o17_hbond": "DS conventional_hbond only",
                               "c18_hydrophobic": "DS singleton C18 alkyl/pi-alkyl calibrated presence",
                               "steric_sandwich": "Separate project-defined paper-reference-calibrated coordinate geometry, not native DS",
                               "pi_mediated": "DS optional descriptor, not required for four-feature pass"},
        "e2_matched_ds_records": matches,
        "e2_hydroxyl_hydrogen_mapping": hydrogen_mapping,
        "e2_provenance": {"pdb": str(pdb_path), "pdb_sha256": hashlib.sha256(pdb_path.read_bytes()).hexdigest(),
                          "ds_csv": str(ds_path) if ds_path else None,
                          "ds_csv_sha256": hashlib.sha256(ds_path.read_bytes()).hexdigest() if available else None,
                          "mapping_params_sha256": role_metadata.get("params_sha256"),
                          "protein_hydrogen_atom_aliases": sorted(protein_aliases.values(), key=lambda item: item["discovery_studio_atom"]),
                          "protein_hydrogen_alias_policy": "Exact names preferred; enumerated ASN/GLN/ARG receptor H aliases only; native DS distance checked within 0.002 A; coordinates unchanged",
                          "reference_paper_doi": "10.1096/fj.02-0524fje"},
    }
    for kind in PI_TYPES.values():
        output[f"e2_ds_{kind}_partner_count"] = len(partners[kind]) if available else None
        output[f"e2_ds_{kind}_residues"] = labels(kind)
    for face, legacy in (("beta", "positive"), ("alpha", "negative")):
        residues = {contact["wall_residue"] for contact in steric["contacts"] if contact["face"] == face}
        output[f"e2_geometric_sandwich_{legacy}_face"] = ";".join(sorted(residues)) or "none"
    return output
