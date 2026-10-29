"""Explicit donor/acceptor and hydrogen geometry for E3 overlap exemptions.

This is a geometric compatibility screen, not a hydrogen-bond energy model.
Only one supplied E3 hydroxyl-H orientation may be used for a joint assessment.
"""
from __future__ import annotations

import math
from typing import Sequence

import numpy as np

from .steric_sandwich import StericAtom


METHOD_VERSION = "e3_explicit_hbond_overlap_compatibility_v1"
MAX_DONOR_ACCEPTOR_A = 4.1
MIN_DONOR_H_ACCEPTOR_DEG = 100.0
MIN_H_ACCEPTOR_A = 1.2
PLIP_CONFIG_REFERENCE = "https://github.com/pharmai/plip/blob/master/plip/basic/config.py"
PLIP_DETECTION_REFERENCE = "https://github.com/pharmai/plip/blob/master/plip/structure/detection.py"

_STANDARD = {"ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE",
             "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL",
             "HID", "HIE", "HIP", "CYX", "CYM", "ASH", "GLH", "LYN"}
_BACKBONE_H = {"H", "HN", "H1", "H2", "H3", "1H", "2H", "3H", "HT1", "HT2", "HT3"}
_DONORS = {
    ("ARG", "NE"): {"HE"},
    ("ARG", "NH1"): {"HH11", "HH12", "1HH1", "2HH1"},
    ("ARG", "NH2"): {"HH21", "HH22", "1HH2", "2HH2"},
    ("ASN", "ND2"): {"HD21", "HD22", "1HD2", "2HD2"},
    ("GLN", "NE2"): {"HE21", "HE22", "1HE2", "2HE2"},
    ("LYS", "NZ"): {"HZ", "HZ1", "HZ2", "HZ3", "1HZ", "2HZ", "3HZ"},
    ("LYN", "NZ"): {"HZ", "HZ1", "HZ2", "HZ3", "1HZ", "2HZ", "3HZ"},
    ("TRP", "NE1"): {"HE1"},
    ("SER", "OG"): {"HG"}, ("THR", "OG1"): {"HG1"}, ("TYR", "OH"): {"HH"},
    ("CYS", "SG"): {"HG"},
}


def _attached(atom, atoms, names):
    """Match standard bonded-H names and covalent distance in this residue."""
    maximum = 1.6 if atom.element == "S" else 1.30
    return [h for h in atoms if h.residue_id == atom.residue_id and h.element in {"H", "D"}
            and h.name in names and .55 <= float(np.linalg.norm(h.xyz-atom.xyz)) <= maximum]


def _chemistry(atom, atoms):
    res, name = atom.residue, atom.name
    result = {"recognized": True, "donor": False, "acceptor": False,
              "hydrogens": [], "missing_donor_hydrogen": False, "typing": ""}
    if res not in _STANDARD:
        return {**result, "recognized": False, "typing": f"unsupported residue {res}"}
    if name in {"O", "OXT"} and atom.element == "O":
        hs = _attached(atom, atoms, {"HXT", "HO", "HOXT"})
        return {**result, "donor": bool(hs), "acceptor": not bool(hs), "hydrogens": hs,
                "typing": "backbone carbonyl/terminal oxygen; protonated OH not an acceptor"}
    if name == "N" and atom.element == "N":
        previous_c = any(a.chain == atom.chain and a.residue_id != atom.residue_id
                         and a.name == "C" and a.element == "C"
                         and float(np.linalg.norm(a.xyz-atom.xyz)) < 1.8 for a in atoms)
        hs = _attached(atom, atoms, _BACKBONE_H)
        if res == "PRO" and previous_c:
            return {**result, "typing": "internal proline amide N is neither donor nor acceptor"}
        return {**result, "donor": True, "hydrogens": hs, "missing_donor_hydrogen": not bool(hs),
                "typing": "backbone/terminal amino N donates but does not accept"}
    if res in {"HIS", "HID", "HIE", "HIP"} and name in {"ND1", "NE2"}:
        hs = _attached(atom, atoms, {"HD1"} if name == "ND1" else {"HE2"})
        ring = [a for a in atoms if a.residue_id == atom.residue_id and a.name in {"ND1", "NE2"}]
        all_ring_h = [h for n in ring for h in _attached(n, atoms, {"HD1"} if n.name == "ND1" else {"HE2"})]
        if not all_ring_h:
            return {**result, "recognized": False, "typing": "histidine protonation lacks explicit ring H"}
        return {**result, "donor": bool(hs), "acceptor": not bool(hs), "hydrogens": hs,
                "typing": "histidine ring N typed from its explicit bonded hydrogen"}
    if (res, name) in _DONORS:
        hs = _attached(atom, atoms, _DONORS[(res, name)])
        hydroxyl = atom.element in {"O", "S"}
        if res == "CYS" and name == "SG" and not hs and any(
            a.residue_id != atom.residue_id and a.element == "S" and a.name == "SG"
            and float(np.linalg.norm(a.xyz-atom.xyz)) <= 2.3 for a in atoms
        ):
            return {**result, "acceptor": True, "typing": "disulfide sulfur is not a donor; possible sulfur acceptor"}
        return {**result, "donor": True, "acceptor": hydroxyl and res not in {"ASH", "GLH"},
                "hydrogens": hs, "missing_donor_hydrogen": not bool(hs),
                "typing": "hydroxyl/thiol donor and possible acceptor" if hydroxyl else "standard side-chain donor N; not an acceptor"}
    if res in {"ASP", "ASH"} and name in {"OD1", "OD2"}:
        hs = _attached(atom, atoms, {"HD1", "HD2"})
        return {**result, "donor": bool(hs), "acceptor": not bool(hs), "hydrogens": hs,
                "typing": "carboxyl oxygen typed from explicit protonation"}
    if res in {"GLU", "GLH"} and name in {"OE1", "OE2"}:
        hs = _attached(atom, atoms, {"HE1", "HE2"})
        return {**result, "donor": bool(hs), "acceptor": not bool(hs), "hydrogens": hs,
                "typing": "carboxyl oxygen typed from explicit protonation"}
    if (res, name) in {("ASN", "OD1"), ("GLN", "OE1"), ("MET", "SD"), ("CYM", "SG"), ("CYX", "SG")}:
        return {**result, "acceptor": True, "typing": "standard carbonyl or sulfur acceptor"}
    if atom.element == "C":
        return {**result, "typing": "carbon is not a conventional hydrogen-bond donor/acceptor"}
    return {**result, "recognized": False, "typing": f"unsupported standard-residue atom {res}:{name}"}


def _trial(donor, hydrogen, acceptor, dxyz, hxyz, axyz, orientation_index=None):
    da = float(np.linalg.norm(dxyz-axyz))
    ha = float(np.linalg.norm(hxyz-axyz))
    dh = dxyz-hxyz
    ah = axyz-hxyz
    length = float(np.linalg.norm(dh)*np.linalg.norm(ah))
    angle = math.degrees(math.acos(float(np.clip(np.dot(dh, ah)/length, -1, 1)))) if length > 1e-12 else None
    criteria = {"donor_acceptor_distance_pass": da < MAX_DONOR_ACCEPTOR_A,
                "hydrogen_acceptor_distance_pass": ha >= MIN_H_ACCEPTOR_A,
                "angle_pass": angle is not None and angle > MIN_DONOR_H_ACCEPTOR_DEG}
    return {"donor": donor, "hydrogen": hydrogen, "acceptor": acceptor,
            "orientation_index": orientation_index, "donor_acceptor_distance_A": da,
            "hydrogen_acceptor_distance_A": ha, "angle_DHA_deg": angle,
            "passes_geometry": bool(all(criteria.values())), "criteria": criteria}


def classify_polar_overlap(protein_atom: StericAtom, protein_all_atoms: Sequence[StericAtom],
                           added_o_xyz_A, added_h_positions_A=None) -> dict:
    """Check one polar contact for supported conventional H-bond geometry.

    The caller evaluates a single joint OH-H orientation at a time. A positions
    list is accepted for API compatibility and each trial retains its index;
    callers must not pool incompatible orientations into one contact network.
    """
    atoms = list(protein_all_atoms)
    oxygen = np.asarray(added_o_xyz_A, dtype=float)
    if oxygen.shape != (3,) or not np.isfinite(oxygen).all():
        raise ValueError("Added oxygen requires three finite coordinates")
    chemistry = _chemistry(protein_atom, atoms)
    trials = []
    if chemistry["donor"]:
        for h in chemistry["hydrogens"]:
            trial = _trial(protein_atom.label, h.label, "E3:O16alpha", protein_atom.xyz, h.xyz, oxygen)
            trial["direction"] = "protein_donor_to_E3_acceptor"
            trials.append(trial)
    positions = [] if added_h_positions_A is None else list(added_h_positions_A)
    if chemistry["acceptor"]:
        for index, position in enumerate(positions):
            xyz = np.asarray(position, dtype=float)
            if xyz.shape != (3,) or not np.isfinite(xyz).all():
                raise ValueError("Added hydroxyl hydrogen requires three finite coordinates")
            oh = float(np.linalg.norm(xyz-oxygen))
            if not .70 <= oh <= 1.20:
                raise ValueError(f"Added hydroxyl H is not bonded to the added oxygen: O-H={oh:.4f} A")
            trial = _trial("E3:O16alpha", "E3:H16alphaOH", protein_atom.label,
                           oxygen, xyz, protein_atom.xyz, index)
            trial["direction"] = "E3_donor_to_protein_acceptor"
            trials.append(trial)
    if any(t["passes_geometry"] for t in trials):
        status, reason = "supported", "Compatible donor/acceptor chemistry and explicit D-H-A geometry pass all gates"
    elif not chemistry["recognized"]:
        status, reason = "unassessed", chemistry["typing"]
    elif chemistry["missing_donor_hydrogen"] or chemistry["acceptor"] and not positions:
        status, reason = "unassessed", "Required donor hydrogen coordinates are missing; failed known directions cannot exclude the missing direction"
    else:
        status, reason = "not_supported", "No chemically compatible direction passes the explicit hydrogen geometry gates"
    return {"status": status, "state": status, "supported": status == "supported", "reason": reason,
            "chemistry": {k: v for k, v in chemistry.items() if k != "hydrogens"},
            "trials": trials, "method": METHOD_VERSION,
            "thresholds": {"donor_acceptor_distance_A": {"operator": "<", "value": MAX_DONOR_ACCEPTOR_A},
                           "angle_DHA_deg": {"operator": ">", "value": MIN_DONOR_H_ACCEPTOR_DEG},
                           "hydrogen_acceptor_distance_A": {"operator": ">=", "value": MIN_H_ACCEPTOR_A},
                           "source": "PLIP donor-angle/donor-acceptor gates; H-A lower bound is a project short-contact guard"},
            "references": [PLIP_CONFIG_REFERENCE, PLIP_DETECTION_REFERENCE]}
