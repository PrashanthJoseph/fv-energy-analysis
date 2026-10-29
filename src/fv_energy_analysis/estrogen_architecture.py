"""E2 recognition features, with Discovery Studio-defined A-ring pi contacts."""

from __future__ import annotations

import json
import csv
import re
from collections import defaultdict
from pathlib import Path

import numpy as np


BACKBONE_ATOMS = {"N", "CA", "C", "O", "OXT"}


def _ds_atoms(group: str) -> set[tuple[str, int, str, str, str]]:
    atoms = set()
    for atom in group.split(";"):
        match = re.fullmatch(r"([^:]+):([A-Za-z0-9]{3})(-?\d+)([A-Za-z]?):([^:]+)", atom.strip())
        if not match:
            raise ValueError(f"Cannot parse Discovery Studio atom: {atom!r}")
        chain, name, number, insertion, atom_name = match.groups()
        atoms.add((chain, int(number), insertion, name.upper(), atom_name))
    return atoms


def face_to_face_pi_geometry(ligand_ring: np.ndarray, partner_ring: np.ndarray) -> dict[str, object]:
    """Explicit project screen, not proof of a distinct pi-stacking energy.

    Require non-degenerate planes, centroid distance <=5.5 A, plane angle <=30
    degrees, and lateral centroid offset <=2 A on both ring planes.
    """
    if len(ligand_ring) < 3 or len(partner_ring) < 3:
        return {"passes": False, "reason": "insufficient_ring_atoms"}
    lc, pc = ligand_ring.mean(axis=0), partner_ring.mean(axis=0)
    _, ls, lv = np.linalg.svd(ligand_ring-lc, full_matrices=False)
    _, ps, pv = np.linalg.svd(partner_ring-pc, full_matrices=False)
    if ls[1] < 1e-6 or ps[1] < 1e-6:
        return {"passes": False, "reason": "degenerate_ring_plane"}
    ln, pn = lv[-1], pv[-1]
    delta = pc-lc
    distance = float(np.linalg.norm(delta))
    angle = float(np.degrees(np.arccos(np.clip(abs(np.dot(ln,pn)),0,1))))
    offset = max(float(np.linalg.norm(delta-np.dot(delta,normal)*normal)) for normal in (ln,pn))
    return {"passes": distance <= 5.5 and angle <= 30.0 and offset <= 2.0,
            "centroid_distance_A": distance, "plane_angle_degrees": angle, "lateral_offset_A": offset,
            "signed_face_offset_A": float(np.dot(delta,ln))}


def _endpoint_label(endpoint: dict[str, object]) -> str:
    chain = str(endpoint.get("auth_asym_id") or endpoint.get("label_asym_id") or "_")
    number = endpoint.get("auth_seq_id") or endpoint.get("label_seq_id")
    insertion = str(endpoint.get("pdbx_PDB_ins_code") or "").strip()
    if insertion in {"?", "."}:
        insertion = ""
    name = str(endpoint.get("label_comp_id") or "")
    return f"{chain}:{number}{insertion}:{name}"


def _read_pdb_atoms(path: str | Path) -> list[dict[str, object]]:
    atoms: list[dict[str, object]] = []
    for line in Path(path).read_text(encoding="utf-8", errors="replace").splitlines():
        if not line.startswith(("ATOM  ", "HETATM")):
            continue
        try:
            atoms.append({
                "name": line[12:16].strip(),
                "resname": line[17:20].strip(),
                "chain": line[21:22].strip() or "_",
                "number": int(line[22:26]),
                "insertion": line[26:27].strip(),
                "coord": np.array([
                    float(line[30:38]), float(line[38:46]), float(line[46:54])
                ]),
                "element": line[76:78].strip() or line[12:16].strip()[0],
            })
        except (ValueError, IndexError):
            continue
    return atoms


def _ligand_entry(entry: dict[str, object], ligand_resname: str) -> tuple[dict, dict] | None:
    bgn, end = entry.get("bgn"), entry.get("end")
    if not isinstance(bgn, dict) or not isinstance(end, dict):
        return None
    if str(bgn.get("label_comp_id")) == ligand_resname:
        return bgn, end
    if str(end.get("label_comp_id")) == ligand_resname:
        return end, bgn
    return None


def _analyze_legacy_estrogen_architecture(
    pdb_path: str | Path,
    arpeggio_json: str | Path,
    *,
    ligand_resname: str = "E2L",
    o3_atom: str | None = None,
    o17_atom: str | None = None,
    c18_methyl_atom: str | None = None,
    a_ring_atoms: tuple[str, ...] | None = None,
    fv_chains: set[str] | None = None,
    face_contact_cutoff_A: float = 4.5,
    minimum_face_offset_A: float = 0.5,
    params_path: str | Path | None = None,
    template_path: str | Path | None = None,
    discovery_studio_csv: str | Path | None = None,
    saved_features: dict | None = None,
) -> dict[str, object]:
    """Return literature-derived E2 recognition metrics for one structural pose."""
    if params_path is not None:
        from .atom_roles import load_atom_roles
        roles = load_atom_roles(params_path, template_path)
        if not roles:
            raise ValueError('Cannot determine E2 atom roles from ligand parameters; refusing a guessed atom-name mapping')
        o3_atom, o17_atom, c18_methyl_atom = roles['o3_atom'], roles['o17_atom'], roles['c18_methyl_atom']
        a_ring_atoms = tuple(roles['a_ring_atoms'])
    if not all((o3_atom, o17_atom, c18_methyl_atom, a_ring_atoms)):
        raise ValueError("Supply ligand params/template or all explicit E2 atom roles; no legacy atom-name defaults are used")
    data = [] if saved_features is not None else json.loads(Path(arpeggio_json).read_text(encoding="utf-8"))
    if isinstance(data, dict):
        data = next((data[key] for key in ("interactions", "contacts", "data") if isinstance(data.get(key), list)), [])
    if not isinstance(data, list):
        raise ValueError(f"Expected an Arpeggio interaction list in {arpeggio_json}")

    end_hbonds: dict[str, set[str]] = {o3_atom: set(), o17_atom: set()}
    c18_hydrophobic: set[str] = set()
    def saved_labels(key):
        value = (saved_features or {}).get(key, "none")
        return set() if value == "none" else set(str(value).split(";"))
    if saved_features is not None:
        end_hbonds[o3_atom] = saved_labels("e2_o3_hbond_residues")
        end_hbonds[o17_atom] = saved_labels("e2_o17_hbond_residues")
        c18_hydrophobic = saved_labels("e2_c18_hydrophobic_residues")
    pi_residues: dict[tuple[str, tuple[str, ...]], dict[str, object]] = {}
    ds_partners = {name: set() for name in ("pi_stack", "pi_t_stack", "pi_cation", "pi_anion")}
    ds_types = {"pi_pi_stacked": "pi_stack", "piPiStackedType": "pi_stack",
                "pi_pi_t_shaped": "pi_t_stack", "piPiTShapedType": "pi_t_stack",
                "pi_cation": "pi_cation", "piCationType": "pi_cation",
                "pi_anion": "pi_anion", "piAnionType": "pi_anion"}
    ligand_chain = ""
    ligand_number = 0
    ligand_insertion = ""

    for entry in data:
        if not isinstance(entry, dict):
            continue
        endpoints = _ligand_entry(entry, ligand_resname)
        if endpoints is None:
            continue
        ligand, partner = endpoints
        partner_chain = str(partner.get("auth_asym_id") or partner.get("label_asym_id") or "_")
        if fv_chains and partner_chain not in fv_chains:
            continue
        ligand_chain = str(ligand.get("auth_asym_id") or ligand.get("label_asym_id") or ligand_chain)
        ligand_number = int(ligand.get("auth_seq_id") or ligand.get("label_seq_id") or ligand_number)
        ligand_insertion = str(ligand.get("pdbx_PDB_ins_code") or "").strip()
        if ligand_insertion in {"?", "."}:
            ligand_insertion = ""
        ligand_atoms = {name.strip() for name in str(ligand.get("auth_atom_id") or "").split(",")}
        contacts = {str(value).lower() for value in entry.get("contact", [])}
        partner_label = _endpoint_label(partner)
        if "hbond" in contacts:
            for atom_name in (o3_atom, o17_atom):
                if atom_name in ligand_atoms:
                    end_hbonds[atom_name].add(partner_label)
        if c18_methyl_atom in ligand_atoms and "hydrophobic" in contacts:
            c18_hydrophobic.add(partner_label)

    ds_available = discovery_studio_csv is not None and Path(discovery_studio_csv).is_file()
    if ds_available:
        with Path(discovery_studio_csv).open(newline="", encoding="utf-8-sig") as handle:
            for row in csv.DictReader(handle):
                interaction = ds_types.get(row.get("interaction_type"))
                if interaction is None:
                    continue
                ligand_group = _ds_atoms(row["ligand_atoms"])
                if not set(a_ring_atoms).issubset({a[4] for a in ligand_group if a[3] == ligand_resname}):
                    continue
                ligand_ids = {a[:4] for a in ligand_group}
                if len(ligand_ids) != 1:
                    raise ValueError("Discovery Studio stack must identify one ligand residue")
                ds_chain, ds_number, ds_insertion, _ = next(iter(ligand_ids))
                if ligand_chain and (ligand_chain, ligand_number, ligand_insertion) != (ds_chain, ds_number, ds_insertion):
                    raise ValueError("Arpeggio and Discovery Studio identify different ligand residues")
                ligand_chain, ligand_number, ligand_insertion = ds_chain, ds_number, ds_insertion
                protein_group = _ds_atoms(row["protein_atoms"])
                identities = {a[:4] for a in protein_group}
                if len(identities) != 1:
                    raise ValueError("Discovery Studio pi-stack partner must identify one protein residue")
                chain, number, insertion, name = next(iter(identities))
                if fv_chains and chain not in fv_chains:
                    continue
                label = f"{chain}:{number}{insertion}:{name}"
                ds_partners[interaction].add(label)
                names = tuple(sorted(a[4] for a in protein_group))
                endpoint = pi_residues.setdefault((label, names), {
                    "auth_asym_id": chain, "auth_seq_id": number,
                    "pdbx_PDB_ins_code": insertion, "auth_atom_id": ",".join(names),
                    "label_comp_id": name, "contact_types": set(),
                })
                endpoint["contact_types"].add(interaction)

    atoms = _read_pdb_atoms(pdb_path)
    ligand_atoms_by_name = {
        str(atom["name"]): atom
        for atom in atoms
        if atom["resname"] == ligand_resname
        and (not ligand_chain or atom["chain"] == ligand_chain)
        and (not ligand_number or atom["number"] == ligand_number)
        and (not ligand_chain or atom["insertion"] == ligand_insertion)
    }
    missing = [name for name in a_ring_atoms if name not in ligand_atoms_by_name]
    if missing:
        raise ValueError(f"Missing E2 A-ring atoms in {pdb_path}: {', '.join(missing)}")
    ring_coords = np.array([ligand_atoms_by_name[name]["coord"] for name in a_ring_atoms])
    ring_centroid = ring_coords.mean(axis=0)
    _, _, vh = np.linalg.svd(ring_coords - ring_centroid, full_matrices=False)
    ring_normal = vh[-1]

    residues: dict[tuple[str, int, str, str], list[dict[str, object]]] = defaultdict(list)
    for atom in atoms:
        if atom["resname"] == ligand_resname or str(atom["element"]).upper() == "H":
            continue
        if fv_chains and atom["chain"] not in fv_chains:
            continue
        residues[(str(atom["chain"]), int(atom["number"]), str(atom["insertion"]), str(atom["resname"]))].append(atom)

    face_residues: dict[str, list[tuple[float, str]]] = {"positive": [], "negative": []}
    for (chain, number, insertion, resname), residue_atoms in residues.items():
        if saved_features is not None:
            break
        sidechain = [atom for atom in residue_atoms if atom["name"] not in BACKBONE_ATOMS]
        if not sidechain:
            continue
        distances = [
            min(float(np.linalg.norm(atom["coord"] - ring_atom)) for ring_atom in ring_coords)
            for atom in sidechain
        ]
        closest_index = int(np.argmin(distances))
        minimum_distance = distances[closest_index]
        if minimum_distance > face_contact_cutoff_A:
            continue
        signed_offset = float(np.dot(sidechain[closest_index]["coord"] - ring_centroid, ring_normal))
        if abs(signed_offset) < minimum_face_offset_A:
            continue
        label = f"{chain}:{number}{insertion}:{resname}"
        side = "positive" if signed_offset > 0 else "negative"
        face_residues[side].append((minimum_distance, label))

    for side in face_residues:
        if saved_features is not None:
            value = saved_features.get(f"e2_geometric_sandwich_{side}_face", "none")
            face_residues[side] = [] if value == "none" else list(enumerate(value.split(";")))
        face_residues[side].sort()
    geometric_sandwich = bool(face_residues["positive"] and face_residues["negative"])

    pi_offsets: dict[str, float] = {}
    stack_offsets: dict[str, float] = {}
    pi_geometries = {}
    for (label, ring_names), endpoint in pi_residues.items():
        chain = str(endpoint.get("auth_asym_id") or endpoint.get("label_asym_id") or "_")
        number = int(endpoint.get("auth_seq_id") or endpoint.get("label_seq_id"))
        insertion = str(endpoint.get("pdbx_PDB_ins_code") or "").strip()
        if insertion in {"?", "."}:
            insertion = ""
        names = {name.strip() for name in str(endpoint.get("auth_atom_id") or "").split(",")}
        coords = [
            atom["coord"] for atom in atoms
            if atom["chain"] == chain and atom["number"] == number
            and atom["insertion"] == insertion and atom["name"] in names
        ]
        if len(coords) != len(names):
            raise ValueError(f"Missing or duplicate coordinates for DS ring {label}: {sorted(names)}")
        geometry = face_to_face_pi_geometry(ring_coords, np.array(coords).reshape(-1,3))
        ring_label = f"{label}[{','.join(ring_names)}]"
        pi_geometries[ring_label] = geometry
        pi_offsets[ring_label] = float(np.dot(np.mean(coords, axis=0) - ring_centroid, ring_normal))
        if "pi_stack" in endpoint["contact_types"]:
            stack_offsets[ring_label] = pi_offsets[ring_label]

    def opposite_faces(offsets):
        items = list(offsets.items())
        return any(first * second < 0
                   and first_label.split("[")[0] != second_label.split("[")[0]
                   and abs(first) >= minimum_face_offset_A and abs(second) >= minimum_face_offset_A
                   for index, (first_label, first) in enumerate(items)
                   for second_label, second in items[index + 1:])
    pi_stack_sandwich = opposite_faces(stack_offsets)
    pi_contact_sandwich = opposite_faces(pi_offsets)

    o3_pass = bool(end_hbonds[o3_atom])
    o17_pass = bool(end_hbonds[o17_atom])
    c18_pass = bool(c18_hydrophobic)
    architecture_score = sum((o3_pass, o17_pass, c18_pass, geometric_sandwich))
    strict_pass = architecture_score == 4 and pi_contact_sandwich if ds_available else None
    return {
        "e2_architecture_status": "complete" if ds_available else "partial: Discovery Studio output missing",
        "e2_atom_role_source": "parameter chemical graph" if params_path else "explicit atom-name arguments",
        "e2_pi_stack_source": "Discovery Studio pi_pi_stacked" if ds_available else "unavailable: Discovery Studio output missing",
        "e2_pi_contact_source": "Discovery Studio: stacked, T-shaped, pi-cation, pi-anion" if ds_available else "unavailable: Discovery Studio output missing",
        "e2_a_ring_pi_contact_count": len({key[0] for key in pi_residues}) if ds_available else None,
        "e2_a_ring_pi_contact_residues": ";".join(sorted({key[0] for key in pi_residues})) or ("none" if ds_available else "unavailable"),
        "e2_pi_contact_sandwich": pi_contact_sandwich if ds_available else None,
        **{f"e2_ds_{kind}_partner_count": len(partners) if ds_available else None
           for kind, partners in ds_partners.items()},
        **{f"e2_ds_{kind}_residues": ";".join(sorted(partners)) or ("none" if ds_available else "unavailable")
           for kind, partners in ds_partners.items()},
        "e2_atom_roles": {"o3_atom": o3_atom, "o17_atom": o17_atom, "c18_methyl_atom": c18_methyl_atom, "a_ring_atoms": list(a_ring_atoms)},
        "e2_o3_hbond": o3_pass,
        "e2_o3_hbond_residues": ";".join(sorted(end_hbonds[o3_atom])) or "none",
        "e2_o17_hbond": o17_pass,
        "e2_o17_hbond_residues": ";".join(sorted(end_hbonds[o17_atom])) or "none",
        "e2_both_end_hbonds": o3_pass and o17_pass,
        "e2_c18_hydrophobic_contact": c18_pass,
        "e2_c18_hydrophobic_residues": ";".join(sorted(c18_hydrophobic)) or "none",
        "e2_a_ring_pi_stack_ring_count": len(stack_offsets) if ds_available else None,
        "e2_a_ring_pi_stack_count": len(ds_partners["pi_stack"]) if ds_available else None,
        "e2_a_ring_pi_stack_residues": ";".join(sorted(ds_partners["pi_stack"])) or ("none" if ds_available else "unavailable"),
        "e2_pi_stack_geometry": pi_geometries,
        "e2_pi_stack_geometry_criteria": "DS defines contact classes; geometry thresholds are diagnostic only. Project pi-contact sandwich includes stacked, T-shaped, pi-cation, pi-anion; requires distinct residues on opposite faces with partner-group centroid offsets >=0.5 A",
        "e2_pi_stack_sandwich": pi_stack_sandwich if ds_available else None,
        "e2_geometric_sandwich": geometric_sandwich,
        "e2_geometric_sandwich_positive_face": ";".join(label for _, label in face_residues["positive"][:3]) or "none",
        "e2_geometric_sandwich_negative_face": ";".join(label for _, label in face_residues["negative"][:3]) or "none",
        "e2_conserved_architecture_score_0_to_4": architecture_score,
        "e2_strict_screen_pass": strict_pass,
    }


def analyze_estrogen_architecture(
    pdb_path: str | Path,
    arpeggio_json: str | Path | None = None,
    *,
    ligand_resname: str = "E2L",
    o3_atom: str | None = None,
    o17_atom: str | None = None,
    c18_methyl_atom: str | None = None,
    a_ring_atoms: tuple[str, ...] | None = None,
    fv_chains: set[str] | None = None,
    params_path: str | Path | None = None,
    template_path: str | Path | None = None,
    discovery_studio_csv: str | Path | None = None,
    saved_features: dict | None = None,
    face_contact_cutoff_A: float = 4.5,
    minimum_face_offset_A: float = 0.5,
) -> dict[str, object]:
    """Current DS-only architecture, with a calibrated coordinate steric proxy.

    The old positional Arpeggio path and cached feature arguments are accepted
    for call-site compatibility but are NEVER read. Old broad distance/face
    settings do not control the two-sided A-ring enclosure definition. The private
    legacy implementation above exists only to reproduce historical reports.
    """
    from .ds_architecture import analyze_ds_architecture

    if params_path is not None:
        from .atom_roles import load_atom_roles
        roles = load_atom_roles(params_path, template_path)
    else:
        if not all((o3_atom, o17_atom, c18_methyl_atom, a_ring_atoms)):
            raise ValueError("Supply ligand params/template or all explicit E2 atom roles; no legacy atom-name defaults are used")
        roles = {"o3_atom": o3_atom, "o17_atom": o17_atom,
                 "c18_methyl_atom": c18_methyl_atom,
                 "a_ring_atoms": list(a_ring_atoms)}
    return analyze_ds_architecture(
        pdb_path, discovery_studio_csv, roles=roles,
        ligand_resname=ligand_resname, fv_chains=fv_chains,
    )
