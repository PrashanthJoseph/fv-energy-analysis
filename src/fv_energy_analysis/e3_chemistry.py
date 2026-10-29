"""Construct a stereochemically checked 16alpha-OH steric probe of bound E2.

This is a counterfactual geometry test, not an estriol docking or relaxation.
The twenty E2 heavy-atom coordinates are retained exactly. A C16 hydrogen is
replaced by oxygen along its existing C-H direction at a model C-O length of
1.43 A. Both C16 directions are tested; only the authoritative estriol identity
is accepted. The neutral E3 probe retains the other 23 ligand hydrogens and
adds one explicit hydroxyl hydrogen. These are starting coordinates, not a
Rosetta hydrogen optimization or an energetic calculation.
"""
from __future__ import annotations

from hashlib import sha256
from functools import lru_cache
from pathlib import Path
import math

import numpy as np

from .ligand_chemistry import ESTRADIOL_INCHIKEY, _params_molecule


ESTRIOL_SMILES = "C[C@]12CC[C@H]3[C@H]([C@@H]1C[C@H]([C@@H]2O)O)CCC4=C3C=CC(=C4)O"
ESTRIOL_INCHIKEY = "PROQIPRRNZUXQM-ZXXIGWHRSA-N"
ESTRIOL_FORMULA = "C18H24O3"
REFERENCE = {
    "name": "Estriol (16alpha,17beta)",
    "pubchem_cid": 5756,
    "url": "https://pubchem.ncbi.nlm.nih.gov/compound/5756",
    "inchikey": ESTRIOL_INCHIKEY,
    "smiles": ESTRIOL_SMILES,
    "chemical_difference": "An additional 16alpha hydroxyl; C18 remains methyl",
}
METHOD = "stereovalidated_c16_hydrogen_direction_full_h_probe_v2"


@lru_cache(maxsize=8)
def _template_topology(params, template, params_digest, template_digest):
    # The file digests bind the cached topology to its actual chemical source.
    return _params_molecule(params, template)[0]


def _name(atom):
    info = atom.GetPDBResidueInfo()
    if info is None:
        raise ValueError("Named E2 template atoms are required")
    return info.GetName().strip()


def _identity_from_geometry(molecule):
    """Use actual 3D geometry, rather than inherited stereochemical tags."""
    from rdkit import Chem

    checked = Chem.Mol(molecule)
    Chem.SanitizeMol(checked)
    Chem.AssignAtomChiralTagsFromStructure(checked, replaceExistingTags=True)
    Chem.AssignStereochemistry(checked, cleanIt=True, force=True)
    return checked, Chem.MolToInchiKey(checked)


def _c16_site(molecule):
    """Locate steroid C16 from graph chemistry, not PDB atom numbering.

    E2's unique aliphatic OH identifies C17. In the carbon-only five-membered
    D ring, C17 neighbors quaternary C13 and methylene C16. C13 is the neighbor
    bonded to the unique angular methyl. This also audits the expected topology.
    """
    carbon_rings = [set(ring) for ring in molecule.GetRingInfo().AtomRings()
                    if len(ring) == 5 and all(molecule.GetAtomWithIdx(i).GetAtomicNum() == 6 for i in ring)]
    oxygens = [atom for atom in molecule.GetAtoms() if atom.GetAtomicNum() == 8
               and any(n.GetAtomicNum() == 6 and not n.GetIsAromatic() for n in atom.GetNeighbors())]
    methyls = [atom for atom in molecule.GetAtoms() if atom.GetAtomicNum() == 6
               and sum(n.GetAtomicNum() == 1 for n in atom.GetNeighbors()) == 3]
    if len(carbon_rings) != 1 or len(oxygens) != 1 or len(methyls) != 1:
        raise ValueError("E2 graph does not uniquely identify the D ring, O17 and C18 methyl")
    c17 = next(n for n in oxygens[0].GetNeighbors() if n.GetAtomicNum() == 6)
    c13 = next(n for n in methyls[0].GetNeighbors() if n.GetAtomicNum() == 6)
    d_ring = carbon_rings[0]
    if c17.GetIdx() not in d_ring or c13.GetIdx() not in d_ring:
        raise ValueError("O17 and angular methyl do not identify the same steroid D ring")
    candidates = [n for n in c17.GetNeighbors() if n.GetIdx() in d_ring and n.GetIdx() != c13.GetIdx()]
    if len(candidates) != 1:
        raise ValueError("Cannot identify one C16 methylene beside C17")
    c16 = candidates[0]
    hydrogens = [n.GetIdx() for n in c16.GetNeighbors() if n.GetAtomicNum() == 1]
    if len(hydrogens) != 2:
        raise ValueError("Two explicit C16 hydrogens are required for the directional probe")
    return c16.GetIdx(), c17.GetIdx(), c13.GetIdx(), hydrogens


def _graft_e3(molecule, bond_length_A=1.43):
    """Return the one C16-H replacement giving the canonical estriol identity."""
    from rdkit import Chem
    from rdkit.Chem import rdMolDescriptors

    if not math.isfinite(bond_length_A) or bond_length_A <= 0:
        raise ValueError("C-O model bond length must be finite and positive")
    molecule, identity = _identity_from_geometry(molecule)
    if identity != ESTRADIOL_INCHIKEY:
        raise ValueError(f"Source geometry is not canonical E2: {identity}")
    c16, c17, c13, hydrogen_indices = _c16_site(molecule)
    xyz = molecule.GetConformer().GetPositions()
    candidates, diagnostics = [], []
    for replace_index in hydrogen_indices:
        direction = xyz[replace_index] - xyz[c16]
        norm = float(np.linalg.norm(direction))
        if norm < 1e-8:
            raise ValueError("Degenerate C16 C-H direction")
        direction /= norm
        oxygen_xyz = xyz[c16] + bond_length_A * direction
        builder = Chem.RWMol(molecule)
        oxygen = Chem.Atom(8)
        oxygen.SetNoImplicit(True)
        info = Chem.AtomPDBResidueInfo()
        info.SetName(" O3 ")
        info.SetResidueName("E3L")
        info.SetIsHeteroAtom(True)
        oxygen.SetMonomerInfo(info)
        builder.ReplaceAtom(replace_index, oxygen, preserveProps=False)
        hydroxyl_h = Chem.Atom(1)
        h_info = Chem.AtomPDBResidueInfo()
        h_info.SetName(" H25")
        h_info.SetResidueName("E3L")
        h_info.SetIsHeteroAtom(True)
        hydroxyl_h.SetMonomerInfo(h_info)
        added_h = builder.AddAtom(hydroxyl_h)
        builder.AddBond(replace_index, added_h, Chem.BondType.SINGLE)
        trial = builder.GetMol()
        conf = trial.GetConformer()
        conf.SetAtomPosition(replace_index, oxygen_xyz)
        # This is only a deterministic starting OH torsion. A geometry-only
        # sampler may replace this coordinate before the all-atom contact test.
        c_direction = -direction
        tangent = xyz[c17] - xyz[c16]
        tangent -= np.dot(tangent, c_direction) * c_direction
        tangent /= np.linalg.norm(tangent)
        angle = math.radians(109.5)
        oh_direction = math.cos(angle) * c_direction + math.sin(angle) * tangent
        conf.SetAtomPosition(added_h, oxygen_xyz + .96 * oh_direction)
        checked, key = _identity_from_geometry(trial)
        diagnostics.append({"replaced_h_atom": _name(molecule.GetAtomWithIdx(replace_index)), "inchikey": key})
        if key == ESTRIOL_INCHIKEY:
            if rdMolDescriptors.CalcMolFormula(checked) != ESTRIOL_FORMULA:
                raise ValueError("Estriol graft has an unexpected formula")
            candidates.append((checked, replace_index, added_h))
    if len(candidates) != 1:
        raise ValueError(f"Expected exactly one stereochemically valid 16alpha-OH direction; found {len(candidates)}: {diagnostics}")
    checked, added_o, added_h = candidates[0]
    return checked, {"c16_index": c16, "c17_index": c17, "c13_index": c13,
                     "added_o_index": added_o, "added_h_index": added_h,
                     "replaced_h_atom": _name(molecule.GetAtomWithIdx(added_o)),
                     "tested_directions": diagnostics}


def build_e3_probe(
    pdb_path: str | Path,
    e2_params: str | Path,
    e2_template: str | Path,
    output_pdb: str | Path | None = None,
    *,
    ligand_chain="X",
    ligand_resname="E2L",
    new_resname="E3L",
    carbon_oxygen_bond_A=1.43,
    added_oh_h_xyz_A=None,
) -> dict:
    """Add correctly directed estriol oxygen to an explicitly protonated E2 pose.

    The params/template supply authoritative bond orders and chemical states;
    actual complex coordinates determine the stereochemistry and graft direction.
    No conformer embedding, minimization, scoring or hydrogen optimization occurs.
    """
    from rdkit import Chem
    from rdkit.Chem import rdMolDescriptors

    pdb_path = Path(pdb_path)
    text = pdb_path.read_text(encoding="utf-8")
    original_sha = sha256(pdb_path.read_bytes()).hexdigest()
    params_digest = sha256(Path(e2_params).read_bytes()).hexdigest()
    template_digest = sha256(Path(e2_template).read_bytes()).hexdigest()
    molecule = Chem.Mol(_template_topology(str(Path(e2_params).resolve()), str(Path(e2_template).resolve()),
                                        params_digest, template_digest))
    template_names = {_name(atom): atom.GetIdx() for atom in molecule.GetAtoms()}
    if molecule.GetNumHeavyAtoms() != 20 or len(template_names) != 44:
        raise ValueError("Expected full canonical E2 template with twenty heavy atoms and twenty-four H")
    if {"O3", "H25"} & template_names.keys():
        raise ValueError("Probe names O3/H25 conflict with source E2 atom names")
    ligand_lines, residues = {}, set()
    in_first_model = []
    for line in text.splitlines():
        if line.startswith("ENDMDL"):
            break
        in_first_model.append(line)
        if line.startswith(("ATOM  ", "HETATM")) and line[21:22] == ligand_chain and line[17:20].strip() == ligand_resname:
            name = line[12:16].strip()
            if name in ligand_lines:
                raise ValueError(f"Duplicate ligand atom or unhandled altloc: {name}")
            if line[16:17].strip():
                raise ValueError("Resolve ligand alternative locations before constructing the E3 probe")
            residues.add(line[21:27])
            ligand_lines[name] = line
    if len(residues) != 1 or set(ligand_lines) != set(template_names):
        missing = sorted(set(template_names)-set(ligand_lines))
        extra = sorted(set(ligand_lines)-set(template_names))
        raise ValueError(f"Expected one complete explicit-H E2 residue; missing={missing}, extra={extra}, residues={len(residues)}")
    conformer = molecule.GetConformer()
    for name, index in template_names.items():
        line = ligand_lines[name]
        xyz = np.array([float(line[start:start+8]) for start in (30, 38, 46)])
        if not np.isfinite(xyz).all():
            raise ValueError(f"Invalid ligand coordinate: {name}")
        conformer.SetAtomPosition(index, xyz)
    graft, site = _graft_e3(molecule, carbon_oxygen_bond_A)
    added_xyz = np.array(graft.GetConformer().GetAtomPosition(site["added_o_index"]))
    c16_xyz = np.array(graft.GetConformer().GetAtomPosition(site["c16_index"]))
    if added_oh_h_xyz_A is not None:
        override = np.asarray(added_oh_h_xyz_A, dtype=float)
        if override.shape != (3,) or not np.isfinite(override).all():
            raise ValueError("Added hydroxyl H override requires three finite coordinates")
        oh = override - added_xyz
        oc = c16_xyz - added_xyz
        distance = float(np.linalg.norm(oh))
        if abs(distance - .96) > .005:
            raise ValueError("Added hydroxyl H override must retain the model 0.96 A O-H bond")
        angle = math.degrees(math.acos(float(np.clip(np.dot(oh, oc)/(distance*np.linalg.norm(oc)), -1, 1))))
        if abs(angle - 109.5) > .25:
            raise ValueError("Added hydroxyl H override must retain the model 109.5 degree C-O-H angle")
        graft.GetConformer().SetAtomPosition(site["added_h_index"], override)
    added_h_xyz = np.array(graft.GetConformer().GetAtomPosition(site["added_h_index"]))
    heavy_names = [_name(atom) for atom in molecule.GetAtoms() if atom.GetAtomicNum() > 1]
    for atom in molecule.GetAtoms():
        if atom.GetAtomicNum() > 1 and not np.array_equal(
            molecule.GetConformer().GetPositions()[atom.GetIdx()], graft.GetConformer().GetPositions()[atom.GetIdx()]
        ):
            raise ValueError("Shared steroid heavy-atom geometry changed during graft")
    exported = None
    if output_pdb is not None:
        output_pdb = Path(output_pdb)
        if output_pdb.resolve() == pdb_path.resolve():
            raise ValueError("Probe output must differ from the source structure")
        output_pdb.parent.mkdir(parents=True, exist_ok=True)
        atom_lines = [line for line in in_first_model if line.startswith(("ATOM  ", "HETATM"))]
        serial = max(int(line[6:11]) for line in atom_lines) + 1
        reference_line = ligand_lines[heavy_names[0]]
        residue_field = reference_line[21:27]
        output = ["REMARK 950 COUNTERFACTUAL NEUTRAL E3 16ALPHA-OH PROBE; NOT DOCKED OR RELAXED",
                  "REMARK 950 SHARED ATOMS FIXED; FULL LIGAND H; NEW OH H TORSION IS MODELLED"]
        for line in atom_lines:
            if line[21:22] == ligand_chain and line[17:20].strip() == ligand_resname:
                if line[12:16].strip() != site["replaced_h_atom"]:
                    output.append(line[:17] + f"{new_resname:>3}" + line[20:])
            else:
                output.append(line)
        x, y, z = added_xyz
        output.append(f"HETATM{serial:5d}  O3  {new_resname:>3} {residue_field}   "
                      f"{x:8.3f}{y:8.3f}{z:8.3f}{1.0:6.2f}{0.0:6.2f}           O")
        x, y, z = added_h_xyz
        output.append(f"HETATM{serial+1:5d}  H25 {new_resname:>3} {residue_field}   "
                      f"{x:8.3f}{y:8.3f}{z:8.3f}{1.0:6.2f}{0.0:6.2f}           H")
        output.append("END")
        output_pdb.write_text("\n".join(output) + "\n", encoding="utf-8")
        exported = str(output_pdb.resolve())
    if sha256(pdb_path.read_bytes()).hexdigest() != original_sha:
        raise ValueError("Source structure changed during probe construction")
    return {
        "method": METHOD, "status": "complete", "source_pdb": str(pdb_path.resolve()),
        "source_pdb_sha256": original_sha, "e3_pdb": exported,
        "e2_params_sha256": params_digest, "e2_template_sha256": template_digest,
        "e3_pdb_sha256": sha256(Path(exported).read_bytes()).hexdigest() if exported else None,
        "reference": REFERENCE, "formula": rdMolDescriptors.CalcMolFormula(graft),
        "inchikey": Chem.MolToInchiKey(graft), "identity_validated": True,
        "source_e2_inchikey": ESTRADIOL_INCHIKEY,
        "c16_atom": _name(molecule.GetAtomWithIdx(site["c16_index"])),
        "c16_xyz_A": c16_xyz.tolist(),
        "c17_atom": _name(molecule.GetAtomWithIdx(site["c17_index"])),
        "c13_atom": _name(molecule.GetAtomWithIdx(site["c13_index"])),
        "added_o_atom": "O3", "added_o_xyz_A": added_xyz.tolist(),
        "added_h_atom": "H25", "added_h_xyz_A": added_h_xyz.tolist(),
        "ligand_chain": ligand_chain, "ligand_resname": new_resname,
        "common_heavy_atom_map": {name: name for name in heavy_names},
        "shared_heavy_atom_count": 20, "shared_heavy_coordinates_changed": False,
        "carbon_oxygen_bond_A": carbon_oxygen_bond_A,
        "coordinate_assumption": "O replaces the stereochemically selected C16-H direction at a model 1.43 A C-O bond; shared E2 framework fixed",
        "replaced_h_atom": site["replaced_h_atom"], "tested_directions": site["tested_directions"],
        "ligand_hydrogens_exported": True, "ligand_hydrogen_count": 24,
        "formal_charge": Chem.GetFormalCharge(graft),
        "protonation_state": "neutral estriol with three hydroxyl hydrogens; no added formal positive charge",
        "hydrogens_optimized": False,
        "added_h_placement": "caller-supplied torsion with fixed model bond geometry" if added_oh_h_xyz_A is not None else "deterministic initial torsion; not optimized",
        "ligand_atoms": [{"name": _name(atom), "element": atom.GetSymbol(),
                           "xyz_A": graft.GetConformer().GetPositions()[atom.GetIdx()].tolist()}
                          for atom in graft.GetAtoms()],
        "ligand_bonds": [{"atom1": _name(bond.GetBeginAtom()), "atom2": _name(bond.GetEndAtom()),
                          "order": "aromatic" if bond.GetIsAromatic() else bond.GetBondTypeAsDouble()}
                         for bond in graft.GetBonds()],
        "added_h_bond_exclusions": ["H25", "O3", _name(molecule.GetAtomWithIdx(site["c16_index"]))],
        "protein_repacked": False, "source_coordinates_changed": False,
        "energetics_assessed": False,
    }
