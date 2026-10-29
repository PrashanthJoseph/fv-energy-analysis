"""Transfer authoritative ligand chemistry; PDB connectivity is not bond order."""
from __future__ import annotations

from pathlib import Path


class LigandIdentityError(ValueError):
    def __init__(self, diagnostics):
        self.diagnostics = diagnostics
        differences = diagnostics["differences"]
        super().__init__(
            f"Ligand identity mismatch after bond-order restoration: {', '.join(differences)}; "
            f"expected InChIKey={diagnostics['expected']['inchikey']}; "
            f"observed InChIKey={diagnostics['observed']['inchikey']}; "
            f"stereochemistry differences={diagnostics['stereocentre_differences']}"
        )

    def __reduce__(self):
        # ProcessPoolExecutor must reconstruct this exception from its
        # diagnostic mapping, rather than from ValueError's formatted message.
        return type(self), (self.diagnostics,)


def identity_diagnostics(observed, reference, original=None):
    """Explain identity differences without changing coordinates or stereochemistry."""
    from rdkit import Chem
    from rdkit.Chem import rdMolDescriptors
    def snapshot(molecule):
        explicit = sum(atom.GetAtomicNum() == 1 for atom in molecule.GetAtoms())
        with_h = Chem.AddHs(Chem.Mol(molecule))
        return {"formula": rdMolDescriptors.CalcMolFormula(molecule),
                "heavy_atom_count": molecule.GetNumHeavyAtoms(),
                "explicit_hydrogen_count": explicit,
                "total_hydrogen_count": sum(a.GetAtomicNum() == 1 for a in with_h.GetAtoms()),
                "formal_charge": Chem.GetFormalCharge(molecule),
                "smiles": Chem.MolToSmiles(Chem.RemoveHs(molecule), isomericSmiles=True),
                "inchikey": Chem.MolToInchiKey(molecule)}
    expected, actual = snapshot(reference), snapshot(observed)
    ref_heavy, obs_heavy = Chem.RemoveHs(reference), Chem.RemoveHs(observed)
    heavy_geometry = Chem.Mol(obs_heavy)
    Chem.AssignAtomChiralTagsFromStructure(heavy_geometry, replaceExistingTags=True)
    Chem.AssignStereochemistry(heavy_geometry, cleanIt=True, force=True)
    connectivity = Chem.MolToSmiles(ref_heavy, isomericSmiles=False) == Chem.MolToSmiles(obs_heavy, isomericSmiles=False)
    differences = [field for field in ("formula", "heavy_atom_count", "total_hydrogen_count", "formal_charge")
                   if expected[field] != actual[field]]
    if not connectivity:
        differences.append("connectivity_or_bond_orders")
    stereo = []
    match = obs_heavy.GetSubstructMatch(ref_heavy, useChirality=False)
    if len(match) == ref_heavy.GetNumAtoms():
        for index, observed_index in enumerate(match):
            old, new = ref_heavy.GetAtomWithIdx(index), obs_heavy.GetAtomWithIdx(observed_index)
            old_cip = old.GetProp('_CIPCode') if old.HasProp('_CIPCode') else None
            new_cip = new.GetProp('_CIPCode') if new.HasProp('_CIPCode') else None
            if old_cip != new_cip:
                info = new.GetPDBResidueInfo()
                stereo.append({"reference_atom_index": index + 1, "observed_heavy_atom_index": observed_index + 1,
                               "atom_name": info.GetName().strip() if info else None,
                               "expected_CIP": old_cip, "observed_CIP": new_cip})
    if stereo or (connectivity and expected["inchikey"] != actual["inchikey"]):
        differences.append("stereochemistry")
    return {"stage": "after_bond_order_restoration", "expected": expected, "observed": actual,
            "heavy_atom_geometry_inchikey": Chem.MolToInchiKey(heavy_geometry),
            "before_restoration": snapshot(original) if original is not None else None,
            "connectivity_and_bond_orders_match": connectivity,
            "stereocentre_differences": stereo, "differences": differences,
            "observed_bonds": [{"atom1": b.GetBeginAtom().GetPDBResidueInfo().GetName().strip() if b.GetBeginAtom().GetPDBResidueInfo() else b.GetBeginAtomIdx() + 1,
                                "atom2": b.GetEndAtom().GetPDBResidueInfo().GetName().strip() if b.GetEndAtom().GetPDBResidueInfo() else b.GetEndAtomIdx() + 1,
                                "order": str(b.GetBondType())} for b in observed.GetBonds()]}


def _params_molecule(params, template=None):
    """Build chemical topology from params, with elements/stereo from template."""
    from rdkit import Chem
    params = Path(params)
    template = Path(template) if template else params.with_name(params.stem + '_0001.pdb')
    if not template.is_file():
        raise FileNotFoundError(f'Ligand template is required to read atom elements: {template}')
    molecule = Chem.MolFromPDBFile(
        str(template), removeHs=False, sanitize=False, proximityBonding=False
    )
    if molecule is None or not molecule.GetNumAtoms():
        raise ValueError(f'Cannot read ligand template {template}')
    builder = Chem.RWMol(molecule)
    for first, second in [(bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()) for bond in builder.GetBonds()]:
        builder.RemoveBond(first, second)
    indices = {}
    for atom in builder.GetAtoms():
        name = atom.GetPDBResidueInfo().GetName().strip()
        if name in indices:
            raise ValueError(f'Duplicate template atom name: {name}')
        indices[name] = atom.GetIdx()
        atom.SetNoImplicit(True)
        atom.SetIsAromatic(False)
        atom.SetFormalCharge(0)
    atoms = []
    orders = {'1': Chem.BondType.SINGLE, '2': Chem.BondType.DOUBLE,
              '3': Chem.BondType.TRIPLE, '4': Chem.BondType.AROMATIC, 'ARO': Chem.BondType.AROMATIC}
    for line in params.read_text(encoding='utf-8').splitlines():
        fields = line.split()
        if not fields:
            continue
        if fields[0] == 'ATOM':
            atoms.append((fields[1], fields[2], float(fields[-1])))
        elif fields[0] in ('BOND', 'BOND_TYPE'):
            kind = fields[3].upper() if len(fields) > 3 else '1'
            builder.AddBond(indices[fields[1]], indices[fields[2]], orders[kind])
        elif fields[0] == 'CHARGE' and fields[2] == 'FORMAL':
            builder.GetAtomWithIdx(indices[fields[1]]).SetFormalCharge(int(fields[3]))
    if len(atoms) != len(indices) or {name for name, _, _ in atoms} != set(indices):
        raise ValueError('Ligand params and template must contain the same atom names')
    molecule = builder.GetMol()
    Chem.SanitizeMol(molecule)
    Chem.AssignAtomChiralTagsFromStructure(molecule, replaceExistingTags=True)
    Chem.AssignStereochemistry(molecule, cleanIt=True, force=True)
    return molecule, atoms


def reference_smiles_from_source(source=None, params=None, template=None):
    """Params/template control analysis; SDF/MOL2 is for generation input only.

    A supplied pair takes priority over source or sidecar metadata. When only
    params is passed, use its standard companion NAME_0001.pdb template.
    """
    from rdkit import Chem
    if params:
        molecule, _ = _params_molecule(params, template)
    elif template:
        raise ValueError('A ligand template must be accompanied by params')
    elif source:
        source = Path(source)
        if source.suffix.lower() == '.sdf':
            molecules = list(Chem.SDMolSupplier(str(source), removeHs=False))
            if len(molecules) != 1:
                raise ValueError('Expected one ligand in SDF')
            molecule = molecules[0]
        elif source.suffix.lower() == '.mol2':
            molecule = Chem.MolFromMol2File(str(source), removeHs=False)
        else:
            raise ValueError('Ligand source must be SDF or MOL2')
        if molecule is None:
            raise ValueError('Cannot read ligand chemical source')
    else:
        return None
    return Chem.MolToSmiles(Chem.RemoveHs(molecule), isomericSmiles=True)

# wwPDB CCD EST, CACTVS stereospecific descriptor; PubChem CID 5757.
# https://www.rcsb.org/ligand/EST
ESTRADIOL_SMILES = 'C[C@]12CC[C@H]3[C@@H](CCc4cc(O)ccc34)[C@@H]1CC[C@@H]2O'
ESTRADIOL_INCHIKEY = 'VOXZDWNPVJITMN-ZBRFXRBCSA-N'


def assign_reference_chemistry(molecule, smiles: str):
    """Preserve all coordinates/atom order while restoring bond orders and stereo.

    Supplied stereochemistry is verified, not silently changed to the reference.
    Explicit hydrogens are retained. Missing hydrogens are not invented here.
    """
    from rdkit import Chem
    from rdkit.Chem import AllChem
    reference = Chem.MolFromSmiles(smiles)
    if reference is None:
        raise ValueError('Invalid authoritative ligand SMILES')
    heavy_indices = [a.GetIdx() for a in molecule.GetAtoms() if a.GetAtomicNum() != 1]
    heavy = Chem.RemoveHs(molecule, sanitize=False)
    if heavy.GetNumAtoms() != reference.GetNumHeavyAtoms():
        raise ValueError('Ligand heavy-atom count differs from authoritative chemistry')
    restored = AllChem.AssignBondOrdersFromTemplate(reference, heavy)
    combined = Chem.RWMol(molecule)
    for atom in restored.GetAtoms():
        target = combined.GetAtomWithIdx(heavy_indices[atom.GetIdx()])
        target.SetIsAromatic(atom.GetIsAromatic())
        target.SetFormalCharge(atom.GetFormalCharge())
        target.SetNumRadicalElectrons(0)
        target.SetNoImplicit(True)
    for bond in restored.GetBonds():
        target = combined.GetBondBetweenAtoms(heavy_indices[bond.GetBeginAtomIdx()], heavy_indices[bond.GetEndAtomIdx()])
        target.SetBondType(bond.GetBondType())
        target.SetIsAromatic(bond.GetIsAromatic())
    result = combined.GetMol()
    Chem.SanitizeMol(result)
    Chem.AssignAtomChiralTagsFromStructure(result, replaceExistingTags=True)
    Chem.AssignStereochemistry(result, cleanIt=True, force=True)
    if Chem.MolToInchiKey(result) != Chem.MolToInchiKey(reference):
        raise LigandIdentityError(identity_diagnostics(result, reference, molecule))
    return result


def inspect_params_chemistry(params: str | Path, reference_smiles: str | None = None,
                             template: str | Path | None = None):
    """Read params chemistry using template elements, never atom-name guesses."""
    from rdkit import Chem
    molecule, atoms = _params_molecule(params, template)
    indices = {atom.GetPDBResidueInfo().GetName().strip(): atom.GetIdx() for atom in molecule.GetAtoms()}
    errors = []
    aromatic_names = [name for name, index in indices.items() if molecule.GetAtomWithIdx(index).GetIsAromatic()]
    atom_types = {name: kind for name, kind, _ in atoms}
    incorrect_aromatic_types = [name for name in aromatic_names if atom_types[name] != 'aroC' and molecule.GetAtomWithIdx(indices[name]).GetAtomicNum() == 6]
    if incorrect_aromatic_types:
        errors.append('Aromatic carbons have non-aromatic Rosetta atom types: ' + ','.join(incorrect_aromatic_types))
    net_charge = sum(q for _, _, q in atoms)
    by_index = {index: name for name, index in indices.items()}
    oxygens = [a for a in molecule.GetAtoms() if a.GetAtomicNum() == 8]
    phenolic = [a for a in oxygens if any(n.GetIsAromatic() for n in a.GetNeighbors())]
    aliphatic_oxygen = [a for a in oxygens if a not in phenolic]
    methyls = [a for a in molecule.GetAtoms() if a.GetAtomicNum() == 6 and sum(n.GetAtomicNum() == 1 for n in a.GetNeighbors()) == 3]
    roles = {}
    if len(aromatic_names) == 6 and len(phenolic) == len(aliphatic_oxygen) == len(methyls) == 1:
        roles = {'o3_atom': by_index[phenolic[0].GetIdx()], 'o17_atom': by_index[aliphatic_oxygen[0].GetIdx()],
                 'c18_methyl_atom': by_index[methyls[0].GetIdx()], 'a_ring_atoms': aromatic_names}
    if reference_smiles:
        reference = Chem.MolFromSmiles(reference_smiles)
        if reference is None:
            raise ValueError('Invalid reference ligand SMILES')
        if Chem.MolToSmiles(Chem.RemoveHs(molecule), isomericSmiles=False) != Chem.MolToSmiles(reference, isomericSmiles=False):
            errors.append('Parameter connectivity/bond orders do not represent authoritative ligand chemistry')
        expected_charge = Chem.GetFormalCharge(reference)
        if abs(net_charge - expected_charge) > 0.005:
            errors.append(f'Partial charges sum to {net_charge:+.6f}, expected {expected_charge:+d}')
    return {'valid': not errors, 'errors': errors, 'aromatic_atom_names': aromatic_names,
            'aromatic_bond_count': sum(b.GetIsAromatic() for b in molecule.GetBonds()),
            'partial_charge_sum': net_charge, 'formal_charge': Chem.GetFormalCharge(molecule),
            'atom_elements': {name: molecule.GetAtomWithIdx(index).GetSymbol() for name, index in indices.items()},
            'atom_types': atom_types, 'estrogen_atom_roles': roles}
