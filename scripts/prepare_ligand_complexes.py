"""Map bound ligands exactly to shared params/template chemistry; retain XYZ."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import re

from rdkit import Chem
from fv_energy_analysis.ligand_chemistry import assign_reference_chemistry, LigandIdentityError


def ligand_fragment(pdb: Path, chain: str) -> Chem.Mol:
    pdb_lines = pdb.read_text(encoding="utf-8").splitlines()
    ligand_lines = [
        line for line in pdb_lines
        if line.startswith("HETATM") and line[21:22].strip() == chain
    ]
    ligand_serials = {int(line[6:11]) for line in ligand_lines}
    # Retain explicit ligand connectivity when it is present.  Omitting CONECT
    # makes RDKit infer bonds solely from distance; tightly packed relaxed poses
    # can then acquire spurious bonds and impossible carbon valences.
    conect_lines: list[str] = []
    for line in pdb_lines:
        if not line.startswith("CONECT"):
            continue
        try:
            serials = [int(value) for value in line[6:].split()]
        except ValueError:
            continue
        if not serials or serials[0] not in ligand_serials:
            continue
        retained = [serials[0], *(value for value in serials[1:] if value in ligand_serials)]
        if len(retained) > 1:
            conect_lines.append("CONECT" + "".join(f"{value:5d}" for value in retained))
    molecule = Chem.MolFromPDBBlock(
        "\n".join(ligand_lines + conect_lines + ["END"]) + "\n",
        removeHs=False,
        sanitize=True,
        proximityBonding=not bool(conect_lines),
    )
    if molecule is None:
        raise ValueError(f"RDKit could not read ligand chain {chain!r} in {pdb}")
    fragments = Chem.GetMolFrags(molecule, asMols=True, sanitizeFrags=True)
    if len(fragments) != 1:
        raise ValueError(f"Expected one ligand fragment for chain {chain!r} in {pdb}, found {len(fragments)}")
    return fragments[0]


def ligand_template(pdb: Path, params: Path | None = None) -> Chem.Mol:
    """Read bonds/charges from params and atom names/stereochemistry from PDB.

    Existing callers pass the standard Rosetta NAME_0001.pdb template; its
    companion is NAME.params. --params supports any other filename/location.
    """
    params = params or pdb.with_name(re.sub(r"_\d{4}$", "", pdb.stem) + ".params")
    molecule = Chem.MolFromPDBFile(str(pdb), removeHs=False, sanitize=False, proximityBonding=False)
    if molecule is None or not molecule.GetNumAtoms():
        raise ValueError(f"Cannot read ligand template {pdb}")
    builder = Chem.RWMol(molecule)
    for first, second in [(bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()) for bond in builder.GetBonds()]:
        builder.RemoveBond(first, second)
    indices = {}
    for atom in builder.GetAtoms():
        name = atom.GetPDBResidueInfo().GetName().strip()
        if name in indices:
            raise ValueError(f"Duplicate template atom name: {name}")
        indices[name] = atom.GetIdx()
        atom.SetNoImplicit(True)
        atom.SetIsAromatic(False)
        atom.SetFormalCharge(0)
    orders = {"1": Chem.BondType.SINGLE, "2": Chem.BondType.DOUBLE,
              "3": Chem.BondType.TRIPLE, "4": Chem.BondType.AROMATIC, "ARO": Chem.BondType.AROMATIC}
    parameter_names = set()
    for line in params.read_text(encoding="utf-8").splitlines():
        fields = line.split()
        if not fields:
            continue
        if fields[0] == "ATOM":
            parameter_names.add(fields[1])
        elif fields[0] in {"BOND", "BOND_TYPE"}:
            kind = fields[3].upper() if len(fields) > 3 else "1"
            builder.AddBond(indices[fields[1]], indices[fields[2]], orders[kind])
        elif fields[0] == "CHARGE" and fields[2] == "FORMAL":
            builder.GetAtomWithIdx(indices[fields[1]]).SetFormalCharge(int(fields[3]))
    if parameter_names != set(indices):
        raise ValueError("Ligand params and template must contain the same atom names")
    molecule = builder.GetMol()
    Chem.SanitizeMol(molecule)
    Chem.AssignAtomChiralTagsFromStructure(molecule, replaceExistingTags=True)
    Chem.AssignStereochemistry(molecule, cleanIt=True, force=True)
    return molecule


def rename_ligand_atoms(
    source: Path,
    destination: Path,
    chain: str,
    output_chain: str,
    template_pdb: Path,
    params: Path | None = None,
) -> None:
    lines = source.read_text(encoding="utf-8").splitlines(keepends=True)
    ligand_indices = [
        index
        for index, line in enumerate(lines)
        if line.startswith("HETATM") and line[21:22].strip() == chain
    ]
    template = ligand_template(template_pdb, params)
    reference_smiles = Chem.MolToSmiles(Chem.RemoveHs(template), isomericSmiles=True)
    try:
        source_ligand = assign_reference_chemistry(ligand_fragment(source, chain), reference_smiles)
    except LigandIdentityError as error:
        diagnostic = {"source": str(source), "chain": chain, **error.diagnostics}
        destination.parent.joinpath('ligand_identity_diagnostic.json').write_text(json.dumps(diagnostic, indent=2) + '\n')
        raise
    if (len(ligand_indices) != template.GetNumAtoms()
            or source_ligand.GetNumAtoms() != template.GetNumAtoms()
            or source_ligand.GetNumBonds() != template.GetNumBonds()):
        raise ValueError("Bound ligand must have the same atoms, explicit hydrogens and bonds as the params/template")
    match = source_ligand.GetSubstructMatch(template, useChirality=True)
    if len(match) != template.GetNumAtoms():
        raise ValueError("Bound ligand does not exactly map to the params/template chemistry and stereochemistry")
    residue_name = template.GetAtomWithIdx(0).GetPDBResidueInfo().GetResidueName()
    mapped_names = [""] * source_ligand.GetNumAtoms()
    for template_index, source_index in enumerate(match):
        mapped_names[source_index] = template.GetAtomWithIdx(template_index).GetPDBResidueInfo().GetName()
    for index, name in zip(ligand_indices, mapped_names):
        raw = lines[index].rstrip("\r\n").ljust(80)
        lines[index] = (
            f"{raw[:12]}{name}{raw[16:17]}{residue_name}{raw[20:21]}"
            f"{output_chain}{raw[22:]}\n"
        )
    last = ligand_indices[-1]
    if last + 1 >= len(lines) or not lines[last + 1].startswith("TER"):
        lines.insert(last + 1, "TER\n")
    for index, line in enumerate(lines):
        if line.startswith("#BEGIN_POSE_ENERGIES_TABLE"):
            lines = lines[:index]
            break
    Chem.Kekulize(template, clearAromaticFlags=True)
    serial_by_name = {line[12:16].strip(): int(line[6:11]) for line in lines
                      if line.startswith('HETATM') and line[21:22].strip() == output_chain}
    ligand_serials = set(serial_by_name.values())
    retained = []
    for line in lines:
        if line[:6].strip() == 'END':
            continue
        if line.startswith('CONECT') and int(line[6:11]) in ligand_serials:
            continue
        retained.append(line)
    adjacency = {serial: [] for serial in ligand_serials}
    for bond in template.GetBonds():
        a = serial_by_name[bond.GetBeginAtom().GetPDBResidueInfo().GetName().strip()]
        b = serial_by_name[bond.GetEndAtom().GetPDBResidueInfo().GetName().strip()]
        order = int(bond.GetBondTypeAsDouble())
        adjacency[a].extend([b] * order)
        adjacency[b].extend([a] * order)
    lines = retained + ['CONECT' + f'{serial:5d}' + ''.join(f'{other:5d}' for other in partners) + '\n'
                        for serial, partners in sorted(adjacency.items()) if partners]
    if not lines or not lines[-1].startswith("END"):
        lines.append("END\n")
    destination.write_text("".join(lines), encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--reference", type=Path, help=argparse.SUPPRESS)
    parser.add_argument("--complex", type=Path, action="append", required=True)
    parser.add_argument("--generated-ligand-pdb", type=Path, required=True)
    parser.add_argument("--params", type=Path, help="Shared ligand params; defaults to the template's companion NAME.params")
    parser.add_argument("--ligand-chain", default="d")
    parser.add_argument("--output-ligand-chain", default="X")
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--ligand-smiles", help=argparse.SUPPRESS)  # Older callers; never overrides params/template.
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)

    for source in args.complex:
        destination = args.output / f"{source.stem}_rosetta.pdb"
        rename_ligand_atoms(
            source,
            destination,
            args.ligand_chain,
            args.output_ligand_chain,
            args.generated_ligand_pdb,
            args.params,
        )
        print(destination)


if __name__ == "__main__":
    main()
