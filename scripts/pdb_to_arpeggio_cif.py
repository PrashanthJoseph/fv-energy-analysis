"""Convert a PDB to Arpeggio's mmCIF input with component annotations."""

from __future__ import annotations

import argparse
from pathlib import Path

import gemmi


AMINO_ACIDS = {
    "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE",
    "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL",
    "MSE",
}


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("pdb", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()

    structure = gemmi.read_structure(str(args.pdb))
    document = structure.make_mmcif_document()
    block = document.sole_block()
    residue_names = sorted(
        {residue.name for model in structure for chain in model for residue in chain}
    )
    loop = block.init_loop("_chem_comp.", ["id", "type", "name"])
    for name in residue_names:
        component_type = "'L-PEPTIDE LINKING'" if name in AMINO_ACIDS else "NON-POLYMER"
        loop.add_row([name, component_type, name])
    args.output.parent.mkdir(parents=True, exist_ok=True)
    document.write_file(str(args.output))


if __name__ == "__main__":
    main()
