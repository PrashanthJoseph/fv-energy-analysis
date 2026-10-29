# Third-party components

The project does not bundle the scientific applications installed by its setup
scripts. Installing this project does not replace those applications' licenses.

## Separately installed scientific tools

PyRosetta/Rosetta and BIOVIA Discovery Studio must be obtained and used under
their own terms. They are not included in this repository or its wheel.
Rosetta's `molfile_to_params.py` is likewise supplied by the user's Rosetta
installation when new ligand parameters are required.

Arpeggio, PLIP, Open Babel, RDKit, AbNumber, ANARCII, Gemmi, NumPy and openpyxl
are separate dependencies. Their source distributions retain their respective
licenses and notices; the installer downloads them rather than incorporating
their source code into this repository. See the dependency guide in the README
for the install commands and external-tool requirements.
