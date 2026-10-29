# Flex ddG for Fv–estradiol

This PC is configured to run the published protein–small-molecule Flex ddG
adaptation through the existing PyRosetta installation in Ubuntu-24.04:
`/home/inito_joey/venvs/fv-pyrosetta/bin/python` (2026.29 quarterly build).
No separate compiled `rosetta_scripts` executable is required. The original
ligand XML is in `vendor/flex_ddg/aldeghi_2018/SI_Files/Rosetta/ddG-backrub.xml`.

Open `run_flex_ddg.ipynb`, edit its mutation dictionary, then run the cells.
The notebook works with a Windows Python kernel or the configured Linux
environment. It runs calculations through the installed Linux interpreter.

## Mutation dictionary

```python
MUTATIONS = {"H": "E99F"}

# These entries form ONE combined mutant:
MUTATIONS = {"H": ["E99F", "R100D", "I102Y"]}

# Entries across chains also form ONE combined mutant:
MUTATIONS = {"H": "E99F", "L": "Y91F"}
```

The notation is **original amino acid + PDB residue number + new amino acid**.
Positions refer to the chain labels and numbers in the input PDB, not automatic
IMGT/Kabat numbering or zero-based sequence indexes. A position with PDB
insertion code A can be written `E99AF`. Values may be a single string or a
list of strings. The runner rejects missing positions, incorrect original
amino acids, duplicate positions, and substitutions that do not change a residue.
To compare separate mutants, run each dictionary in its own output folder.

## Run directly from Windows PowerShell

```powershell
cd D:\INITO\myProjects\Energy_analysis
./run_flex_ddg.ps1 -Mutations @{ H = 'E99F' }

# Combined mutant, custom output folder:
./run_flex_ddg.ps1 -Mutations @{ H = @('E99F','R100D','I102Y') } -Output flex_ddg/runs/triple

# Fast execution check (results are not affinity predictions):
./run_flex_ddg.ps1 -Mutations @{ H = 'E99F' } -Output flex_ddg/runs/smoke -Smoke

# Validate chosen mutations without calculating:
./run_flex_ddg.ps1 -Mutations @{ H = 'E99F' } -Output flex_ddg/runs/check -DryRun
```

Use a fresh output folder when changing a mutation or settings. Completed
replicas are reused only when input-file hashes, settings, protocol/runner
hashes, and the installed PyRosetta distribution version match the stored job.
Windows paths with spaces are supported. The launcher selects Ubuntu-24.04
explicitly, regardless of the default WSL distribution.

## Run from Linux/Python

```python
from fv_energy_analysis.flex_ddg import run_flex_ddg

result = run_flex_ddg(
    pdb="flex_ddg/inputs/parent.pdb",
    params="flex_ddg/inputs/E2L.params",
    mutations={"H": "E99F"},
    output="flex_ddg/runs/H_E99F",
    ligand_chain="X",
    replicas=35,
    steps=35000,
    workers=1,
)
```

Use the installed PyRosetta interpreter and set `PYTHONPATH=src` when running
from this project. Alternatively invoke `scripts/run_flex_ddg.py` using that
interpreter and supply `--pdb`, `--params`, `--mutations-file`, and `--output`.
The example dictionary is in `examples/flex_ddg_mutations.json`.

## Default inputs and estradiol chemistry

The default parent is the project's previously prepared
`e2pico_scfv_m7_estradiol_model_0_relaxed_rosetta.pdb`, copied into
`flex_ddg/inputs/parent.pdb` with freshly parameterized estradiol.
Its antibody chains are H and L; estradiol is chain X.
In this input, H99=E, H100=R and H102=I.

The previous ligand parameters in the existing analysis results encoded every
bond as single. RDKit interpreted that connectivity as C18H30O2 with zero
aromatic atoms. The separate new Flex ddG input is verified against
[PubChem estradiol CID 5757](https://pubchem.ncbi.nlm.nih.gov/compound/5757):
C18H24O2, six aromatic carbons, and the expected 17beta-estradiol InChIKey
`VOXZDWNPVJITMN-ZBRFXRBCSA-N`. All bound heavy-atom coordinates are preserved
within PDB rounding; ligand hydrogens are regenerated. The earlier result
files have not been recalculated with this corrected chemistry.

To prepare another estradiol-bound PDB, run in the Linux environment:

```bash
/home/inito_joey/venvs/fv-pyrosetta/bin/python scripts/prepare_flex_ddg_estradiol.py \
  --pdb /path/to/complex.pdb --ligand-chain X --output flex_ddg/inputs_new
```

This requires an intact connected estradiol structure and checks stereochemistry
from its 3D coordinates. For a different input, select the resulting PDB and
params together; atom names must match. The ligand should follow the antibody
chains and occupy its own chain with one residue. Resolve alternate locations
and supply one structural model.

## Protocol and outputs

Defaults: `beta_nov16` scoring, 35 independent replicas, 35,000 backrub moves,
5,000 minimization iterations per cycle, score convergence threshold 1.0,
extra chi1/chi2 rotamers, one replica at a time. Increase `workers` only when
resources allow; each replica uses a separate process and deterministic seed.
`ref2015` and `talaris2014` are also available with `--weights`/`-Weights`.
Smoke mode uses two replicas, ten backrub moves, and shortened minimization.

The runner uses the authors' XML for constrained minimization, local backrub
ensembles, paired WT/mutant repacking and minimization, and the original
`InterfaceDdGMover` four-state score reporter. It makes these explicit changes:

- Sets the selected score function, ligand chain, mutations and sampling settings.
- Freezes ligand packing/torsions and jump minimization, retaining its internal
  bound conformation. It does not explicitly sample alternate ligand poses.
  Upstream protein motion can carry a ligand rigidly without changing its conformation.
- Replaces CLI-only structure database reporting with WT and mutant PDB snapshots.
- Uses a separate Linux temporary working directory per replica for SQLite performance.

The final directory contains:

- `summary.json`: mean/median raw ddG, standard deviation, standard error,
  average component changes and historical GAM-reweighted estimate.
- `replicas.csv`: per-replica WT binding score, mutant binding score and ddG.
- `job.json`, `mutations.resfile`, `ligand_validation.json`, `status.json`:
  exact settings, validated mutation selection and execution status.
- `replica_###/`: actual XML, execution log, original four-state `ddG.db3`,
  `wt_bound.pdb`, `mutant_bound.pdb`, final backrub structure and numerical result.

`ddG = (mutant_bound − mutant_separated) − (WT_bound − WT_separated)`.
Negative favors binding. Raw scores are Rosetta energy units (REU).
The GAM calculation applies the original seven-term nonlinear transformation
to the averaged component changes, matching the tutorial's aggregation order.
It is reported separately and is not an estradiol-specific calibration.
The sampling SEM measures variation among simulations, not experimental accuracy.
Neither score predicts association/dissociation rates separately.

The paper's benchmark used Rosetta 2017.52. This setup uses the existing 2026
build, so it is an adaptation of the published protocol, not an exact reproduction
of its benchmark. Smoke outputs verify execution only and must not be used to
rank mutations. A full chosen-mutation ensemble is a separate calculation.

## Sources

- [Aldeghi, Gapsys & de Groot, 2018: ligand adaptation](https://doi.org/10.1021/acscentsci.8b00717).
- [Authors' original ligand input archive](https://acs.figshare.com/articles/journal_contribution/7462244), downloaded from `https://ndownloader.figshare.com/files/13820498`.
- [Kortemme Lab Flex ddG tutorial](https://github.com/Kortemme-Lab/flex_ddG_tutorial).
- [Original Flex ddG paper](https://doi.org/10.1021/acs.jpcb.7b11367).

The existing licensed PyRosetta installation is reused; its license still applies.
