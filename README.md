# Fv energy analysis

The reusable package compares parent and mutant Fv–ligand complexes, caches
the parent, processes mutants in parallel, and writes per-candidate evidence
plus a sortable grouped CSV/XLSX. Scientific inputs must use the same validated
ligand chemistry and equal-length antibody chains.

## Complete production folder workflow

`fv-energy-full` is the installed controller for the complete folder workflow:
one whole-complex OptH pass → RosettaLigand scoring and interaction analysis →
native Discovery Studio → DS-based E2 recognition features → optional
fixed-pocket E3 discrimination → CSV and formatted XLSX. It keeps completed
work resumable, reuses compatible parent artifacts and runs blocking child
stages without background polling. The historical `tmp` helpers and the Codex
desktop app are not needed.

See [accepted inputs and configurable controls](docs/configuration.md) for
score weights, analysis switches, structural/network thresholds and selectable
CSV/XLSX sections or columns. List defaults with `fv-energy-full --describe-options`.

The complete four-tool setup uses a Windows reporting/controller environment
and a WSL scientific environment, because native Discovery Studio runs on
Windows. Install the scientific environment once under WSL:

```bash
cd /path/to/Energy_analysis
bash scripts/bootstrap_wsl.sh
source .venv/bin/activate
fv-energy-doctor --profile analysis
```

Then install the lightweight Windows controller from the same checkout:

```powershell
py -3.12 -m venv .venv-windows
.\.venv-windows\Scripts\python -m pip install ".[reporting]"
.\.venv-windows\Scripts\fv-energy-doctor --profile reporting
```

Copy `examples/full_pipeline_windows.json`, supply your actual input, output,
ligand-parameter, WSL Python and Discovery Studio paths, then run:

```powershell
.\.venv-windows\Scripts\fv-energy-full examples/full_pipeline_windows.json --preflight
.\.venv-windows\Scripts\fv-energy-full examples/full_pipeline_windows.json
```

The example uses the agreed production settings: three workers, one OptH pass,
`ligand.wts` with the legacy pre-Talaris behavior, DS-only interaction labels
for E2 recognition, and the current two-sided-enclosure geometry. E3 is an
optional steric assessment of the added **C16α OH**, not a binding-free-energy
calculation. It does not establish experimental E2/E3 selectivity.

Dependencies are declared in `pyproject.toml`; the bootstrap installs ordinary
Python dependencies, PLIP and the authorized PyRosetta build. Native Discovery
Studio remains a separate licensed installation. See
[dependency handling and environment records](docs/dependencies.md).

To regenerate only the workbook from saved results:

```bash
fv-workbook /path/to/completed/batch
fv-workbook /path/to/completed/batch --verify-only
```

This validates the saved data and uses Python/openpyxl, with no scientific
rescoring, Node.js or Codex dependency.

For a two-complex integration check, use
[the Jupyter notebook](examples/full_pipeline_two_complexes.ipynb). Copy
`examples/two_complex_settings.example.json` to `two_complex_settings.json`,
fill in the input and environment paths, and launch the notebook with a Windows
Python kernel (`pip install notebook ipykernel`). Set `SETTINGS_PATH` in its
first cell, or the `FV_NOTEBOOK_SETTINGS` environment variable. It clones the
public repository at a recorded commit, installs its Windows controller, runs
all analysis sections, and verifies the resulting CSV/XLSX. This example expects
already hydrogen-optimized inputs and does not repeat OptH. A new output folder
is used by default; set `resume` to true to reuse a compatible completed test.

## Flex ddG for selectable Fv–estradiol mutations

The optional `fv-flex-ddg` module accepts mutations such as `{"H": "E99F"}`
or a combined mutant such as `{"H": ["E99F", "R100D"]}` using input PDB
numbering. It requires an installed PyRosetta environment. See
[the setup and usage guide](docs/flex_ddg.md) for commands, settings and
interpretation.

## Individual scientific workflow

Run the individual scientific comparison inside the activated WSL environment. This
module invocation performs ligand preparation, PyRosetta scoring, H-bond and
Arpeggio and PLIP extraction, ANARCII numbering, RMSDs, mutation chemistry/rotamers, and
report merging:

```bash
source .venv/bin/activate
cd /path/to/Energy_analysis
export PYTHONPATH="$PWD/src"

python -m fv_energy_analysis.workflow \
  --parent /path/to/parent.pdb \
  --mutant /path/to/mutant.pdb \
  --ligand-source /path/to/ligand.sdf \
  --weights ligand \
  --output results/mutant_vs_parent
```

The module writes `workflow.log` while it runs and sets `status` to `complete`
or `failed` in `workflow_manifest.json`.  Final tables are placed in `final/`,
with concise machine-readable and Markdown summaries at
`automatic_summary.json` and `automatic_summary.md`.
Use a new or empty output directory: the individual workflow refuses to overwrite
existing files. Discovery Studio remains a separate Windows post-processing
stage; DS-dependent architecture criteria are unavailable until its results are supplied.
This lower-level command does not run OptH first; use `fv-energy-full` for the
production preparation and all-stage folder workflow.
`ligand_atom_validation.json` verifies ligand atom names, explicit hydrogen
completeness, parameterized connectivity, and bonded distances before scoring.
PLIP is run as an independent protein-ligand interaction cross-check on those
same validated coordinates with `--nohydro`; its XML/TXT reports are saved in
`plip_parent/` and `plip_mutant/`, and `plip_ligand_interactions.json` contains
compact counts and contributing residues.

The final directory also contains `mutation_network_changes.csv` (every
mutation-partner edge, formed/broken interaction types, contact-count changes,
pair-energy components and strength) and `mutation_network_summary.csv`.
`mutation_network_type_energy_summary.csv` names the contributing residue pairs and aggregates pair energies and score
terms for every formed/broken interaction type and records nominal charge
changes. Because contact types overlap, a multi-type pair contributes its pair
energy to each applicable type sum.
`mutation_interpretation_summary.csv` provides a conservative per-mutation
assessment with direct-ligand evidence, net local-network energy, strongest
supporting and opposing pair changes, chemistry, rotamer evidence, and the
non-additive multi-mutant context.
CDR rotamers are evaluated only when the relevant RMSD after the global fit is
at or below 3.0 A: anchor-only RMSD for HCDR3, overall region RMSD for other CDRs.
Change the threshold with `--rotamer-rmsd-cutoff`. Framework mutations are not gated.

To reuse an already validated ligand definition, supply both
`--params /path/to/E2L.params` and
`--ligand-template-pdb /path/to/E2L_0001.pdb`.
These two files are sufficient: `--ligand-source` is only needed when generating
new parameters. On reuse, chemical bond orders come from the parameters and
stereochemistry from the matching template; bound coordinates stay in the input
complex. An adjacent `ligand_definition.json` is optional provenance, not a
required third input.

If the input PDBs use different ligand chain IDs, add, for example,
`--parent-ligand-chain d --mutant-ligand-chain B`.  Both are normalized to the
Rosetta ligand chain internally before comparison.

PDB and mmCIF inputs are accepted. If an input has swapped or nonstandard
antibody chain labels, remap them during conversion, for example
`--mutant-chain-map H:L,L:H`.

## Batch analysis of a relaxed-complex folder

The batch module runs the complete workflow independently for every matching
complex, places each analysis in its own subfolder, and writes ranked master
tables. Completed child analyses are reused when the command is run again;
add `--force` to rescore them. A failed candidate is recorded without stopping
the remaining folder unless `--stop-on-error` is supplied.

Failures are appended to `batch_issues.csv` with the candidate, failing stage,
actual exception message and diagnostic-log path. Ligand identity/chemistry
rejections are labelled `ligand_validation`; they are not bypassed or scored as
successful cases. Retrying an incomplete candidate appends to its `workflow.log`
instead of erasing the previous attempt. Final summaries retain failed cases.

Ligand preparation saves a validated `NAME.atom_roles.json` alongside the
parameters, bound to SHA-256 hashes of both params and template. Reporting reads
this mapping without RDKit. Native Discovery Studio supplies production E2
O3/O17 H-bond, C18 packing and pi-mediated interaction labels; the steric
enclosure is separately labelled coordinate geometry. Arpeggio remains an
independent fingerprint section.
Missing or stale mappings are explicit errors, not guessed atom identities.
`fv-discovery-studio-batch --reports-only` verifies and reuses saved detector
outputs; `--skip-fingerprint-reports` avoids duplicate fingerprint generation
when a following reporting stage owns it.

Model-score reporting accepts either the merged Boltz-2/GNINA table or a native
Boltz affinity table keyed by `job_name`. The latter supplies only its six
affinity/probability fields; missing confidence, GNINA and parent values stay
blank. Exact model matching is required, including validation of repeated-job
PDB filename wrappers. For metadata supplied after a batch starts, a
`report_sources.json` in its output folder can specify `scores_csv` and its
`scores_csv_sha256`; the reporting stage reads it without restarting scientific
analysis. The same file can supply `cif_folder` and a verified `parent_cif` for
original ModelCIF pLDDT. `--include-empty-metadata-sections` preserves the full grouped layout.

Identity rejections now write `ligand_identity_diagnostic.json` with expected
and observed formulas, H counts, charges, SMILES, InChIKeys and named
stereocentre differences. Stereochemistry is also assessed from heavy-atom
geometry independently. These diagnostics do not change the input coordinates
or force a rejected stereoisomer to match the reference.

```bash
python -m fv_energy_analysis.batch_workflow \
  --parent /path/to/parent_relaxed.pdb \
  --input-folder /path/to/relaxed_complexes \
  --output results/my_batch \
  --jobs 2 \
  --weights ligand \
  --mutant-chain-map H:L,L:H \
  --params results/reference/prepared/E2L.params \
  --ligand-template-pdb results/reference/prepared/E2L_0001.pdb
```

The batch implementation is optimized for large folders. It builds a
configuration-keyed parent cache once, reuses the parent's prepared structure,
Rosetta pair/interface scores, H-bonds, Arpeggio, PLIP, E2 architecture and
numbering work, and scores only each mutant thereafter. Long-lived worker
processes initialize PyRosetta once, load each mutant pose once for both energy
and H-bond output, and retain the ANARCII sequence cache across their assigned
mutants. `--jobs 2` is the conservative default; increase it only when the
machine has enough RAM for additional PyRosetta processes. Completed child
folders and a matching parent cache remain resumable.

During cache creation, the shared workflow recognizes the identical parent
inputs: preparation, combined Rosetta scoring/H-bond extraction, the configured
Arpeggio stage, PLIP and pocket analysis each run once. Parent/mutant filenames
in this self-comparison are aliases of those results, not duplicate calculations.

If only the parent reference changes, the completed mutant calculations can be
reused instead of rescoring every mutant:

```bash
python -m fv_energy_analysis.batch_workflow \
  --parent /path/to/new_parent_relaxed.pdb \
  --input-folder /path/to/relaxed_complexes \
  --output results/batch_vs_new_parent \
  --reuse-mutants-from results/completed_old_batch \
  --params results/reference/prepared/E2L.params \
  --ligand-template-pdb results/reference/prepared/E2L_0001.pdb \
  --jobs 2
```

This evaluates the new parent once and regenerates only parent-relative RMSDs,
substitutions, pair-energy deltas, formed/broken interaction networks, rankings,
and reports. The saved batch must contain the same candidates and use the same
ligand parameterization.

The batch manifest records the cache path, worker count, worker process IDs and
whether each result was reused. A parent cache key includes the parent file
timestamp, ligand parameter/template paths, score weights, chain mapping and
antibody numbering settings, so scientifically different configurations do not
silently share cached results.

By default, the module scans for `*_relaxed.pdb`, `*_relaxed.cif`, and
`*_relaxed.mmcif`. Use repeatable `--glob` arguments and `--recursive` when
needed. The parent output folder contains:

- `batch_summary.csv` and `batch_summary.md`: one ranked row per complex;
- `batch_mutation_summary.csv`: one row per substitution across all complexes;
- `batch_interaction_network_summary.csv`: formed/broken interaction classes, the contributing residue pairs, and their energy sums;
- `e2_architecture_screen.csv`: compact E2-specific ranking and the residues satisfying each conserved pocket feature;
- E2 recognition columns: DS O3/O17 hydrogen bonds, DS C18 hydrophobic recognition, geometric A-ring enclosure, a 0–4 architecture score and strict pass/fail; DS pi-mediated partners are a separate descriptor;
- `COLUMN_DEFINITIONS.md`: definitions, units, sign conventions, and interpretation guidance for every batch-table column;
- `batch_manifest.json` and `batch.log`: resumable status and execution log;
- one named subfolder per complex with the full individual analysis.

### Optional Discovery Studio cross-check

On Windows systems with BIOVIA Discovery Studio Visualizer installed, the
saved batch can be augmented with BIOVIA's complete ligand nonbond vocabulary,
including pi-anion, pi-cation, pi-pi, pi-alkyl, H-bond, salt-bridge and
unfavorable-contact calls. This is a report-only post-processing stage and does
not rerun Rosetta, Arpeggio, PLIP, ANARCII or structural alignment:

```powershell
$env:PYTHONPATH = "D:\path\to\Energy_analysis\src"
python -m fv_energy_analysis.discovery_studio `
  --batch-results "D:\path\to\results\my_batch"
```

The module auto-detects the standard Discovery Studio 2024 installation at
`C:\Program Files\BIOVIA\Discovery Studio 2024`. Use
`--discovery-studio-root` for another version or location. It invokes the
bundled BIOVIA Perl API in one persistent process, calculates the parent once,
and uses short temporary filenames to avoid Discovery Studio's legacy Windows
path-length limit. Hydrogens are not re-added.

See `docs/discovery_studio_columns.md` for the generated files and column
definitions. BIOVIA Discovery Studio is optional third-party software and is
not installed by this project.

The post-processor also creates deliberately compact batch tables:

- `batch_summary_grouped_fingerprints.xlsx`: the recommended summary, with one
  row per mutant and separate numeric metric columns under four grouped tool
  headings: Rosetta, Arpeggio, PLIP and Discovery Studio. A pinned parent
  reference row makes the baseline visible during scrolling;
- `batch_summary_grouped_fingerprints.csv`: the corresponding numeric mutant
  values, with flat tool-prefixed headers for analysis software;
- `batch_interaction_fingerprint_numeric.csv`: the corresponding one-row-per-
  candidate numeric table with separate parent, mutant and delta columns for
  filtering or calculations.

Every candidate folder receives `interaction_fingerprint_counts.csv` plus
`interaction_fingerprint_all_detector_counts.csv`. The latter retains every
native Rosetta, Arpeggio, PLIP and Discovery Studio interaction label. Both
tables are rendered into `interaction_fingerprint_summary.md` and appended to
the candidate's `automatic_summary.md`.

Regenerate these reports from saved analyses with
`python -m fv_energy_analysis.fingerprints --batch-results PATH`.
Install `python -m pip install '.[reporting]'` once and add `--export-xlsx` to
create the formatted workbook. Excel export uses the optional `openpyxl`
Python dependency and does not require Node.js or this desktop app.

This package compares a parent and mutant Fv–ligand complex. It requires the Fv chains to have identical residue counts and allows only amino-acid substitutions. It produces:

- one global rigid-body fit and RMSDs for the full Fv, combined framework, and every heavy/light CDR;
- both RMSD after the shared global fit and an explicitly separate local-fit diagnostic RMSD;
- AbNumber region definitions backed by ANARCII;
- mutation chemical classes and shared/lost/gained physicochemical properties;
- canonical side-chain chi angles and explicitly coarse `p/m/t` torsion bins (not Dunbrack assignments);
- Rosetta residue-pair score terms, interaction counts from PDBe Arpeggio, mutant-minus-parent deltas, and mutation-partner subsets.

## Installation

The source can be installed as a normal Python package with
`python -m pip install .`. For optional Excel export use `.[reporting]`.
See [CONTRIBUTING.md](CONTRIBUTING.md) for module responsibilities, unit tests,
package builds and the remaining license decisions before public release.

The supported setup is Linux or WSL. From a freshly copied project folder,
run the bootstrap script; it creates a project-local `.venv` and installs the
analysis dependencies:

```bash
cd /path/to/Energy_analysis
bash scripts/bootstrap_wsl.sh
source .venv/bin/activate
python -m fv_energy_analysis.doctor
```

The bootstrap installs AbNumber/ANARCII, scientific dependencies, PDBe
Arpeggio, and PLIP 3.0.1 with binary Open Babel Python bindings. PyRosetta is installed from RosettaCommons' official quarterly
package index and remains subject to the Rosetta licence. If its download is
not available to a recipient, install an authorised PyRosetta build in the
same environment and rerun the doctor.

Rosetta's ligand-parameter generator is an external licensed dependency. When
generating new ligand parameters, supply `--molfile-to-params PATH` or set
`FV_MOLFILE_TO_PARAMS`. This is not needed when reusing `--params` together
with `--ligand-template-pdb`.

PLIP is installed automatically by the bootstrap and does not need to be
installed again. To add it once to an environment that was created before
PLIP support was added, run:

```bash
bash scripts/install_plip_wsl.sh /path/to/existing/venv
source /path/to/existing/venv/bin/activate
python -m fv_energy_analysis.doctor
```

The installer is idempotent: rerunning it verifies the pinned installation
without creating a second environment. It reads the required binding versions
from this package's installed dependency metadata.

AbNumber 0.4.4 or newer exposes the ANARCII-backed API through
`abnumber.future.Chain` or `Chain(..., use_anarcii=True)`. The code tries both
interfaces.

The source folder is self-contained and has no dependency on ABB4. Keep
coordinate files, generated results, and virtual environments outside a source
archive (they are ignored by `.gitignore`). Inputs and outputs may be anywhere
and are always selected with command-line paths.

On Windows, invoke the workflow through WSL after setup:

```powershell
wsl bash -lc "cd /path/to/Energy_analysis && source .venv/bin/activate && python -m fv_energy_analysis.workflow --help"
```

## Inputs

- Parent and mutant coordinate files must be PDB files with matching Fv chain IDs.
- Each Fv chain must have the same number of resolved polymer residues in both structures.
- Ligands should have identical chain IDs, residue identities, protonation/tautomer states and Rosetta params in both structures.
- Rosetta tables are produced by `residue_energy_breakdown` with the same score function for both structures.
- Interaction files are JSON outputs from `pdbe-arpeggio`. Arpeggio itself takes mmCIF input.

## Run Rosetta decomposition

### Directly with the local WSL PyRosetta environment

The local quarterly PyRosetta build can generate both required tables without standalone Rosetta executables:

```bash
source .venv/bin/activate
cd /path/to/Energy_analysis
export PYTHONPATH="$PWD/src"

python -m fv_energy_analysis.cli score-pyrosetta \
  --pdb parent_relaxed.pdb \
  --params ligand.params \
  --ligand-chain X \
  --pair-output parent_energy.sc \
  --interface-output parent_interface.sc

python -m fv_energy_analysis.cli score-pyrosetta \
  --pdb mutant_relaxed.pdb \
  --params ligand.params \
  --ligand-chain X \
  --pair-output mutant_energy.sc \
  --interface-output mutant_interface.sc
```

The pair terms are weighted exactly once using the requested score function. Backbone hydrogen bonds are explicitly decomposed into pair energies. Interface values are calculated as bound minus a state in which the ligand chain is translated 1000 Å along its Rosetta rigid-body jump.

### With standalone Rosetta

```bash
fv-energy run-rosetta \
  --executable /path/to/residue_energy_breakdown.linuxgccrelease \
  --pdb parent_relaxed.pdb \
  --params ligand.params \
  --weights ref2015 \
  --output parent_energy.sc

fv-energy run-rosetta \
  --executable /path/to/residue_energy_breakdown.linuxgccrelease \
  --pdb mutant_relaxed.pdb \
  --params ligand.params \
  --weights ref2015 \
  --output mutant_energy.sc
```

Rosetta pose numbering is mapped to PDB residues using the coordinate-file residue order. Do not reorder or remove residues between generating the PDB and score table.

For global and component-wise ligand interface scores, edit the ligand chain in `examples/interface_score_X.xml` and run the same protocol on both structures:

```bash
rosetta_scripts.linuxgccrelease \
  -s parent_relaxed.pdb \
  -extra_res_fa ligand.params \
  -parser:protocol examples/interface_score_X.xml \
  -out:file:scorefile parent_interface.sc \
  -out:nooutput

rosetta_scripts.linuxgccrelease \
  -s mutant_relaxed.pdb \
  -extra_res_fa ligand.params \
  -parser:protocol examples/interface_score_X.xml \
  -out:file:scorefile mutant_interface.sc \
  -out:nooutput
```

## Run Arpeggio

Both automated workflows default to `--arpeggio-scope neighbourhood`. They keep
the complete prepared structures, but query mutated residues, whole protein
residues within 6 A of a mutation in **either** state, and the complete ligand.
The same union of residue IDs is used in both states. Thus lost/new contacts,
inter-/intra-chain mutation partners and contacts between queried neighbours
are retained without scanning unrelated protein regions. The ligand is always
queried, so ligand fingerprints cover the entire pocket, not only mutations.

`--arpeggio-neighbourhood-radius 6.0` changes the protein-neighbour radius, not
Arpeggio's native interaction cutoffs. Arpeggio also expands its query to nearby
atoms internally. No structure is cropped and no chemical parameters change.
Default exclusion of sequence-adjacent nonbonded contacts remains unchanged.

Batch execution calculates the parent's full contacts once and filters that
immutable cache to each candidate's query. It does **not** rerun parent Arpeggio
per mutant. `arpeggio_selection` in each workflow manifest records queried
residues and scope. Full saved contacts remain usable; narrower saved contacts
cannot be reused for a wider/different query unless they cover it. Existing
corrected Rosetta scores need not be recalculated for a contact-scope change.
Use `--arpeggio-scope full` when whole-protein geometric coverage is required.

Convert the same relaxed structures to mmCIF without changing author chain/residue IDs, then run:

```bash
pdbe-arpeggio -sf common_residue_selection.txt -o arpeggio_parent parent_relaxed.cif
pdbe-arpeggio -sf common_residue_selection.txt -o arpeggio_mutant mutant_relaxed.cif
```

Pass the resulting JSON files to the comparison command.

## Compare

```bash
fv-energy compare \
  --parent parent_relaxed.pdb \
  --mutant mutant_relaxed.pdb \
  --fv-chains H,L \
  --ligand-chains X \
  --scheme imgt \
  --cdr-definition imgt \
  --fit-region fv \
  --rmsd-atoms backbone \
  --parent-rosetta parent_energy.sc \
  --mutant-rosetta mutant_energy.sc \
  --parent-arpeggio parent.json \
  --mutant-arpeggio mutant.json \
  --parent-hbonds parent_hbonds.csv \
  --mutant-hbonds mutant_hbonds.csv \
  --parent-interface parent_interface.sc \
  --mutant-interface mutant_interface.sc \
  --output comparison
```

Use `--fit-region framework` if you want the single global transform determined only from framework atoms. All reported `rmsd_after_global_fit_A` values still use that one transform. `local_fit_rmsd_A` removes rigid-body displacement of each individual region and should not be substituted for the global-frame result.

## Add original CIF confidence to an existing batch

Install the reporting dependencies once with `pip install -e ".[reporting]"`.
This command only reads saved analyses and the original prediction CIFs; it does
not rerun Rosetta, relaxation, Arpeggio, PLIP or Discovery Studio:

```bash
fv-interaction-fingerprints --batch-results /path/to/results \
  --cif-folder /path/to/original_cifs --parent-cif /path/to/parent.cif \
  --export-xlsx
```

The grouped summary gains an **HCDR3 pLDDT** section containing H97–H113,
their mean and their minimum (0–100). The original CIF must explicitly declare
a local pLDDT metric. The complete canonical H-chain sequence is matched to
the CIF chain before transferring scores by sequence position; a CIF chain
letter is not assumed to identify the heavy chain. These residue numbers refer
to the normalized analysis PDB, not a newly assigned antibody-numbering scheme.
Relaxed-PDB B factors are not used. pLDDT describes the original prediction's
local confidence, not binding energy or experimental affinity.

To append confidence without rebuilding any fingerprint reports, use
`fv-hcdr3-confidence --batch-results /path/to/results --cif-folder /path/to/original_cifs`
(optionally `--parent-cif /path/to/parent.cif`). This updates CSV/JSON/Markdown
reports only. Each candidate gets `hcdr3_plddt.csv` and `hcdr3_plddt.json` with
source-chain and residue mappings. A configured CIF folder is remembered in the
batch manifest so later fingerprint-report regeneration retains the section.

### Add Boltz-2 and GNINA scores

Append scores already calculated in a merged CSV without rerunning any analysis:

```bash
fv-model-scores --batch-results /path/to/results \
  --scores-csv /path/to/boltz2_gnina_merged.csv
```

To regenerate the full report and Excel file with both score sections, add
`--scores-csv /path/to/boltz2_gnina_merged.csv` to
`fv-interaction-fingerprints --batch-results /path/to/results --export-xlsx`.
The source CSV is remembered in the batch manifest for subsequent regeneration.
Matching uses the exact model/CIF basename, ignoring only the relaxed suffix and
file extension. Model indices are retained. Ambiguous or missing matches fail
explicitly. The new columns retain the CSV's original score values and scales.
Parent scores remain unavailable when the source does not include the parent.

## Outputs

- `structure_rmsd.csv`: global Fv, framework, individual CDR, HCDR3-anchor and HCDR3-tip RMSDs.
- `mutant_aligned_to_parent.pdb`: mutant after the global fit.
- `substitutions.csv`: sequence validation, AbNumber positions/regions, chemical classes and rotamers.
- `rosetta_pairs_parent.csv`, `rosetta_pairs_mutant.csv`: all pairwise Rosetta terms.
- `interactions_parent.csv`, `interactions_mutant.csv`: Arpeggio contact types and counts.
- `pair_comparison.csv`: parent/mutant values and mutant-minus-parent deltas.
- `mutation_partner_impacts.csv`: pair-comparison rows involving a substituted residue.
- `summary.json`: run settings and interpretation warnings.

If supplied with `--parent-interface` and `--mutant-interface`, `interface_metrics_comparison.csv` contains all numeric `InterfaceScoreCalculatorMover` columns—including `interface_delta_X` and `if_X_<scoreterm>`—plus mutant-minus-parent deltas. Both scorefiles must be supplied together.

## Optional individual Rosetta H-bonds

If PyRosetta is available, extract donor–acceptor-level H-bond energies:

```bash
fv-energy extract-hbonds --pdb parent_relaxed.pdb --params ligand.params --output parent_hbonds.csv
fv-energy extract-hbonds --pdb mutant_relaxed.pdb --params ligand.params --output mutant_hbonds.csv
```

Add `--parent-hbonds parent_hbonds.csv --mutant-hbonds mutant_hbonds.csv` to `compare`. The merged pair table then includes the H-bond count, summed and mean individual energy, and strongest H-bond energy. Arpeggio remains the detailed geometric-label source and Rosetta the energy source; the full workflow also runs PLIP independently and reports agreement or disagreement in the automatic summary.

### Rotamer interpretation

Each available chi angle is assigned to `p` (0° to 120°, gauche+), `m` (-120° to 0°, gauche-), or `t` (the remaining trans range). The report distinguishes:

- `same_chi1_bin`: meaningful for substitutions when both side chains have chi1;
- `same_common_chi_bins`: compares all chi positions shared by the two residue types;
- `same_full_rotamer_bin`: only for the same amino-acid type with all expected chi angles present; different amino-acid types do not share a residue-specific rotamer identity.

This is a transparent coarse-bin comparison, not a lookup of the backbone-dependent rotamer probability. Rosetta's `fa_dun` term remains the energetic measure of side-chain conformational preference.

## Ligand preparation and validation (2026-10-19 correction)

Do not generate ligand chemistry from PDB coordinates alone. Supply an explicit
SDF/MOL2 file. The module uses installed Open Babel to convert SDF to aromatic
MOL2, then runs the licensed Rosetta parameter generator:

```bash
python -m fv_energy_analysis.ligand_parameters \
  --input /path/to/ligand.sdf --output /path/to/new_ligand_parameters
```

No coordinates are minimized, aligned, or moved during the Open Babel conversion;
explicit input hydrogens are retained. The generated Rosetta template may be
translated by Rosetta's generator, but its coordinates are not substituted for
the bound ligand: the bound complex keeps its own coordinates.

The output includes `ligand.mol2`, `E2L.params`, `E2L_0001.pdb`, an atom mapping,
a bond mapping, and `ligand_definition.json`. MOL2 partial charges are retained
at eight decimals rather than accepting the generator's two-decimal rounding.
For SDF input, Open Babel's default MOL2 output uses Gasteiger charges; these
are not AM1-BCC/RESP charges and this model choice must be kept consistent.
No arbitrary charge redistribution or aromatic atom-type patching is applied.

Pass `--params` and `--ligand-template-pdb` from this output to the workflow.
The adjacent definition file supplies the actual reference chemistry; or pass
`--ligand-source` explicitly. E2 atom roles are determined from the parameter
chemical graph, not hardcoded names: a different ligand atom order must not
silently change which ring/oxygen is being analyzed.

Legacy scores produced with the former all-single-bond/aliphatic E2 parameters
are not compatible with corrected parameters. They must be preserved as legacy
results and not mixed into a corrected parent/mutant comparison. Corrected
scoring/reuse records a protocol version, score weights, and parameter hash.

### Additional scoring and reporting corrections

- Pair export includes long-range Rosetta energy containers, including
  `rama_prepro` and `dslf_fa13`, and separate onebody records. Self plus all
  unique pair records reconstruct the score for the validated `ref2015` case.
- `total` can contain sequence-coupled/conformational and disulfide terms;
  `nonbonded_total` excludes `rama_prepro`, `dslf_fa13`, `cart_bonded`, and
  `pro_close`. These terms must not be called independent noncovalent contacts.
- Individual H-bond output distinguishes raw `hbond_energy` from
  `hbond_weighted_energy_REU`; only allowed/scored H-bonds are counted.
- Switching from an Arpeggio H-bond label to a Rosetta H-bond label is not, by
  itself, formation/loss of an H-bond residue-pair family.
- Network changed-pair summaries include retained pairs with energetic/count
  changes, not only pairs with newly appearing/disappearing type labels.
- The strict A-ring metric additionally checks non-degenerate aromatic planes,
  centroid distance <=5.5 A, plane angle <=30 degrees, and lateral displacement
  <=2 A on both planes. These are explicit project screening criteria, not a
  universal definition or an isolated pi-stacking energy.

### HCDR3 overall, anchors and tip

The existing overall HCDR3 RMSD is retained. Two additional metrics use the
same global Fv alignment and backbone atom set:

- `H_CDR3_TIP`: the positions of the unique `YYGSS` motif in the **parent** HCDR3;
- `H_CDR3_ANCHORS`: every other position within that parent HCDR3.

For m7, the tip is H104-H108; anchors are H97-H103 plus H109-H113. Selection is
transferred by matched residue positions, so mutations in the motif do not change
the selection. No separate tip/anchor fit is used for the main metrics. The
anchor-only HCDR3 RMSD is the <=3 A gate for coarse torsion comparisons throughout
HCDR3, including tip mutations. The overall HCDR3 RMSD remains reported separately.
Missing or ambiguous parent motifs yield an explicit unavailable result, not zero,
and HCDR3 rotamer comparisons are skipped without a whole-CDR fallback.
Use `--hcdr3-tip-sequence` to specify a different parent tip motif.

### Developer-only single-comparison numerical audit

Audit tooling lives under `tests/`, not in the runtime package or batch workflow.
From a source checkout:

```bash
python tests/manual/single_complex_sanity.py \
  --result-dir /path/to/completed_comparison \
  --params /path/to/new_ligand_parameters/E2L.params \
  --ligand-template-pdb /path/to/new_ligand_parameters/E2L_0001.pdb \
  --output /path/to/new_validation_folder
```

The checker compares direct two-body evaluation, score reconstruction,
independent Cartesian interface separation, weighted H-bond totals, Rosetta's
own canonical chi angles, an independent BioPython fit, tip/anchor squared-
displacement partitioning, and raw PLIP XML counts. It writes `validation.json`
and `validation.md`. Passing validates those implementations on that input;
it does not prove affinity accuracy, protonation state, ensemble robustness, or
that no possible implementation bug remains. Non-pair-decomposable score terms
in arbitrary custom score functions are outside this tested `ref2015` scope.

### Scoring-only rescore of a completed batch

To produce a **separate one-pass hydrogen-optimized RosettaLigand rescore**:

```bash
python -m fv_energy_analysis.hydrogen_rescore /path/to/completed_ligand_rescore \
  /path/to/new_output --jobs 3 --expected-count 1397
```

This performs one deterministic whole-complex proton-chi optimization per
complex, preserves all heavy atoms and chemical states, and retains the original
if the optimized bound score is worse. The matched parent is optimized once.
It writes `rosettaligand_hydrogen_optimized_rescore.csv` (one row per mutant),
`hydrogen_optimized_structures/` (mutants plus `parent_m7.pdb`), fresh Rosetta
H-bond records, per-case validation and a resumable manifest. Final energies are
calculated after reloading the exported PDB, not from unrounded coordinates.
Saved input hashes, scoring options, parameter hashes and output hashes must
match before completed results are reused. It does not regenerate Arpeggio,
PLIP, Discovery Studio or their network/fingerprint summaries.

For a separate parent/mutant hydrogen-placement diagnostic, run:

```bash
python -m fv_energy_analysis.hydrogen_benchmark \
  --rescore-root /path/to/completed_rescore --candidate candidate.pdb \
  --params /path/to/E2L.params --output /path/to/new_benchmark --repeats 3
```

This benchmarks standard whole-complex proton-chi rotamer sampling with the saved
score function, fixed heavy atoms and fixed chemical states. It is not continuous
hydrogen minimization or full relaxation. It requires canonical Ser H:35 and the
current estradiol O2/H2 mapping for its local diagnostic. It verifies saved baseline
scores, records startup and per-trial timings separately, and selects by the lowest
whole-complex bound score, retaining the original if no trial improves it. Results
are written to a new folder; production scores, fingerprints and workbooks are not
updated. Changing hydrogen positions requires regenerating hydrogen-dependent
fingerprints before producing a consistent new production report.

```bash
python -m fv_energy_analysis.rescore_saved /path/to/completed_batch \
  /path/to/separate_rescore_output --weights ligand --jobs 3
```

This reuses saved prepared structures and ligand parameters, with no relaxation,
docking, numbering or geometric-detector reruns. `ligand` selects the installed
hard-repulsive `ligand.wts` with `restore_pre_talaris_2013_behavior`; no tutorial
overrides are added. Legacy and modern scores use separate PyRosetta processes.
The supplied ligand partial charges are unchanged. This is fixed-coordinate
rescoring, not a reproduction of a benchmark that also minimizes structures.

The parent is scored once. Persistent workers process completed candidates only;
source failures remain excluded. Every result records the input, parameter and
weights hashes. A rerun in the same output reuses matching completed results and
refuses incompatible caches. The source batch is never modified.

Outputs include `rosettaligand_rescore.csv`, `rescore_manifest.json`, and per-case
pair-energy, interface and scored H-bond tables. The exporter validates pair-score,
interface-term and H-bond energy reconstruction. Use `--limit 3 --jobs 1` in a
separate validation output for a small initial numerical check.

To prepare the same grouped report with the new scores, first rebuild its
score-dependent mutation-network summaries, then merge the saved geometry and
model metrics (neither step invokes a scoring engine or fingerprint detector):

```bash
python -m fv_energy_analysis.rescore_network \
  --source /path/to/completed_batch --rescore /path/to/separate_rescore_output
python -m fv_energy_analysis.rescore_report \
  /path/to/completed_batch /path/to/separate_rescore_output
```

The new grouped CSV/JSON preserves the original sections and static metrics,
updates parent/mutant energies, scored H-bonds, ranks, network deltas and quality
warnings, and exposes the additional active legacy terms `fa_pair` and
`hbond_bb_sc`. Rejected source ligands remain listed as exclusions. The original
batch is unchanged. The grouped JSON is the input to the XLSX renderer.

### Energy strength

`pair_energy_strength` is only a convenience label: ≤ -2 REU strong favorable, ≤ -1 moderate favorable, < 0 weak favorable, 0–1 neutral/weak unfavorable, and ≥ 1 unfavorable. Retain and interpret the numeric score terms; these thresholds are not experimental affinity cutoffs.

### DS-only E2 recognition refresh

Reuse a completed single-OptH analysis without repeating expensive calculations:

```bash
python -m fv_energy_analysis.architecture_batch \
  --source /path/to/completed_single_OptH_batch \
  --output /path/to/new_architecture_output \
  --params /path/to/E2L.params --template /path/to/E2L_0001.pdb \
  --opth-results /path/to/single_pass_OptH_results \
  --c18-calibration /path/to/reference_calibration.json --jobs 3 --xlsx
```

This checks saved OptH hashes, named coordinates and DS artifact associations,
evaluates parent m7 once, and preserves all non-architecture scores and metrics.
O3/O17 and C18 use DS only. The steric sandwich uses the calibrated coordinate
definition of **two-sided A-ring enclosure**, replacing the old normal-ray rule.
It requires qualifying side-chain atoms from distinct residues on opposite
faces, surface gaps 0–1.0 A, and a joining segment crossing the ring's projected
carbon van der Waals footprint. Each atom must lie at least its own radius from
the ring plane. All qualifying atoms are searched, residue pairs are counted
once, and overlaps are retained separately. The numerical conditions are
project-defined, compatible with the sixteen examined paper copies and m7;
they do not establish C19 exclusion or E2 specificity. DS A-ring pi-mediated
contacts are a separate unique-residue descriptor. Version checks prevent reuse
of historical ray-based architecture results. No additional H replicates,
scoring or detector runs occur in this artifact-only refresh.
See [E2 recognition definitions and limitations](docs/e2_architecture.md).

### E3 discrimination: fixed-pocket sterics

Add a separate E3 section to a new copy of a completed combined report:

```bash
python -m fv_energy_analysis.e3_discrimination \
  --source /path/to/combined_batch \
  --params /path/to/E2L.params --template /path/to/E2L_0001.pdb \
  --output /path/to/new_e3_steric_output --jobs 3
```

This tests estriol's extra **16alpha-OH**, with its stereochemistry checked
against the authoritative molecular identity. The protein and twenty shared
steroid heavy atoms stay fixed. E3 has all twenty-four hydrogens explicit;
the new OH hydrogen is sampled over twelve torsions without repacking or
minimization. Primary oxygen-overlap measurements and possible polar contacts
remain separate. Outputs include per-case protonated probes, contact details,
reports and an updated combined CSV/JSON. Existing ranks, energies, detector
results and E2 recognition scores are preserved. No energetic or affinity
discrimination is inferred. Inputs must be the named, explicit-H E2L complexes
from the existing preparation pipeline. `--pdb` and `--input-folder` also allow
standalone assessments. See [E3 definitions and limitations](docs/e3_discrimination.md).

To correct a previously generated E3 batch using its saved protonated probes
and twelve saved OH-H orientations, without rebuilding ligands or rescoring:

```bash
python -m fv_energy_analysis.e3_refresh \
  --source /path/to/existing_e3_batch \
  --output /path/to/new_geometry_checked_e3_batch --jobs 3
```

Polar overlaps are exempted only when donor/acceptor chemistry and explicit
hydrogen geometry support the contact. All simultaneous exemptions must use
one common OH-H orientation. The selected saved H position is exported; all
other coordinates stay fixed. Missing geometry is distinct from a negative
steric warning. Non-E3 metrics and rankings remain unchanged.
