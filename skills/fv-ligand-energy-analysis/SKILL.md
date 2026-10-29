---
name: fv-ligand-energy-analysis
description: Run and interpret the bundled Rosetta/PyRosetta parent-versus-mutant Fv-ligand analysis for one relaxed complex or a folder of relaxed complexes, including ANARCII/AbNumber CDR RMSDs, mutation chemistry and rotamers, Rosetta energy decomposition, Arpeggio/PLIP interactions, mutation-network changes, H-bonds, salt bridges, and estradiol-pocket screening. Use when asked to analyze, compare, rank, debug, resume, or explain outputs from this Energy_analysis project; do not use for unrelated docking or experimental affinity claims.
---

# Fv-ligand energy analysis

Use the modules in this repository as the computational source of truth. Do not recreate the analysis with ad hoc scripts when the packaged workflow can perform it.

## Locate the project and read only what is needed

Find the project root containing `pyproject.toml`, `src/fv_energy_analysis/`, and `scripts/bootstrap_wsl.sh`. When this skill remains in the repository, the root is two directories above this file.

- For installation, command options, or an unfamiliar input format, read `README.md`.
- For a column or interaction-label question, read `docs/batch_output_columns.md`.
- For a failed run, inspect the relevant manifest and log before changing code.

Keep this project separate from ABB4 or other antibody-design packages. Inputs and outputs may live elsewhere, but do not copy their source into this repository.

## Current E2 production protocol

The agreed protocol uses one whole-complex OptH pass before scientific analysis,
then RosettaLigand `ligand.wts` with
`restore_pre_talaris_2013_behavior` for both parent and mutants. Preserve heavy
atoms and chemical states during OptH. Use the exported optimized PDBs as input
to all four interaction tools; do not combine pre-OptH fingerprints with
post-OptH energies. Use the same validated ligand parameters throughout.

E2 O3/O17 H-bonds, C18 packing and A-ring pi-mediated contact labels come from
native Discovery Studio only. The steric sandwich is separately labelled
coordinate geometry using the current two-sided-enclosure definition in
`docs/e2_architecture.md`. Arpeggio and PLIP remain independent fingerprint
sections; they do not supply E2 architecture flags. Reuse compatible parent
artifacts rather than repeatedly optimizing or scoring the parent.

## Choose the execution mode

Use the individual workflow for one parent-mutant comparison:

```bash
cd /path/to/Energy_analysis
source .venv/bin/activate
export PYTHONPATH="$PWD/src"
python -m fv_energy_analysis.workflow \
  --parent /path/to/parent_relaxed.pdb \
  --mutant /path/to/mutant_relaxed.pdb \
  --output /path/to/results/mutant_vs_parent \
  --weights ligand \
  --params /path/to/E2L.params \
  --ligand-template-pdb /path/to/E2L_0001.pdb
```

Use the batch workflow for a folder. It creates one child folder per complex and ranked master reports:

```bash
cd /path/to/Energy_analysis
source .venv/bin/activate
export PYTHONPATH="$PWD/src"
python -m fv_energy_analysis.batch_workflow \
  --parent /path/to/parent_relaxed.pdb \
  --input-folder /path/to/relaxed_complexes \
  --output /path/to/results/batch_vs_parent \
  --jobs 2 \
  --weights ligand \
  --params /path/to/E2L.params \
  --ligand-template-pdb /path/to/E2L_0001.pdb
```

The batch module is parent-cached and parallel by default. It computes the
parent Rosetta/H-bond/Arpeggio/PLIP/numbering evidence once, reuses it for every
mutant, keeps PyRosetta and the ANARCII sequence cache warm inside persistent
workers, and processes `--jobs` mutants concurrently. Keep the default of two
workers unless available RAM supports more PyRosetta processes or the user has
chosen a different count (three is the current local run setting). Do not replace
this optimized module with a shell loop over the individual workflow.

Add chain maps or explicit ligand-chain flags only when the structures require them. Common options are:

```text
--mutant-chain-map H:L,L:H
--parent-ligand-chain d
--mutant-ligand-chain B
--recursive
--glob "*.pdb"
```

Use WSL/Linux for the supported environment. Convert Windows paths to mounted WSL paths when invoking from WSL. Do not hardcode a particular user's virtual-environment path in portable instructions.

## Environment and input gates

Before the first scientific run in a copied or changed environment:

```bash
bash scripts/bootstrap_wsl.sh
source .venv/bin/activate
python -m fv_energy_analysis.doctor
```

The doctor must confirm PyRosetta, AbNumber/ANARCII, Gemmi, RDKit, Arpeggio, PLIP, and Open Babel. PyRosetta remains subject to its licence.

Require these invariants:

- Parent and mutant Fv chains have identical polymer lengths and differ only by amino-acid substitutions.
- Parent and mutant use the same ligand identity, atom naming, protonation/tautomer state, and Rosetta parameterization.
- Compare relaxed with relaxed structures, using the same score function and settings.
- Reuse the same validated ligand `.params` and template PDB for every candidate in a batch.
- Treat missing or inconsistent ligand hydrogens, connectivity, or atom names as a blocking validation error; do not silently score mismatched ligands.

PDB and mmCIF are accepted. Remap nonstandard antibody chains through the workflow rather than manually rewriting scored structures.

## Run behavior

- Prefer the single module invocation; it owns all stages and records status in its manifest.
- Do not use `--force` unless the user asks to recompute completed analyses or inputs/settings changed materially.
- Re-running a batch safely reuses child analyses whose manifests are already `complete` and retries incomplete or failed children.
- Reuse the configuration-keyed `_parent_cache_*` folder. Do not delete it
  between compatible runs; a changed parent timestamp, score weights, ligand
  parameters, chain map, or numbering settings automatically selects a new key.
- Do not repeatedly poll a CPU-bound batch. Start it through a durable foreground execution/session, allow it to run, and use the manifest/log for a later status request. Never rely on a transient detached Windows/WSL shell that may terminate when its launcher closes.
- If the user wants to run it themselves, give one complete command with their actual paths; do not require a notebook.
- Do not claim success from a process launch. Completion requires `status: complete` in `workflow_manifest.json`, or `complete`/`complete_with_failures` in `batch_manifest.json` with the failure count stated.

## Failure handling

For an individual failure, read:

1. `workflow_manifest.json`
2. `workflow.log`
3. the last successfully written stage files

For a batch failure, read `batch_manifest.json` and `batch.log`, then inspect only the failed child folders. Diagnose the actual exception or exit code. Do not fix a crash by suppressing errors, splitting the scientifically required stages, or substituting a different scoring method without explicit approval.

After a code change, run the relevant tests and a focused regeneration from saved intermediate tables when possible. Do not rerun PyRosetta merely to change report wording or aggregation.

## What the workflow measures

The full workflow performs:

- ligand preparation and explicit-hydrogen/connectivity validation;
- PyRosetta pair energies and ligand-interface scoring with the selected protocol;
- Rosetta donor-acceptor H-bond extraction;
- Arpeggio geometric contact classification;
- PLIP protein-ligand cross-checking on the same validated coordinates;
- native Discovery Studio interactions on those same optimized coordinates;
- ANARCII-backed AbNumber numbering and region assignment;
- one global Fv fit followed by global-frame RMSDs for Fv, framework, and all six CDRs;
- amino-acid chemical-class comparison and rotamer-bin comparison;
- mutation-centered interaction-network formation, loss, retention, and energy decomposition;
- intra-chain, inter-chain, and protein-ligand H-bond and salt-bridge changes;
- literature-inspired estradiol-pocket architecture metrics.

Rotamers are comparable only when the relevant CDR RMSD after the global Fv fit is at most the configured cutoff, 3.0 A by default. Otherwise report that rotamer analysis was skipped. Do not substitute local-fit RMSD for this decision.

## Energy interpretation

All deltas are `mutant - parent`:

- negative is favorable in the Rosetta model;
- positive is unfavorable;
- REU are score-function units, not experimentally measured kcal/mol;
- multi-mutant energy changes are non-additive and cannot be assigned causally to one substitution.

The interface energy is the selected score function's score of the bound pose
minus the score after rigidly translating the ligand 1000 A away without
unbound-state repacking. Current E2 production uses legacy RosettaLigand
`ligand.wts`; `ref2015` is available for separately requested protocols and must
not be mixed with RosettaLigand results. This is a fixed-conformation interaction
score, not a complete binding free energy. It omits or approximates conformational
entropy, alternate states, explicit solvent reorganization, protonation
equilibria, and unbound relaxation.

Characterize interface results on separate axes:

1. **Pose-level effect size:** magnitude and sign of `interface_delta_REU`.
2. **Relative performance:** rank and percentile within candidates scored by the identical protocol.
3. **Term consistency:** whether attraction, solvation, H-bond, electrostatic, and repulsive terms support or oppose the total.
4. **Structural confidence:** CDR RMSDs, clashes, ligand validation, and whether the binding mode is conserved.
5. **Robustness:** replicate/ensemble distribution, not the single best score.

Do not call a Rosetta delta statistically significant from one pose. Close-ranked candidates are the same practical tier until replicate relaxation or ensemble scoring separates them. For an ensemble, report the median delta, IQR or SD, 95% interval, fraction of replicates with delta below zero, rank stability, and recurrence of the proposed interaction network.

The package's pair-energy strength labels are convenience descriptors only:

```text
<= -2 REU   strong favorable
<= -1 REU   moderate favorable
< 0 REU     weak favorable
0 to <1     neutral or weak unfavorable
>= 1 REU    unfavorable
```

Never present these thresholds as experimental affinity cutoffs.

## Interaction interpretation

Keep these quantities distinct:

- **Atom-level annotations:** overlapping atom-pair labels; 14 annotations do not mean 14 independent interactions.
- **Residue-pair presence/count:** one protein-residue/partner pair, regardless of how many atoms or labels contribute.
- **Interaction class:** geometry/type annotation such as H-bond, ionic, hydrophobic, pi-stack, or proximal.
- **Rosetta pair energy:** energetic score for the complete residue pair.

For every claimed gain or loss, name the exact pair, scope, interaction type, parent energy, mutant energy, and delta. A formed geometric contact is not automatically favorable, and a broken geometric contact is not automatically an energetic loss.

Arpeggio, Rosetta, and PLIP use different definitions. Report disagreement rather than forcing counts to match. `proximal` means close in space, not favorable. `vdw_clash` may coexist with favorable annotations because different atom pairs contribute to the same residue pair. `energy_only` means the Rosetta threshold changed without a newly assigned geometric class.

H-bond pairs are deduplicated when Rosetta and Arpeggio recognize the same residue pair. Salt bridges are based on Arpeggio `ionic` annotations and are protonation-dependent. Estradiol is neutral, so protein-ligand salt bridges should be treated skeptically unless the modeled chemistry explicitly supports charge.

Do not add energy totals across interaction-class rows: one multi-annotated residue pair is repeated in every applicable class. Use the deduplicated changed-pair total for an overall network energy.

## Estradiol-specific screen

Use four recognition features: DS conventional O3 H-bond, DS conventional O17
H-bond, DS singleton-C18 alkyl/pi-alkyl contact, and two-sided A-ring enclosure.
The latter is `e2_geometric_sandwich`: all qualifying side-chain heavy atoms are
searched, with nearest ring-carbon Bondi surface gaps 0–1.0 A and absolute plane
height at least the atom's own radius. A pair from distinct residues on opposite
faces qualifies when its joining segment crosses the projected union of the
six carbon van der Waals disks (radius 1.70 A). Count unique residue pairs and
report overlaps separately. The former carbon-normal-ray rule is retired.

DS A-ring stacked/T-shaped pi-pi, pi-cation and pi-anion partners form a separate
unique-residue pi-mediated descriptor. They are not an additional pass
requirement. The paper's enclosure need not contain two aromatic residues or
two pi stacks. Numerical enclosure conditions are project screening choices,
compatible with the sixteen examined reference copies and m7, not published
thresholds or proof of E2 specificity. C19 obstruction is a separate diagnostic
and must not be inferred from the enclosure flag. Use method-version checks
when refreshing saved architecture so historical ray results are not reused.

## Reporting

For an individual comparison, use `automatic_summary.md` and its CSV/JSON evidence. For a batch, start with:

- `batch_summary.md` and `batch_summary.csv`;
- `batch_mutation_summary.csv`;
- `batch_interaction_network_summary.csv`;
- `batch_hbond_saltbridge_changes.csv`;
- `e2_architecture_screen.csv`;
- `COLUMN_DEFINITIONS.md`.

Lead with the result, then state limitations. A candidate verdict should include:

- interface parent, mutant, delta, and batch rank/percentile;
- bound-total and major interface-term changes;
- global, framework, and mutated-CDR RMSDs;
- each substitution's chemical-class change and rotamer result or skip reason;
- exact formed/broken ligand, intra-chain, and inter-chain interactions with energy deltas;
- lost or gained salt bridges;
- ligand H-bond counts from each detector without conflating atom events and unique residue partners;
- conserved E2 architecture and strict-screen result;
- a verdict separating effect size from confidence.

Use “supports,” “suggests,” or “is consistent with,” not “proves.” Recommend experimental validation. Recommend an ensemble especially when CDR RMSD exceeds 3 A, the ranking margin is small, a conclusion depends on one pose, or detector disagreement affects the mechanism.
