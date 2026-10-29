# Interaction fingerprint output definitions

## Purpose

These files keep interaction fingerprints readable without pretending that
Rosetta, Arpeggio, PLIP and Discovery Studio use identical definitions.

## Grouped summary (recommended)

`batch_summary_grouped_fingerprints.xlsx` contains one row per mutant, with
four visually grouped fingerprint sections. Each essential metric is its own
numeric column. The pinned reference row contains parent counts and scores.
The companion `batch_summary_grouped_fingerprints.csv` contains the same mutant
values with one flat, tool-prefixed header row. The reference values and the
section definitions are also saved in `batch_summary_grouped_fingerprints.json`.

| Section | Essential columns |
|---|---|
| Rosetta | Interface energy and delta, scored ligand H-bonds, VdW attraction/repulsion, electrostatics, solvation, `hbond_sc` energy only (not all H-bond terms). |
| Arpeggio | Native names: `hbond`, `weak_hbond`, `vdw`, `vdw_clash`, `hydrophobic`, `ionic`, `aromatic`, `CATIONPI`, `CARBONPI`, `DONORPI`; `plane-plane` counts the native record type. |
| PLIP | H-bonds, hydrophobic, salt bridges, pi-stacks, pi-cation, water bridges. |
| Discovery Studio | Conventional/carbon H-bonds, pi-stacked/T-shaped, pi-cation/pi-anion, pi-alkyl, alkyl, pi-sigma, attractive charge, salt bridges, unfavorable interactions. |
| HCDR3 pLDDT (when original CIFs are supplied) | Mean, minimum and one numeric column for each canonical heavy-chain residue H97–H113. |
| Boltz-2 (when the merged score CSV is supplied) | Confidence score, pTM, iPTM, ligand/protein iPTM, complex pLDDT/iPLDDT, complex PDE/iPDE, affinity value and binder probability including the source's suffixed 1/2 outputs. |
| GNINA (when the merged score CSV is supplied) | Rank1 and best values for minimizedAffinity, CNNscore, CNNaffinity, CNN_VS and CNNaffinity_variance. |

The two model-score sections import numeric values directly from the supplied
CSV. No score is recalculated or converted to Rosetta units. The Boltz-2 complex
pLDDT/iPLDDT values retain their source 0–1 scale, distinct from the H97–H113
residue confidence section's 0–100 scale. Affinity suffixes `1` and `2` are
retained as source labels without assuming that they identify structural
models. GNINA rank1 and best values remain separate, even when some happen to
coincide. Each `best_*` column is imported as supplied, not assumed to identify
one common pose. Absent parent scores remain blank, never zero. Original
chain/pair PTM matrix scores and source-row identifiers are retained in each
candidate's `model_scores.json` rather than widening the summary further.

`hcdr3_plddt_97` through `hcdr3_plddt_113` are the original ModelCIF local
pLDDT values (0–100), not relaxed-PDB B factors. `hcdr3_plddt_mean` is the
arithmetic mean of those 17 values; `hcdr3_plddt_min` is their minimum. The
complete canonical H sequence must uniquely match a CIF protein chain. The
source chain need not be named H. The interval uses the analysis PDB's residue
numbers, not a fresh IMGT or other numbering assignment. Missing or ambiguous
sources raise an error rather than producing zero scores. These values assess
prediction confidence and do not measure binding affinity.

Arpeggio `plane-plane` counts records once each, regardless of
the number of orientation labels on a record. Aromatic atom annotations are
counted separately. Arpeggio does not report pi-anion or pi-alkyl; no synthetic
zero columns are supplied for them. Weak H-bonds are not equated to carbon
H-bonds. Generic favorable VdW counts are unavailable from the DS monitor.
Arpeggio labels are retained as emitted, including case. `ionic` is not aliased
to salt bridges or attractive charge; `aromatic` is not renamed to pi-stacking.
Missing Rosetta/PLIP summary data are unavailable, not zero: JSON uses `null`,
CSV/Excel cells are blank (read as NaN by standard numeric-data loaders), and
Markdown displays `NaN`. A genuine recorded empty result retains zero. Unsupported
classes remain `NA`. Missing files still raise an error rather than being treated
as an empty detector run. Ligand identity comes from recorded configuration,
PLIP metadata or params, never a hard-coded E2L filter.
DS "Conventional H-bonds" counts only `conventional_hbond` records. The internal
derived `hbond_total` combines multiple categories, but is not the conventional
column and is not included as a DS total-H-bond column in the grouped XLSX.
Rosetta `hbond_sc` is only one weighted interface component. Backbone-side-chain
H-bond contributions use `hbond_bb_sc`; backbone-backbone terms are separately
reported as `hbond_sr_bb` and `hbond_lr_bb`. No combined H-bond energy is implied
by the explicitly labelled `hbond_sc only (REU)` column.

Earlier tidy, compact-text and four-text-column summaries are superseded by
this grouped layout. Existing saved files are preserved; new report generation
does not produce those redundant variants.

## `batch_interaction_fingerprint_numeric.csv`

This machine-readable table also has exactly one row per candidate. Each
supported tool/interaction combination has three numeric columns:

- `parent_<tool>_<interaction>`
- `mutant_<tool>_<interaction>`
- `delta_<tool>_<interaction>` (mutant minus parent)

For example, `parent_ds_pi_anion`, `mutant_ds_pi_anion` and
`delta_ds_pi_anion` describe Discovery Studio's pi-anion annotations.

| Column | Meaning |
|---|---|
| `candidate` | Mutant complex filename. |
| `result_folder` | Individual analysis folder. |
| `parent_<tool>_<interaction>` | Parent count from the named tool. |
| `mutant_<tool>_<interaction>` | Mutant count from the named tool. |
| `delta_<tool>_<interaction>` | Mutant minus parent count. |

Interaction classes unsupported by a tool are omitted from the flattened
table, avoiding columns that contain only `NA`. A recorded `0` means the tool
supports the class but found none. Counts must be compared within a tool, not
added across tools.

Rosetta is an energy function rather than a general interaction-fingerprint
detector. Its count row therefore reports scored ligand H-bond records only.
Its `fa_atr`, `fa_rep`, electrostatic, solvation and H-bond energies remain in
the Rosetta energy-decomposition outputs.

Arpeggio and Discovery Studio counts are atom/group-level annotations and can
overlap for the same residue pair. PLIP counts are PLIP interaction entries.

## Individual candidate outputs

- `interaction_fingerprint_counts.csv`: one flattened row containing the
  standardized counts for all four tools for one candidate.
- `interaction_fingerprint_all_detector_counts.csv`: every native interaction
  label actually emitted by each detector, with parent, mutant and delta.
- `interaction_fingerprint_summary.md`: readable standardized and native-count
  tables; this section is also included in `automatic_summary.md`.
- `interaction_fingerprint_counts.json`: machine-readable wide representation,
  including Rosetta interface energy components.

The original `batch_summary.csv` from the full scientific workflow remains
available for complete scoring, RMSD and mutation-network details.
# E2 architecture section

O3/O17 hydrogen bonds and C18 methyl hydrophobic recognition use Arpeggio.
Two-sided A-ring packing is a coordinate-based side-chain proximity screen, not a pi-stack assignment.
DS `pi_stack` / "Pi-stacks (unique residues)" counts each protein residue partner once, even when both fused rings of a Trp generate records. Native ring-pair records remain unchanged in detailed detector outputs. Other DS interaction classes retain their existing counting units.

The architecture pi-contact definition includes Discovery Studio stacked, T-shaped, pi-cation and pi-anion records involving the chemically mapped E2 A ring. Count each unique protein partner once across all four classes. Class-specific counts remain in the existing DS section, not duplicate columns in the architecture section.
The pi-contact sandwich requires two distinct partners on opposite A-ring faces, each with a partner-group centroid normal offset of at least 0.5 A. For charged partners use the DS interacting atom group, not an invented aromatic ring. Multiple groups/rings in one residue are retained for geometry but count as one partner. This broader project definition is not a classical pi-stack sandwich.
The recognition score counts four features (O3, O17, C18 and two-sided packing). The strict screen additionally requires the broader DS pi-contact sandwich; this is a project-defined filter, not a binding-energy or affinity prediction. Detailed outputs retain a separate stacked-only diagnostic.
Missing DS output leaves stack counts and the strict-screen result unavailable, not zero/false. The architecture status column identifies missing prerequisites. Parent values use the same method as mutants.
