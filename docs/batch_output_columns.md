# Batch output column definitions

This guide describes the batch-level CSV files. It is copied into each batch output folder as `COLUMN_DEFINITIONS.md`.

## Conventions used throughout

- **REU** means Rosetta Energy Units. REU are scoring-function units, not experimentally measured kcal/mol.
- Every energy **delta** is `mutant - parent`. A negative energy delta is favorable in the Rosetta model; a positive delta is unfavorable.
- Distances and RMSDs are in angstroms (**A**). A lower RMSD means greater structural similarity, not automatically stronger binding.
- Interaction counts describe the analyzed structural pose. More contacts are not necessarily better: geometry and energy must also be considered.
- Results are pose-level computational evidence, not statistical significance or experimental proof.

## `batch_summary.csv`

### Identity and run status

| Column | Meaning |
|---|---|
| `candidate` | Filename of the mutant complex. |
| `candidate_path` | Full path to the mutant input structure. |
| `result_folder` | Folder containing the individual analysis for this candidate. |
| `status` | `complete` when analysis and aggregation succeeded; otherwise `failed`. |
| `error` | Failure message when `status` is `failed`; blank for a successful analysis. |
| `substitution_count` | Number of amino-acid differences between mutant and parent Fv sequences. |
| `substitutions` | Semicolon-separated substitutions in structure numbering, such as `E99A;R100D`. |

### Rosetta energies

| Column | Meaning and interpretation |
|---|---|
| `interface_parent_REU` | Parent ligand-interface energy, calculated as bound minus separated partners. More negative is more favorable. |
| `interface_mutant_REU` | Mutant ligand-interface energy, calculated in the same way. More negative is more favorable. |
| `interface_delta_REU` | `interface_mutant_REU - interface_parent_REU`. Negative favors the mutant; positive favors the parent. This is the primary ranking metric. |
| `bound_total_delta_REU` | Mutant minus parent total Rosetta score for the complete bound complex. Negative is lower-energy, but absolute totals across different sequences must be treated cautiously. |
| `interface_rep_delta_REU` | Mutant minus parent change in interface `fa_rep`, the short-range steric-repulsion term. Positive means more steric repulsion/clashing; negative means less. |

### Structural similarity and ligand contacts

| Column | Meaning and interpretation |
|---|---|
| `global_fv_rmsd_A` | Backbone or CA RMSD for the complete Fv after the requested global fit. |
| `framework_rmsd_A` | Framework-region RMSD measured after that same global fit. |
| `parent_ligand_contacts` | Total Arpeggio atom-level contact annotations between antibody residues and ligand in the parent. This is not a count of unique residues or unique residue pairs. |
| `mutant_ligand_contacts` | Equivalent atom-level contact-annotation count for the mutant. |
| `parent_ligand_residue_pairs` | Number of unique antibody-residue/ligand pairs with a structural contact in the parent. Multiple atom-level annotations collapse to one pair. |
| `mutant_ligand_residue_pairs` | Equivalent unique residue-pair count for the mutant. |
| `parent_rosetta_ligand_hbonds` | Ligand hydrogen bonds recognized by Rosetta in the parent pose. |
| `mutant_rosetta_ligand_hbonds` | Ligand hydrogen bonds recognized by Rosetta in the mutant pose. |
| `parent_arpeggio_ligand_hbonds` | Strict geometric ligand hydrogen-bond annotations from Arpeggio for the parent. |
| `mutant_arpeggio_ligand_hbonds` | Strict geometric ligand hydrogen-bond annotations from Arpeggio for the mutant. Rosetta and Arpeggio counts may differ because their definitions differ. |
| `parent_plip_ligand_hbonds` | Protein-ligand hydrogen bonds independently detected by PLIP in the parent, using the validated explicit hydrogens without reprotonation. |
| `mutant_plip_ligand_hbonds` | Protein-ligand hydrogen bonds independently detected by PLIP in the mutant. PLIP, Arpeggio, and Rosetta use different definitions, so disagreement is reported rather than silently reconciled. |
| `parent_plip_pi_stacks` | Parent protein-ligand pi stacks detected by PLIP. |
| `mutant_plip_pi_stacks` | Mutant protein-ligand pi stacks detected by PLIP. A PLIP/Arpeggio disagreement should trigger direct geometry inspection. |
| `parent_arpeggio_pi_stack_detected` | Binary flag: `True` when Arpeggio reports at least one plane-plane contact between a protein aromatic ring and the estradiol A-ring in the parent. |
| `mutant_arpeggio_pi_stack_detected` | Equivalent binary Arpeggio flag for the mutant. This does not imply a two-sided sandwich or PLIP agreement. |
| `parent_arpeggio_pi_stack_sandwich` | Binary flag: `True` when the parent has Arpeggio A-ring plane-plane contacts on opposite faces of the estradiol A-ring. |
| `mutant_arpeggio_pi_stack_sandwich` | Equivalent opposite-face Arpeggio sandwich flag for the mutant. |
| `parent_plip_pi_stack_detected` | Binary flag derived from `parent_plip_pi_stacks > 0`. |
| `mutant_plip_pi_stack_detected` | Binary flag derived from `mutant_plip_pi_stacks > 0`. |
| `parent_arpeggio_plip_pi_stack_consensus` | `True` when both Arpeggio and PLIP independently detect at least one parent pi-stack. It does not by itself establish an opposite-face sandwich. |
| `mutant_arpeggio_plip_pi_stack_consensus` | Equivalent cross-method binary for the mutant. |
| `comparison_quality_warning` | `True` when the comparison has a large total-score mismatch (`abs(bound_total_delta_REU) > 100`) or a large repulsive-interface increase (`interface_rep_delta_REU > 10`). Such rows require inspection before interpretation. |

### Ligand-hydrogen validation

| Column | Meaning |
|---|---|
| `parent_ligand_validation_status` | Parent ligand hydrogen validation result, such as `complete_explicit_hydrogens`, `missing_hydrogens_rosetta_rebuild_required`, or `invalid`. |
| `parent_explicit_ligand_hydrogens` | Explicit hydrogen atoms found on the parent ligand in the input structure. |
| `parent_expected_ligand_hydrogens` | Hydrogens expected for the parent ligand according to the ligand template/parameters used by the workflow. |
| `mutant_ligand_validation_status` | Mutant ligand hydrogen validation result. |
| `mutant_explicit_ligand_hydrogens` | Explicit hydrogen atoms found on the mutant ligand in the input structure. |
| `mutant_expected_ligand_hydrogens` | Hydrogens expected for the mutant ligand. |

### CDR RMSDs

All six values are measured **after one shared global Fv alignment**. They are not independently fitted local RMSDs.

| Column | Meaning |
|---|---|
| `h_cdr1_rmsd_A` | Heavy-chain CDR1 RMSD. |
| `h_cdr2_rmsd_A` | Heavy-chain CDR2 RMSD. |
| `h_cdr3_rmsd_A` | Heavy-chain CDR3 RMSD. |
| `l_cdr1_rmsd_A` | Light-chain CDR1 RMSD. |
| `l_cdr2_rmsd_A` | Light-chain CDR2 RMSD. |
| `l_cdr3_rmsd_A` | Light-chain CDR3 RMSD. |

### Mutation and interaction-network summary

An interaction **pair-type event** is one unique residue pair plus one interaction class that is gained or lost near a mutated residue. A single residue pair can carry several classes. Shared mutation-mutation pairs are deduplicated at the candidate level.

| Column | Meaning and interpretation |
|---|---|
| `mutation_assessments` | Semicolon-separated compact assessment for every substitution. The detailed evidence is in `batch_mutation_summary.csv`. |
| `network_added_pair_type_events` | Number of unique residue-pair/interaction-class combinations present in the mutant but not the parent. |
| `network_removed_pair_type_events` | Number of unique residue-pair/interaction-class combinations present in the parent but not the mutant. |
| `network_added_unique_pairs` | Number of unique residue pairs with at least one added interaction class. |
| `network_removed_unique_pairs` | Number of unique residue pairs with at least one removed interaction class. A pair may be in both added and removed sets if its interaction type changed. |
| `network_added_interaction_classes` | Added classes and their pair counts, formatted as `class:count`; for example, `hbond:2;vdw:4`. |
| `network_removed_interaction_classes` | Removed classes and their pair counts in the same format. |
| `network_added_pairs_parent_energy_REU` | Sum of parent pair energies for all unique pairs that gained at least one class. Each pair is summed once. |
| `network_added_pairs_mutant_energy_REU` | Sum of mutant pair energies for those same added-class pairs. |
| `network_added_pairs_energy_delta_REU` | Mutant minus parent energy for the union of added-class pairs. Negative means that set of pairs became energetically more favorable. |
| `network_removed_pairs_parent_energy_REU` | Sum of parent pair energies for all unique pairs that lost at least one class. Each pair is summed once. |
| `network_removed_pairs_mutant_energy_REU` | Sum of mutant pair energies for those same removed-class pairs. |
| `network_removed_pairs_energy_delta_REU` | Mutant minus parent energy for the union of removed-class pairs. Negative means the affected pairs are more favorable overall despite losing one or more geometric classes; positive means they worsened. |
| `network_changed_unique_pairs` | Number of unique residue pairs in the union of added and removed sets. A pair is counted once even if it both gained and lost classes. |
| `network_changed_unique_pairs_energy_delta_REU` | Mutant-minus-parent summed energy for that deduplicated union. This is the safest batch-level network energy total because it avoids overlap between added and removed sets. |

The summary also includes formed/broken H-bond and salt-bridge fields named
`{event}_{scope}_{family}_{measure}`, where:

- `event` is `formed` or `broken`;
- `scope` is `intra_chain`, `inter_chain`, or `protein_ligand`;
- `family` is `hbond` or `salt_bridge`;
- `measure` is `pair_count`, `pairs`, or `energy_delta_REU`.

H-bond pairs are deduplicated when both Rosetta and Arpeggio annotate the same
pair. Salt bridges use Arpeggio `ionic` annotations and are protonation-dependent.

### Conserved estradiol-binding architecture screen

These fields implement the four conserved elements described for estrogen-binding proteins and antibodies: hydrogen bonding at the O3 and O17 ends, hydrophobic recognition of the C18 methyl, and a two-sided clamp around the aromatic A-ring. The paper's general "sandwich" is geometric and may contain aliphatic, aromatic, or charged side chains; it is not automatically a pi-stack sandwich. The strict project screen additionally requires two opposite-face Arpeggio pi-stack contacts.

Atom roles are read from the validated params/template hash-bound mapping, never guessed from names. In the current corrected E2L definition: `O1` = phenolic O3, `O2` = O17, `C10` = C18 methyl, and `C1,C2,C3,C4,C5,C6` = aromatic A-ring. Historical tables used a different mapping and architecture protocol; do not mix their flags with the current DS-only screen. See [E2 architecture](e2_architecture.md) for calibration and limitations.

| Column | Meaning and interpretation |
|---|---|
| `e2_o3_hbond` | `True` when DS reports at least one conventional H-bond involving the mapped O3 or its hydroxyl hydrogen. |
| `e2_o3_hbond_residues` | Residues contributing the O3 hydrogen bond(s), or `none`. |
| `e2_o17_hbond` | `True` when DS reports at least one conventional H-bond involving the mapped O17 or its hydroxyl hydrogen. |
| `e2_o17_hbond_residues` | Residues contributing the O17 hydrogen bond(s), or `none`. |
| `e2_both_end_hbonds` | `True` only when both the O3 and O17 ends have a detected hydrogen bond. |
| `e2_c18_hydrophobic_contact` | `True` when at least one unique residue has a native DS alkyl/pi-alkyl contact to the singleton mapped C18 endpoint. Not SASA or burial. |
| `e2_c18_hydrophobic_partner_count` | Unique DS C18 partner residues, unioned across alkyl and pi-alkyl. |
| `e2_c18_hydrophobic_residues` | Residues forming those C18 hydrophobic contacts, or `none`. |
| `e2_a_ring_pi_stack_count` | Unique protein residues with a DS stacked pi-pi interaction to the mapped A-ring. T-shaped and charged pi contacts remain distinct classes. |
| `e2_a_ring_pi_contact_count` | Unique DS A-ring pi-mediated partners across stacked/T-shaped pi-pi, pi-cation and pi-anion; not an energy sum or pass requirement. |
| `e2_a_ring_pi_stack_residues` | Identities of the A-ring pi-stacking residues, or `none`. |
| `e2_pi_stack_sandwich` | Historical diagnostic only; not used by the current four-feature screen. |
| `e2_geometric_sandwich` | Six-ray side-chain wall proxy: at least two near paths per face, at least one distinct-residue paired position, gaps 0–1.5 A, zero overlapping rays. |
| `e2_steric_alpha_coverage`, `e2_steric_beta_coverage` | Carbon-normal paths with qualifying near walls, each out of six. Not a surface-burial fraction. |
| `e2_steric_paired_coverage` | Same-carbon paths bounded by distinct residues on opposite faces, out of six. |
| `e2_steric_overlap_rays` | Negative-gap paths, out of twelve; flagged rather than rewarded. |
| `e2_conserved_architecture_score_0_to_4` | One point each for O3 H-bond, O17 H-bond, C18 hydrophobic recognition, and geometric A-ring sandwich. Higher indicates more of the literature-derived architecture is present. |
| `e2_strict_screen_pass` | Current all-four-feature flag: O3, O17, C18 and steric features in the same pose. No Arpeggio/PLIP input and no pi-contact prerequisite. Not proof of improved affinity. |
| `e2_metric_status` | `complete`, `failed`, or `not_available` for the architecture calculation. |
| `e2_metric_error` | Diagnostic message when the architecture metric could not be calculated. |

### Workflow bookkeeping

| Column | Meaning |
|---|---|
| `reused_existing_result` | `True` if a previously completed individual result was aggregated instead of recalculated. |
| `rank_by_interface_delta` | Rank among completed candidates by `interface_delta_REU`, with the most negative value ranked 1. |
| `rank_by_e2_screen` | Screening rank: strict-pass candidates first, then higher architecture score, then more favorable interface-energy delta. |

## `batch_interaction_network_summary.csv`

Each row groups one candidate, one event direction, and one interaction class. Use `contributing_pairs` to see exactly which residue pairs were formed or broken.

| Column | Meaning and interpretation |
|---|---|
| `candidate` | Mutant complex filename. |
| `result_folder` | Folder containing the individual candidate analysis. |
| `event` | `added` means the class formed in the mutant; `removed` means the class present in the parent was lost in the mutant. |
| `interaction_class` | Interaction category that formed or broke; see the glossary below. |
| `unique_pair_count` | Number of unique residue pairs contributing to this event/class row. |
| `contributing_pairs` | Semicolon-separated pair identities, such as `H:99<->X:1`. `X:1` is typically the ligand after Rosetta chain normalization. |
| `sum_parent_pair_energy_REU` | Sum of parent Rosetta pair energies for the listed pairs. |
| `sum_mutant_pair_energy_REU` | Sum of mutant Rosetta pair energies for the same listed pairs. |
| `sum_pair_energy_delta_REU` | `sum_mutant_pair_energy_REU - sum_parent_pair_energy_REU`; negative is energetically favorable. |
| `energy_overlap_warning` | Reminder that a multi-annotated pair appears in more than one interaction-class row. Therefore class-row energy sums must not be added together to obtain a grand total. Use `network_changed_unique_pairs_energy_delta_REU` for a deduplicated total. |

## `batch_hbond_saltbridge_changes.csv`

Each row is one unique formed or broken H-bond/salt-bridge residue pair near a
mutation. This is the readable source for intra-chain, inter-chain, and
protein-ligand bond changes.

| Column | Meaning and interpretation |
|---|---|
| `candidate` | Mutant complex filename. |
| `result_folder` | Individual analysis folder. |
| `event` | `formed` in the mutant or `broken` from the parent. |
| `interaction_family` | Deduplicated `hbond` or Arpeggio `salt_bridge`. |
| `scope` | `intra_chain`, `inter_chain`, or `protein_ligand`. |
| `pair` | Exact residue pair, such as `H:99<->L:91`. |
| `mutation` / `antibody_number` | Mutation through which this local-network pair was reported. |
| `source_interaction_types` | Underlying `hbond`, `rosetta_hbond`, or `ionic` annotations. |
| `parent_pair_energy_REU` | Parent Rosetta pair energy. |
| `mutant_pair_energy_REU` | Mutant Rosetta pair energy. |
| `pair_energy_delta_REU` | Mutant minus parent pair energy; negative is favorable. |
| `parent_pair_nominal_charge` / `mutant_pair_nominal_charge` / `delta_pair_nominal_charge` | Recorded nominal charge context for interpreting ionic contacts. |

## `batch_mutation_summary.csv`

Each row is one mutation in one candidate.

| Column | Meaning and interpretation |
|---|---|
| `candidate` | Mutant complex filename. |
| `result_folder` | Folder containing the individual candidate analysis. |
| `mutation` | Substitution in structure numbering: parent one-letter amino acid, residue number, mutant one-letter amino acid (for example, `E99A`). |
| `antibody_number` | ANARCII/AbNumber antibody position in the selected numbering scheme, such as `H107`. |
| `region` | Antibody region assigned from that numbering, such as `CDR3` or `FR3`. |
| `assessment` | Conservative rule-based interpretation of direct ligand and local-network evidence for this single pose. It is not experimental proof or a calculated single-mutant free energy. |
| `direct_ligand_evidence` | Parent-to-mutant Rosetta energy for the mutation-ligand pair, its delta, and parent-to-mutant atom-contact annotation counts. `No direct ligand edge detected` means neither contacts nor a meaningful pair energy was detected. |
| `local_network_delta_REU` | Sum of mutant-minus-parent Rosetta pair-energy changes between this mutated residue and its analyzed partners. Negative is favorable. The same mutation-mutation pair can appear in both mutations' rows, so do not sum this column across mutations as an independent total. |
| `strongest_supporting_pair_changes` | Up to two partner pairs that improved by at least 0.25 REU, with old/new interaction classes and energy delta. |
| `strongest_opposing_pair_changes` | Up to two partner pairs that worsened by at least 0.25 REU, with old/new classes and energy delta. |
| `chemistry_interpretation` | Whether the substitution retains the same primary amino-acid chemical class, followed by the parent and mutant classes. |
| `rotamer_interpretation` | Old-versus-new coarse torsion-bin result. HCDR3 uses anchor-only RMSD after global fitting at most the configured cutoff (default 3 A); other CDRs use their overall RMSD. Missing gate RMSD skips comparison. Framework mutations are not gated. |
| `cdr_global_fit_rmsd_A` | Overall RMSD of the CDR containing the mutation after global fitting. Retained for context, but not the HCDR3 rotamer gate. |
| `rotamer_gate_region` | In the substitution table: HCDR3 anchors for heavy-chain CDR3; the corresponding overall region for other CDRs; blank for framework mutations. |
| `rotamer_gate_rmsd_A` | In the substitution table: global-fit RMSD actually used for the rotamer cutoff. No whole-HCDR3 fallback if anchors are unavailable. |
| `combination_context` | Overall multi-mutant interface-energy change and a warning that this global score is non-additive and cannot be assigned to one mutation. |
| `evidence_limit` | Scope limitation: one relaxed pose, rather than a conformational ensemble, causal single-mutant scan, or experimental binding measurement. |
| `comparison_quality` | Human-readable quality-check result. A warning means the energy comparison should be inspected before biological interpretation. |

## Interaction-class glossary

| Class | Meaning |
|---|---|
| `hbond` | Conventional hydrogen bond identified by Arpeggio geometry. |
| `rosetta_hbond` | Hydrogen bond recognized/scored by Rosetta. |
| `weak_hbond` | Weaker geometric hydrogen bond, commonly involving a weak donor such as C-H. |
| `polar` | Polar atom contact satisfying Arpeggio's polar-contact definition. |
| `weak_polar` | Weaker or less stringent polar contact. |
| `ionic` | Oppositely charged ionic contact/salt bridge. |
| `hydrophobic` | Hydrophobic contact. |
| `vdw` | Favorable or permissible van der Waals contact annotation. |
| `vdw_clash` | Van der Waals overlap/clash annotation; generally unfavorable geometrically. |
| `pi_stack` | Aromatic pi-pi stacking interaction. |
| `aromatic_contact` | Aromatic-ring contact not necessarily satisfying the stricter pi-stacking geometry. |
| `carbon_pi` | Carbon-to-pi-system interaction. |
| `cation_pi` | Positively charged group interacting with an aromatic pi system. |
| `donor_pi` | Hydrogen-bond donor oriented toward an aromatic pi system. |
| `amidering` | Amide group interacting with an aromatic ring. |
| `carbonyl` | Carbonyl-associated contact classification. |
| `metsulphur_pi` | Methionine sulfur-to-aromatic-pi interaction. |
| `proximal` | Atoms are within Arpeggio's proximity cutoff, but this alone does not establish a specific favorable interaction. |
| `energy_only` | The residue pair crossed the configured Rosetta energetic presence/change threshold without an added or removed geometric class annotation. |
| `none` | No interaction class of the tracked types was detected. |

The individual candidate folder contains the pair-level tables used to produce these summaries. For mechanistic interpretation, examine both the interaction labels and the pair-energy delta rather than treating counts alone as proof of improved binding.
