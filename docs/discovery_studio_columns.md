# Discovery Studio batch-output columns

The optional Windows-native Discovery Studio stage runs BIOVIA Discovery
Studio Visualizer's `LigandNonbondMonitor` on the same prepared PDB structures
used by the main workflow. It sets `AddHydrogens => False`, so BIOVIA reads the
validated explicit coordinates without adding another hydrogen model.

## Files

- `discovery_studio_summary.csv`: one parent-versus-mutant summary row per candidate.
- `discovery_studio_interaction_changes.csv`: exact formed and broken atom/group interactions for all candidates.
- `discovery_studio_parent_interactions.csv/json`: parent interaction records calculated once.
- `discovery_studio_status.csv`: one status row for the parent and every candidate.
- `batch_summary_with_discovery_studio.csv`: the main batch table augmented with Discovery Studio columns.
- Each candidate folder receives `discovery_studio_interactions.csv/json` and
  `discovery_studio_interaction_changes.csv`.

## Naming convention

For every Discovery Studio interaction class, the summary contains:

- `parent_ds_<class>_annotations` and `mutant_ds_<class>_annotations`: BIOVIA
  nonbond records of that class.
- `delta_ds_<class>_annotations`: mutant minus parent annotation count.
- `parent_ds_<class>_unique_residues` and
  `mutant_ds_<class>_unique_residues`: unique protein residue partners.
- `parent_ds_<class>_unique_atom_group_pairs` and
  `mutant_ds_<class>_unique_atom_group_pairs`: unique protein-group/ligand-group pairs.

Classes include conventional and carbon H-bonds, salt bridges, attractive
charge, pi-cation, pi-anion, pi-pi stacked, pi-pi T-shaped, pi-alkyl, pi-sigma,
donor-pi, lone-pair-pi, pi-sulfur, pi-amide, alkyl, halogen, metal and
unfavorable-contact categories.

`ds_formed_atom_group_interactions` and `ds_broken_atom_group_interactions`
count exact atom/group interaction records present on only one side of the
comparison. The detailed changes table gives the type, protein residue, both
atom groups and the BIOVIA distance.

## Cross-detector fields

- `mutant_ds_pi_cation_detected`, `mutant_arpeggio_pi_cation_detected`, and
  `mutant_plip_pi_cation_detected` are detector-specific binary results.
- `mutant_pi_cation_positive_detectors` lists the positive detectors.
- `mutant_pi_cation_detector_count` is the number of positive detectors.
- `mutant_pi_cation_consensus_2plus` requires at least two independent detectors.
- `mutant_ds_pi_anion_detected` is Discovery Studio-specific because PLIP and
  PDBe Arpeggio do not expose a corresponding standard interaction class.

These are geometric classifications, not energy estimates. Interpret them
alongside the Rosetta pair and interface-energy terms.
