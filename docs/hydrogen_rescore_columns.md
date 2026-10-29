# Hydrogen-optimized RosettaLigand rescore

Each CSV row represents one mutant. Its parent values come from the same m7
reference, optimized once with the same hydrogen-only protocol. Original files
and production workbooks are not modified.

## Protocol

One standard optimize-H packing pass samples supported proton-chi rotamers across
the intact antibody–ligand complex. Heavy atoms, residue identities, protonation
states and histidine tautomers are fixed. There is no continuous minimization,
heavy-atom repacking or docking. A fixed seed makes results independent of worker
scheduling. The original is retained if the trial does not improve the bound
whole-complex energy. Authoritative final scores are calculated after reloading
the saved PDB, including its coordinate rounding.

The historical `ligand.wts` function and pre-Talaris settings, ligand parameters
and PyRosetta version are identical to the source rescore. Energies are weighted
Rosetta energy units (REU), not experimental binding free energies.

## Columns

| Column or family | Meaning |
| --- | --- |
| `candidate`, `substitutions` | Original filename and identified substitutions. |
| `status` | `complete` or explicit failure; never interpret missing failed scores as zero. |
| `optimized_pdb`, `parent_optimized_pdb` | Saved coordinates used for the final comparison. |
| `optimized_pdb_sha256` | Content hash of the exported mutant structure. |
| `weights` | `ligand`, the historical RosettaLigand weight set. |
| `interface_parent_REU`, `interface_mutant_REU` | Optimized bound-minus-rigidly-separated ligand interface energy. Hydrogens are not reoptimized in the separated state. |
| `interface_delta_REU` | Optimized mutant minus optimized parent interface energy; negative is more favorable within this scoring protocol. |
| `before_interface_*` | Original, unoptimized scores and their matched original-parent difference, reproduced from the previous rescore. |
| `hydrogen_optimization_interface_change_REU` | Mutant interface after hydrogen optimization minus its own original interface. It is not the mutant-versus-parent difference. |
| `parent_bound_total_REU`, `mutant_bound_total_REU` | Optimized whole-complex energies, including protein internal terms; not interface-only energies. |
| `before_parent_bound_total_REU`, `before_mutant_bound_total_REU` | Original whole-complex energies. |
| `parent_<term>_REU`, `mutant_<term>_REU`, `delta_<term>_REU` | Optimized weighted interface contribution, and mutant minus parent contribution, for the named term. |
| `parent_ligand_hbonds`, `mutant_ligand_hbonds` | Newly calculated scored Rosetta donor–acceptor H-bond records crossing the protein–ligand interface. |
| `before_parent_ligand_hbonds`, `before_mutant_ligand_hbonds` | Scored Rosetta H-bond counts from the original rescore. |
| `selected_trial` | `1`: optimized pass selected; `0`: original retained because the trial's bound score did not improve. |
| `optimization_seconds` | Time to build the hydrogen-only task, score and run the single optimization pass; excludes initialization and final exports. |
| `max_heavy_displacement_A` | Largest heavy-atom displacement after exporting and reloading the selected structure; must be at most 0.000001 Å. |
| `moved_hydrogen_count` | Number of hydrogens differing from the original loaded structure above the same coordinate tolerance. Includes minor PDB rounding effects. |
| `rank_by_interface_delta` | Ascending rank by optimized mutant-minus-parent interface energy, not an affinity probability. |

The exported interface terms are `fa_atr` (van der Waals attraction), `fa_rep`
(steric repulsion), `fa_sol` (implicit solvation), `fa_elec` (electrostatics),
`fa_pair` (historical residue-pair statistical term), `hbond_bb_sc` (backbone–side
chain H-bond energy) and `hbond_sc` (side-chain H-bond energy). All active terms
are retained in each case's `interface.sc`, including terms that are zero across
this interface.

Arpeggio, PLIP and Discovery Studio fingerprints are **not regenerated** by this
rescore. Their original H-bond or network counts must not be interpreted as
measurements on the new hydrogen-optimized structures.
