# E3 discrimination: fixed-pocket steric assessment

Estriol (E3) differs from estradiol (E2) by an additional **16alpha-OH**.
Its C18 methyl group is unchanged. C18 burial or hydrophobic contacts therefore
cannot establish whether the extra estriol substituent fits.
See the [NCBI estriol definition](https://www.ncbi.nlm.nih.gov/mesh/68004964).

## Literature precedent and scope

The [human SHBG structural study by Grishkovskaya et al. (2000)](https://pmc.ncbi.nlm.nih.gov/articles/PMC305588/)
discusses reduced estriol affinity and proposes steric hindrance between its
16alpha-OH and the Asp65-Trp66 peptide bond. This is a mechanism proposed for
SHBG, not a validated m7 specificity criterion. It motivates inspecting the
actual added oxygen in its correct stereochemical direction.

The present section assesses **fixed-pocket sterics with explicit H-bond
geometry checks**. Energetic discrimination, local accommodation and solvent penalties
are not computed in this stage. A steric warning is not proof that E3 cannot
bind, and a warning-free pose is not proof of cross-reactivity.

## Operational definition

Construct estriol with the correct 16alpha stereochemistry while retaining the
protein and all shared E2/E3 heavy-atom coordinates. Construct the **neutral,
fully hydrogenated E3 ligand**: C18H24O3, comprising 21 heavy atoms and 24
hydrogens (45 ligand atoms in total). All source protein hydrogens remain
present with unchanged coordinates. Validate the molecular
identity and shared-atom correspondence before measuring its added oxygen.
Standard steroid numbering must come from chemical mapping; a file atom named
`C16` is not automatically steroid C16.

Inspect the added oxygen against all selected protein heavy atoms, including
backbone atoms. For each pair record:

- Protein atom and residue identity, element and backbone/side-chain status.
- Atom-centre distance in angstroms.
- Van der Waals surface gap: distance minus the sum of atomic radii.
- Overlap: the sum of atomic radii minus the distance, when positive.

The scientific result retains its exact radius set and overlap threshold.
The current default for a **large overlap** is at least **0.4 A**. This is a
project screening choice, not a cutoff published in the SHBG paper. Count each
protein residue once even when several of its atoms meet the criterion.

A large overlap against carbon or a chemically incompatible donor/acceptor
pair raises a steric warning. **N/O/S atoms are not automatically exempt.**
For a polar overlap to be counted as H-bond-supported, its chemically typed
donor/acceptor roles and explicit donor-H-acceptor geometry must satisfy the
recorded geometric criteria. A polar overlap that was fully tested and fails
the criteria contributes to the non-Hbond overlap category. Missing atom
typing or hydrogen geometry is retained as unassessed, not silently called a
clash or an H-bond. Protein backbone atoms are included.

The detailed report retains the atom identities, donor-hydrogen-acceptor roles,
donor-acceptor and hydrogen-acceptor distances, donor-H-acceptor angle and the
reason for accepting, rejecting or leaving a contact unassessed. The recorded
scientific protocol gives the exact cutoffs and their source; these are
geometric criteria, not measured hydrogen-bond energies.

The current conventional H-bond geometry requires donor-acceptor distance
**<4.1 A** and donor-H-acceptor angle **>100 degrees**, following the PLIP
geometric convention ([official thresholds](https://github.com/pharmai/plip/blob/master/plip/basic/config.py),
[official detection implementation](https://github.com/pharmai/plip/blob/master/plip/structure/detection.py)).
An additional project safeguard requires
hydrogen-acceptor distance **>=1.2 A**, so an unrealistically short contact is
not excused as a hydrogen bond. The safeguard is not a PLIP threshold.
Chemical donor/acceptor eligibility is checked first. These custom checks
are labelled project geometry, not native DS or Rosetta interaction counts.

The warning cannot account for side-chain changes, ligand movement or changes
to steroid-ring conformation that may accommodate E3. It is a diagnostic of
the specified fixed pose, not a calibrated selectivity prediction.

## Added hydroxyl hydrogen

Oxygen-to-protein-heavy-atom distances and overlap sizes are independent of the
orientation of the hydroxyl hydrogen. Polar-contact classification additionally
consults explicit-H geometry. The hydrogen diagnostic samples
**12 evenly spaced OH torsions** to avoid treating a clash in one arbitrarily
chosen OH-H orientation as decisive.

For each sampled orientation, determine the maximum non-H-bond-compatible
overlap of the added OH hydrogen with protein atoms. Report the **minimum of
those 12 maxima** as `e3_oh_h_best_sampled_protein_overlap_A`, displayed as
"Best OH-H nonpolar protein overlap (A)". This is a finite torsion sample,
not the minimum over every possible continuous angle. Intraligand contacts
contribute to orientation selection but are excluded from this protein-only
summary. This diagnostic can refer to a different sampled orientation from the
one selected to classify the complete oxygen-overlap network. Details retain the sampled orientations,
contact partners and radius conventions. The added polar hydrogen uses an
idealized **1.0 A radius**, a project diagnostic choice. This calculation is
not a MolProbity clashscore, PyRosetta OptH, repacking or energy minimization.
An overlap by itself does not assign an H-bond; the separate geometric test
must support it. No hydrogen-bond energy is calculated.

If the optional hydrogen diagnostic has not been assessed, its value remains
blank rather than zero; the primary oxygen assessment can still be complete.

The entire oxygen-overlap network is classified in **one consistent H25
orientation**. Among saved candidates the selection minimizes definite
non-Hbond overlaps, then unassessed contacts, and prefers the previously
persisted orientation on ties. Contacts supported only in different,
mutually incompatible hydroxyl orientations are not combined into one
apparently clash-free network. The exported E3 probe updates only H25 to
the selected saved orientation, matching the geometry in the report. All
protein atoms and other ligand atoms remain unchanged.

## Summary section

The optional, contiguous **E3 discrimination** section contains one row per
candidate and these independently sortable columns:

| Column | Meaning |
|---|---|
| E3 assessment status | Complete, failed, not assessed, or excluded because the input failed the existing chemistry gate. |
| 16alpha-O nearest distance (A) | Closest protein-heavy-atom distance to the added oxygen. |
| 16alpha-O nearest residue | Residue containing that closest atom. A neighbour identity is not a clash verdict. |
| Maximum vdW overlap (A) | Largest measured overlap, with polar contacts still requiring interpretation. |
| Large-overlap residues | Unique residues with any protein atom meeting the configured overlap threshold. |
| Non-Hbond overlap residues | Unique residues with a large overlap that lacks H-bond support after chemical typing and geometry assessment; includes incompatible or geometry-failing polar contacts. |
| Hbond-supported overlap residues | Unique residues whose large overlaps have chemically valid, geometrically supported H-bond explanations. |
| Unassessed overlap residues | Unique residues with large polar overlaps whose atom typing or required explicit-H geometry remains unresolved. |
| Large-overlap residue identities | Actual partner residues associated with the overlap count. |
| Fixed-pocket assessment | Steric overlap detected, H-bond-supported overlaps only, no large overlap in the fixed pose, or unassessed polar geometry. |
| Heuristic steric warning | TRUE for a definite non-Hbond-supported overlap; FALSE for a fully assessed pose without one; blank when geometry remains unassessed without a definite clash. Not an E3 specificity pass/fail. |
| Best OH-H nonpolar protein overlap (A) | Smallest maximum non-H-bond-compatible protein overlap across the 12 sampled OH-H torsions; intraligand contacts are excluded from this summary. |

The categorical assessment accompanies the three-state warning. `FALSE` is
not an E3-specificity conclusion, and an unassessed polar overlap cannot be
represented by a reassuring `FALSE`. Category counts refer to unique residues;
one residue can contain different contacts in more than one category.

Complete per-atom contacts, protein/ligand coordinates, chemical construction,
validation and method provenance remain in each candidate's
`e3_discrimination.json` and detailed report. Missing measurements are blank
in the table and `null` in JSON, never numerical zero or `FALSE`.

## Reporting and preservation of existing results

`fv_energy_analysis.e3_reports.append_e3_to_payload` adds the section to a
**deep copy** of an existing grouped payload. It preserves candidate identity,
row order, ranks, all non-E3 values and the separate four-feature E2 score. It
accepts saved results indexed by candidate plus a separately assessed parent.
The parent is not an additional mutant.
When upgrading an older E3 section, obsolete E3-only columns are replaced by
the geometry-resolved descriptors; they are not retained as hidden stale fields.

The standard fingerprint generator can include saved E3 artifacts only when
its explicit `include_e3_discrimination=True` option is supplied. Reporting
loads existing `e3_discrimination.json` files; it never constructs ligands,
measures contacts, runs Rosetta or launches a detector. A new output export
should be used for adding the section to an established production dataset.
Existing workbooks remain unchanged unless a new export is explicitly made.

## Run the module

Run these examples in the installed Linux analysis environment from the project
root, with this package installed or `src` on `PYTHONPATH`. Paths below are
placeholders; use the validated E2 parameters and matching ligand template
from the source analysis. The module performs no Rosetta rescoring.

One prepared E2 complex:

```bash
python -m fv_energy_analysis.e3_discrimination \
  --pdb /path/to/parent_m7.pdb \
  --params /path/to/E2L.params \
  --template /path/to/E2L_0001.pdb \
  --output /path/to/new_e3_parent_output
```

A folder of prepared E2 complexes:

```bash
python -m fv_energy_analysis.e3_discrimination \
  --input-folder /path/to/optimized_complexes \
  --params /path/to/E2L.params \
  --template /path/to/E2L_0001.pdb \
  --output /path/to/new_e3_folder_output \
  --jobs 3
```

Add the section to a **new copy** of a completed combined dataset:

```bash
python -m fv_energy_analysis.e3_discrimination \
  --source /path/to/completed_combined_dataset \
  --params /path/to/E2L.params \
  --template /path/to/E2L_0001.pdb \
  --output /path/to/new_combined_dataset_with_e3 \
  --jobs 3
```

Source mode follows each candidate's saved prepared-structure provenance,
preserves the original rows and ranks, and leaves prior chemistry rejections
explicitly excluded. The resulting grouped payload can be passed to the
established formatted-workbook exporter without running the scientific tools
again. `--limit` selects a small validation panel; unassessed rows are never
assigned zero values. `--resume` reuses compatible completed E3 artifacts.
Run `python -m fv_energy_analysis.e3_discrimination --help` for current options.

## Future energetic stage

A separate future assessment could compare matched E2/E3 interface energies
after consistent ligand preparation and local accommodation, and examine
the energies of any compensating polar interactions made by the extra OH. The general
energetic motivation for inspecting buried unsatisfied groups is described by
[Coventry and Baker (2021)](https://pmc.ncbi.nlm.nih.gov/articles/PMC7971855/).
Those energetic mechanisms are deliberately not inferred from the fixed-pocket
steric descriptors in the present section.
