# Inputs and configuration controls

Run `fv-energy-full CONFIG.json`. Relative paths resolve beside that JSON file.
`fv-energy-full --describe-options` lists inputs, defaults and switches without
starting a job. `fv-energy-full CONFIG.json --preflight` checks configuration
and dependencies in the interpreters that will actually run the stages.

## Accepted inputs

| Field | Accepted input |
|---|---|
| `parent` | One parent Fv–ligand complex |
| `input_folder` | Folder of candidates, default pattern `*.pdb` |
| `mutant` | One candidate file instead of `input_folder`; selects only that exact file |
| `output` | Output folder for individual analyses and master results |
| `params`, `ligand_template_pdb` | Matching Rosetta ligand parameters and reference PDB, required together |
| `cif_folder`, `parent_cif` | Optional matching original ModelCIF sources for H97–H113 pLDDT |
| `scores_csv` | Optional matching Boltz-2/GNINA table, or a native Boltz affinity table |
| `ligand_source`, `ligand_smiles` | Optional authoritative ligand identity; SDF/MOL2 is a ligand definition, not a complex |

The full OptH stage accepts PDB-format complexes. With
`optimize_hydrogens: false`, the scientific workflow also accepts `.cif` and
`.mmcif` complexes through Gemmi; set an appropriate `input_glob`. Convert
mmCIF to PDB before using the full OptH stage. Extract ZIP archives first.
Fv chains must have equal polymer lengths with substitutions only, and parent
and mutant must use the same ligand identity and validated chemistry. Original
confidence and external scores are joined to the exact matching model.

## Scoring and hydrogen preparation

| Knob | Default | Meaning |
|---|---|---|
| `weights` | `ligand` | Installed Rosetta weight-set name or custom `.wts` file |
| `legacy_pre_talaris` | `null` | Automatic for built-in ligand weights; explicit true/false required for a custom file |
| `optimize_hydrogens` | `true` | One full-complex OptH pass, preserving heavy atoms and chemical states |
| `hydrogen_weights` | Same as `weights` | Optional separate OptH weights, recorded independently |
| `hydrogen_legacy_pre_talaris` | Inherits scoring override | Historical energy-method setting for OptH |
| `hydrogen_seed` | `20261022` | Reproducible OptH seed |

Use `"optimize_hydrogens": true` to optimize the parent and every candidate
before scoring and interaction detection. All downstream analyses consume the
saved optimized structures. Use `false` for direct analysis of the supplied
coordinates, including PDBs optimized earlier. The Jupyter notebook exposes the
same toggle in `two_complex_settings.json`, with a default of `true`. Changing
it requires a new output folder (`run_root` in the notebook); completed results
are not silently reused across these two modes.

Examples are `ligand`, `ref2015`, or `D:/project/custom.wts`. For a custom
file derived from historical ligand weights, explicitly set
`legacy_pre_talaris: true`; use `false` for a modern function. Weight files
are hashed. A changed score function requires recalculation in a new output,
not numerical reweighting of an old CSV. These options apply to the full
runner; not every historical standalone rescore CLI accepts custom files.

## Which analyses run

Put the following booleans inside `analyses`. Ligand validation, Rosetta
pair/interface scoring and H-bonds, antibody numbering, and the structural
comparison remain the mandatory core of this energy-analysis workflow.

| Key | Default | Controls |
|---|---|---|
| `arpeggio` | `true` | Geometric contact classification |
| `plip` | `true` | Independent protein–ligand fingerprint |
| `discovery_studio` | `true` | Native Windows DS fingerprint |
| `e2_architecture` | DS enabled | DS O3/O17/C18/pi features and current enclosure geometry |
| `e3_discrimination` | `false` | Added C16α-OH fixed-pocket steric screen, without E3 energies |
| `network` | `true` | Mutation-centered formed/broken interaction networks |
| `rotamers` | `true` | Coarse torsion bins, subject to the RMSD gate |
| `hcdr3_plddt` | `true` | Confidence from matching original CIFs |
| `boltz2_scores` | `true` | Supplied Boltz-2 scores |
| `gnina_scores` | `true` | Supplied GNINA scores |

Disabled analyses are marked `not_run`, with blank unavailable values rather
than measured zero counts. Missing metadata sources are also unavailable.
With Arpeggio disabled, networks can still use Rosetta energies/H-bonds,
but Arpeggio-derived ionic/aromatic annotations are unavailable.

E2 requires DS. Disabling DS defaults E2 off; explicitly enabling E2 with DS
off is rejected. The integrated E3 report currently requires the completed
DS-based E2 report; standalone E3 supports a separate DS-free diagnostic.
The E2/E3 screens are estradiol-specific: disable them for another ligand.

## Structural, network and execution knobs

| Knob | Default | Choices / meaning |
|---|---|---|
| `scheme` | `imgt` | `imgt`, `chothia`, `kabat`, `aho` |
| `cdr_definition` | `imgt` | `imgt`, `chothia`, `kabat`, `north` |
| `fit_region` | `fv` | Whole Fv or `framework` fit |
| `rmsd_atoms` | `backbone` | N/CA/C/O or `ca` |
| `rotamer_rmsd_cutoff` | `3.0` Å | Global-frame gate; HCDR3 uses anchors |
| `hcdr3_tip_sequence` | `YYGSS` | Unique parent tip motif; all other HCDR3 residues are anchors |
| `network_energy_presence_cutoff` | `0.1` REU | Pair-energy presence threshold |
| `network_energy_change_cutoff` | `0.25` REU | Changed-pair threshold |
| `arpeggio_scope` | `neighbourhood` | Mutation neighbourhood or `full` |
| `arpeggio_neighbourhood_radius` | `6.0` Å | Union of parent/mutant neighbourhoods |
| `e3_overlap_warning` | `0.4` Å | Project steric-overlap warning threshold |
| `jobs` | `2` | Persistent workers; production example uses 3 |
| `input_glob` | `*.pdb` | Candidate filename pattern |
| `recursive` | `false` | Include nested folders |
| `stop_on_error` | `false` | Stop after a candidate failure |
| `include_empty_metadata_sections` | `true` | Retain blank requested metadata sections without sources |

Chain controls are `fv_chains`, `parent_chain_map`, `mutant_chain_map`,
`ligand_chain`, `parent_ligand_chain`, `mutant_ligand_chain` and `ligand_resname`.
Examples: `"fv_chains": ["H", "L"]`, `"mutant_chain_map": "A:H,B:L"`.
OptH retains protein chain labels; the optimized ligand is normalized to X.

Compatible results resume from their caches. Changed scientific options need
a new output folder. Report selection alone does not require new analysis.

## Which results appear in CSV/XLSX

| Knob | Effect |
|---|---|
| `report_sections` | Whole sections to show; `null` selects all |
| `report_columns` | Additional column keys; a union with the selected sections |
| `report_formats` | `["csv"]`, `["xlsx"]`, or both |
| `report_output` | Optional destination for the selected view |

Sections are Identity, Rosetta/RosettaLigand, Arpeggio, PLIP, Discovery Studio,
Structure & network, E2 architecture, E3 discrimination, HCDR3 pLDDT, Boltz-2,
GNINA and Analysis status, when available. Aliases include `rosetta`, `ds`,
`structure`, `e2`, `e3`, `plddt`. Identity always remains, with one candidate
per row. Unknown names and column keys are rejected before a full run starts.

Explicit selection creates one view in `selected_report/`. Complete canonical
scientific JSON/CSV remain for provenance, detail and reuse. The selected view
is bound to those source values. Hiding a section does not disable its analysis;
showing a section cannot make a disabled analysis computed.

For example, add to the full-run JSON:

```json
{
  "report_sections": ["rosetta", "structure", "e2", "plddt", "Boltz-2", "GNINA"],
  "report_formats": ["csv", "xlsx"]
}
```

Or change only the view of existing results:

```bash
fv-workbook BATCH --sections "rosetta,e2,plddt,Boltz-2,GNINA" --formats csv,xlsx
fv-workbook BATCH --columns "rosetta_interface_delta_REU,h_cdr3_anchors_rmsd_A,e2_strict_screen_pass" --formats csv
```

See [the editable example](../examples/full_pipeline_windows.json) for all
analysis switches and the standard Windows/WSL input layout.
