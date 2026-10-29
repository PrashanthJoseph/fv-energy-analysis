# DS-only E2 recognition screen

This screen uses one saved OptH pose per complex. It does not generate an
ensemble, rerun Rosetta, or change coordinates. RosettaLigand scores and the
four independent tool fingerprint sections remain available separately.
Arpeggio and PLIP do not supply any of the E2 recognition flags.

## Features

1. **O3 and O17 hydrogen bonds:** Discovery Studio conventional H-bond records
   involving each chemically mapped oxygen, or its bonded hydroxyl hydrogen.
   Carbon H-bonds do not qualify. These operational flags do not require the
   same amino-acid identities in every protein family.
2. **C18 packing:** reference-calibrated DS hydrophobic contact recognition.
   At least one unique protein residue must have a DS alkyl or pi-alkyl contact
   whose ligand endpoint is the mapped C18 atom alone. A multi-carbon ligand
   group merely containing C18 is retained as a diagnostic, not counted as a
   C18-specific pass. The rule and reference-copy measurements are retained
   with the output.
   A DS alkyl/pi-alkyl annotation concerns atom groups: group-centroid distances
   must not be relabelled as C18-to-protein-atom distances. The latter are
   recorded separately. This flag is not SASA or percent burial.
3. **Steric sandwich (two-sided A-ring enclosure):** a coordinate calculation
   around the six chemically mapped A-ring carbons. It is not a native DS
   interaction type. Fit a least-squares plane through those carbons and search
   all protein side-chain heavy atoms. A qualifying atom must have a nearest
   A-ring-carbon van der Waals surface gap between 0 and 1.0 A, inclusive, and
   lie at least its own van der Waals radius from the plane. This excludes
   atoms confined to the ring edge. Require at least one pair of qualifying
   atoms from distinct residues on opposite faces whose joining segment
   crosses the plane within the projected union of the six ring-carbon van der
   Waals disks (carbon radius 1.70 A). The crossing need not lie inside the
   smaller polygon joining the carbon centres. Count each distinct residue pair
   once, however many atom pairs qualify. Backbone atoms, waters and hydrogens
   are excluded. Negative surface gaps are recorded separately as overlapping
   atom pairs; they do not qualify as favourable close contacts and do not
   invalidate a different nonoverlapping enclosure pair. Full contact atom,
   residue, distance, gap, face and pair measurements are saved.

   Surface gap is atom-centre distance minus the two van der Waals radii.
   Bondi radii are C 1.70, N 1.55, O 1.52 and S 1.80 A. Face labels use the
   chemically mapped C13-to-C18 direction when available; swapping the normal
   sign cannot change the binary result. The former six-carbon normal-ray
   definition is replaced by this contact-pair definition.

These are four features: O3, O17, C18 and steric geometry. An all-four pass means
all four coexist in the same pose. It is not proof of affinity or specificity;
the paper itself discusses an antibody lacking the entrance interaction.

## DS A-ring pi-mediated contacts

The union of DS stacked pi-pi, T-shaped pi-pi, pi-cation and pi-anion contacts
with the chemically mapped A-ring is counted by unique protein residue. Two
fused rings of one tryptophan therefore contribute one partner, not two.
Individual DS classes remain separate in the fingerprint section and detailed
JSON. Pi-alkyl is not silently relabelled pi-stacking. Pi-mediated contacts are
an additional descriptor, not a mandatory fifth feature or an energy estimate.

## Calibration and limitations

The source concept is Nahoum et al. (2003), *How estrogen-specific proteins
discriminate estrogens from androgens: A common steroid-binding site
architecture*, [doi:10.1096/fj.02-0524fje](https://doi.org/10.1096/fj.02-0524fje).
Figures 3–4 and the antibody discussion identify the reference pockets and C18
contacts. The antibody discussion supports at least one hydrophobic C18 partner;
it does not establish a universal two-residue requirement or a DS cutoff.

Seven representative systems and sixteen crystallographic bound copies are
used: 1IOL, 1ERE, 1AQU, 2BFV, 1JGL, 1JNH and 1JNN. Crystal copies are not
independent protein families. 2BFV and 1JNH have modified estrogens with an
unchanged C18 methyl; these are local-geometry references, not unmodified-E2
affinity controls. Reference hydrogen preparation is documented separately
and must not move heavy atoms. It is not an OptH replicate campaign.

For the enclosure rule, the nearest side-chain contact in each paper-named
wall group has a largest surface gap of 0.844 A across the sixteen examined
copies. The project rounds that allowance to 1.0 A. All sixteen copies have
named-wall pairs satisfying the new rule. The current m7 reference also
qualifies through Asp112H CB and Tyr49L CB, with gaps approximately 0.34 and
0.46 A. Both m7 and reference antibody 1JNN illustrate offset enclosure whose
contact segment crosses the carbon van der Waals footprint outside the
carbon-centre polygon. Perfect centring is therefore not required.

The 1.0 A allowance, face-height and footprint criteria are the project's
operational definition, not numerical cutoffs supplied by the paper. These
positive examples establish reference compatibility, not independent validation
or E2 specificity. This feature is not SASA, favourable packing energy, or a
requirement for two aromatic residues or pi stacks. The paper's proposed C19
exclusion is a separate mechanism: the enclosure flag does not test a
hypothetical C19 methyl position or establish androgen discrimination.

For C18, native DS with explicit-H reference preparation detects 1–4 unique
partners in all sixteen copies; a two-partner requirement misses 1JNN. Removing
all detected partner residues causes all sixteen synthetic controls to lose
the flag. The paper's at-least-one concept is therefore used, without adding a
fitted distance cutoff on top of DS's own classification. All reference heavy
coordinates remain unchanged. Heavy-only reference inputs gave no methyl DS
calls and are not used as production-negative calibration examples.

Native DS does not reproduce every paper-named C18 residue (for example HSD1
Val143 and 17E12E5 Glu99). It recognizes other hydrophobic C18 partners in
those pockets. The criterion is thus a DS operationalization of the paper's
local recognition concept, not exact reproduction of every figure annotation.

All steroid names are chemically mapped. In standard reference PDBs the A-ring
is C1,C2,C3,C4,C5,C10 and the methyl is C18. In this project's validated E2L
parameters the A-ring is file C1–C6 and the methyl is file C10 (bonded to file
C16). Literal file atom C18 is not steroid C18.

The paper addresses estrogen-versus-androgen recognition. These metrics do
not establish E2-versus-E3 discrimination, quantitative binding energy, or
statistical significance from one static pose.

## Refresh a completed batch

Run `python -m fv_energy_analysis.architecture_batch --help` for the single
command entry point. Supply the completed source batch, a new output folder,
the matching params/template, original one-pass OptH result folder and the
reference calibration JSON. Three reporting workers are the default; the
parent is checked once. Input hashes, all named coordinates (including H),
detector row/status associations and unchanged non-architecture metrics are
validated. Reuse also requires current architecture and steric-method versions,
so historical ray results are recalculated from saved coordinates rather than
silently reused. Errors retain the exact file, candidate, exception and traceback.

The historical DS run did not save its own input digest. Reuse is checked
through retained prepared inputs, completed detector associations and matching
OptH coordinates; it is not claimed to be cryptographically DS-input-bound.

Outputs are a flat full-metrics CSV/JSON, architecture-only CSV, parent JSON,
and individual architecture reports with links to the original full reports.
The `--xlsx` reporting option can create the matching sortable XLSX
without re-executing scientific analysis. Existing results are not overwritten.
