"""DS-only classification tests; all inputs are synthetic temporary artifacts."""
import csv
import json
import math
from pathlib import Path
from tempfile import TemporaryDirectory
import unittest

from fv_energy_analysis.ds_architecture import analyze_ds_architecture


ROLES = {"o3_atom": "O1", "o17_atom": "O2", "c18_methyl_atom": "C10",
         "a_ring_atoms": [f"C{i}" for i in range(1, 7)]}
CALIBRATION = {"production_rule": {"interaction_types": ["alkyl", "pi_alkyl"],
               "ligand_endpoint": "singleton_mapped_C18", "minimum_unique_residues": 1,
               "maximum_atom_distance_A": None}}


def line(serial, name, xyz, residue="E2L", chain="X", number=1, element="C"):
    x, y, z = xyz
    record = "HETATM" if residue == "E2L" else "ATOM"
    return (f"{record:<6}{serial:5d} {name:>4} {residue:>3} {chain}{number:4d}    "
            f"{x:8.3f}{y:8.3f}{z:8.3f}{1:6.2f}{0:6.2f}          {element:>2}\n")


class DSArchitectureTests(unittest.TestCase):
    def setUp(self):
        self.temporary = TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        root = Path(self.temporary.name)
        self.pdb, self.csv = root / "pose.pdb", root / "ds.csv"
        atoms = [(f"C{i+1}", (math.cos(i * math.pi / 3), math.sin(i * math.pi / 3), 0), "C") for i in range(6)]
        atoms += [("O1", (-4, 0, 0), "O"), ("H1", (-4.96, 0, 0), "H"),
                  ("O2", (4, 0, 0), "O"), ("H2", (4.96, 0, 0), "H"),
                  ("C10", (0, 0, 7), "C")]
        text = "".join(line(i + 1, name, xyz, element=element) for i, (name, xyz, element) in enumerate(atoms))
        protein = [("CB", (0, 0, 3.7), "H", "ALA", 1, "C"),
                   ("CB", (0, 0, -3.7), "L", "ALA", 2, "C"),
                   ("OE1", (-7, 0, 0), "L", "GLU", 55, "O"),
                   ("HG", (7, 0, 0), "H", "SER", 35, "H"),
                   ("OD1", (0, 3, 3), "H", "ASP", 112, "O"),
                   ("NZ", (0, -3, 3), "H", "LYS", 100, "N")]
        for i, (name, xyz, chain, residue, number, element) in enumerate(protein, 100):
            text += line(i, name, xyz, residue, chain, number, element)
        for i, name in enumerate(("CG", "CD1", "CD2", "NE1", "CE2", "CE3", "CZ2", "CZ3", "CH2"), 200):
            text += line(i, name, (i / 100, 0, 4), "TRP", "H", 99, "N" if name == "NE1" else "C")
        self.pdb.write_text(text, encoding="utf-8")
        self.rows = []

    def row(self, kind, protein, ligand):
        self.rows.append({"interaction_type": kind, "protein_atoms": protein, "ligand_atoms": ligand,
                          "distance_A": "4.0"})

    def run_analysis(self, **kwargs):
        with self.csv.open("w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(handle, fieldnames=["interaction_type", "protein_atoms", "ligand_atoms", "distance_A"])
            writer.writeheader()
            writer.writerows(self.rows)
        return analyze_ds_architecture(self.pdb, self.csv, roles=ROLES, c18_calibration=CALIBRATION, **kwargs)

    def ring(self):
        return ";".join(f"X:E2L1:C{i}" for i in range(1, 7))

    def test_conventional_both_h_donor_endpoint_and_oxygen_acceptor(self):
        self.row("conventional_hbond", "L:GLU55:OE1", "X:E2L1:H1")
        self.row("conventional_hbond", "H:SER35:HG", "X:E2L1:O2")
        result = self.run_analysis()
        self.assertTrue(result["e2_o3_hbond"])
        self.assertTrue(result["e2_o17_hbond"])
        self.assertTrue(result["e2_both_end_hbonds"])
        self.assertEqual(result["e2_hydroxyl_hydrogen_mapping"]["o3"]["hydrogen_distances_A"], {"H1": .96})

    def test_carbon_hbond_is_not_conventional(self):
        self.row("carbon_hbond", "L:GLU55:OE1", "X:E2L1:H1")
        self.assertFalse(self.run_analysis()["e2_o3_hbond"])

    def test_c18_requires_singleton_group_not_merely_membership(self):
        self.row("alkyl", "H:ALA1:CB", "X:E2L1:C10;X:E2L1:C1")
        result = self.run_analysis()
        self.assertFalse(result["e2_c18_hydrophobic_contact"])
        self.assertEqual(len(result["e2_matched_ds_records"]["c18_group_diagnostics"]), 1)
        self.row("pi_alkyl", "H:ALA1:CB", "X:E2L1:C10")
        result = self.run_analysis()
        self.assertTrue(result["e2_c18_hydrophobic_contact"])
        self.assertEqual(result["e2_c18_hydrophobic_partner_count"], 1)
        self.assertAlmostEqual(result["e2_matched_ds_records"]["c18"][0]["nearest_protein_heavy_distance_A"], 3.3)

    def test_native_c18_call_not_filtered_by_unpublished_distance(self):
        self.row("alkyl", "L:ALA2:CB", "X:E2L1:C10")
        result = self.run_analysis()
        self.assertTrue(result["e2_c18_hydrophobic_contact"])
        self.assertGreater(result["e2_matched_ds_records"]["c18"][0]["nearest_protein_heavy_distance_A"], 10)

    def test_pi_unique_residue_across_fused_rings_and_contact_classes(self):
        self.row("pi_pi_stacked", "H:TRP99:CG;H:TRP99:CD1;H:TRP99:NE1;H:TRP99:CE2;H:TRP99:CD2", self.ring())
        self.row("pi_pi_stacked", "H:TRP99:CD2;H:TRP99:CE2;H:TRP99:CZ2;H:TRP99:CH2;H:TRP99:CZ3;H:TRP99:CE3", self.ring())
        self.row("pi_pi_t_shaped", "H:TRP99:CG;H:TRP99:CD1;H:TRP99:NE1", self.ring())
        self.row("pi_anion", "H:ASP112:OD1", self.ring())
        self.row("pi_cation", "H:LYS100:NZ", self.ring())
        result = self.run_analysis()
        self.assertEqual(result["e2_a_ring_pi_stack_count"], 1)
        self.assertEqual(result["e2_a_ring_pi_stack_ring_count"], 2)
        self.assertEqual(result["e2_a_ring_pi_contact_count"], 3)
        self.assertEqual(result["e2_ds_pi_anion_partner_count"], 1)

    def test_incomplete_a_ring_and_pi_alkyl_not_pi_mediated(self):
        self.row("pi_cation", "H:LYS100:NZ", "X:E2L1:C1")
        self.row("pi_alkyl", "H:ALA1:CB", self.ring())
        self.assertEqual(self.run_analysis()["e2_a_ring_pi_contact_count"], 0)

    def test_four_feature_pass_does_not_require_any_pi_contact(self):
        self.row("conventional_hbond", "L:GLU55:OE1", "X:E2L1:H1")
        self.row("conventional_hbond", "H:SER35:HG", "X:E2L1:O2")
        self.row("alkyl", "H:ALA1:CB", "X:E2L1:C10")
        result = self.run_analysis()
        self.assertTrue(result["e2_geometric_sandwich"])
        self.assertTrue(result["e2_strict_screen_pass"])
        self.assertEqual(result["e2_conserved_architecture_score_0_to_4"], 4)
        self.assertEqual(result["e2_a_ring_pi_contact_count"], 0)
        self.assertEqual(result["e2_steric_method"], result["e2_steric_sandwich"]["method"])
        self.assertGreaterEqual(result["e2_steric_opposing_residue_pair_count"], 1)
        self.assertNotIn("e2_steric_overlap_rays", result)
        self.assertNotIn("rays", result["e2_steric_sandwich"])
        json.dumps(result)

    def test_missing_ds_returns_unknown_not_negative_interactions(self):
        result = analyze_ds_architecture(self.pdb, None, roles=ROLES, c18_calibration=CALIBRATION)
        for key in ("e2_o3_hbond", "e2_o17_hbond", "e2_both_end_hbonds", "e2_c18_hydrophobic_contact",
                    "e2_a_ring_pi_contact_count", "e2_strict_screen_pass"):
            self.assertIsNone(result[key])
        self.assertIsInstance(result["e2_geometric_sandwich"], bool)

    def test_wrong_csv_atom_or_ligand_rejected(self):
        for ligand in ("X:E2L2:O1", "X:E2L1:NOT_REAL"):
            self.rows = []
            self.row("conventional_hbond", "L:GLU55:OE1", ligand)
            with self.assertRaises(ValueError):
                self.run_analysis()

    def test_irrelevant_multiresidue_group_ignored_but_relevant_group_rejected(self):
        protein = "L:GLU55:OE1;H:SER35:HG"
        self.row("pi_amide", protein, self.ring())
        self.assertEqual(self.run_analysis()["e2_architecture_status"], "complete")
        self.row("pi_anion", protein, self.ring())
        with self.assertRaisesRegex(ValueError, "exactly one protein residue"):
            self.run_analysis()

    def test_missing_hydroxyl_hydrogen_fails_instead_of_silent_false(self):
        self.pdb.write_text("".join(row for row in self.pdb.read_text().splitlines(keepends=True) if row[12:16].strip() != "H1"))
        with self.assertRaisesRegex(ValueError, "exactly one explicit hydroxyl hydrogen"):
            self.run_analysis()

    def test_no_external_interaction_reader_or_rdkit_needed(self):
        import fv_energy_analysis.ds_architecture as module
        source = Path(module.__file__).read_text()
        self.assertNotIn("import rdkit", source)
        self.assertNotIn("import arpeggio", source)
        self.assertEqual(self.run_analysis()["e2_architecture_status"], "complete")

    def test_canonical_protein_hydrogen_alias_has_exact_distance_audit(self):
        self.pdb.write_text(self.pdb.read_text() + line(500, "2HD2", (-4, 0, 2), "ASN", "H", 99, "H"))
        self.row("conventional_hbond", "H:ASN99:HD22", "X:E2L1:O1")
        self.rows[-1]["distance_A"] = "2.0"
        result = self.run_analysis()
        self.assertTrue(result["e2_o3_hbond"])
        alias = result["e2_provenance"]["protein_hydrogen_atom_aliases"][0]
        self.assertEqual(alias["actual_pdb_atom"], "H:99:ASN:2HD2")
        self.assertEqual(result["e2_matched_ds_records"]["o3"][0]["protein_hydrogen_alias_distance_check"]["actual_pdb_distance_A"], 2)
        self.rows[-1]["distance_A"] = "2.01"
        with self.assertRaisesRegex(ValueError, "alias distance disagrees"):
            self.run_analysis()

    def test_protein_hydrogen_alias_never_applies_to_heavy_atom_or_ligand(self):
        self.pdb.write_text(self.pdb.read_text() + line(500, "2HD2", (-4, 0, 2), "ASN", "H", 99, "N"))
        self.row("conventional_hbond", "H:ASN99:HD22", "X:E2L1:O1")
        self.rows[-1]["distance_A"] = "2.0"
        with self.assertRaisesRegex(ValueError, "atoms missing"):
            self.run_analysis()

    def test_exact_protein_hydrogen_name_preferred_over_alias(self):
        self.pdb.write_text(self.pdb.read_text() + line(500, "2HD2", (-4, 0, 3), "ASN", "H", 99, "H")
                            + line(501, "HD22", (-4, 0, 2), "ASN", "H", 99, "H"))
        self.row("conventional_hbond", "H:ASN99:HD22", "X:E2L1:O1")
        self.rows[-1]["distance_A"] = "2.0"
        self.assertEqual(self.run_analysis()["e2_provenance"]["protein_hydrogen_atom_aliases"], [])


if __name__ == "__main__":
    unittest.main()
