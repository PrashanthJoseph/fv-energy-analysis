"""Pure coordinate/state controls for whole-complex hydrogen optimization."""
from copy import deepcopy
import unittest

from fv_energy_analysis.hydrogen_benchmark import compare_snapshots


def _snapshot():
    # Include both antibody chains and ligand, including a fixed His tautomer.
    return {
        "residue_types": ["SER", "HIS_D", "THR", "E2L"],
        "atoms": {
            (1, "OG"): {"xyz": (0.0, 0.0, 0.0), "is_hydrogen": False, "residue": "H:35:SER"},
            (1, "HG"): {"xyz": (1.0, 0.0, 0.0), "is_hydrogen": True, "residue": "H:35:SER"},
            (2, "ND1"): {"xyz": (3.0, 0.0, 0.0), "is_hydrogen": False, "residue": "H:99:HIS"},
            (2, "HD1"): {"xyz": (4.0, 0.0, 0.0), "is_hydrogen": True, "residue": "H:99:HIS"},
            (3, "OG1"): {"xyz": (0.0, 3.0, 0.0), "is_hydrogen": False, "residue": "L:35:THR"},
            (3, "HG1"): {"xyz": (1.0, 3.0, 0.0), "is_hydrogen": True, "residue": "L:35:THR"},
            (4, "O2"): {"xyz": (0.0, 0.0, 3.0), "is_hydrogen": False, "residue": "X:1:E2L"},
            (4, "H2"): {"xyz": (1.0, 0.0, 3.0), "is_hydrogen": True, "residue": "X:1:E2L"},
        },
    }


class HydrogenBenchmarkSnapshotTests(unittest.TestCase):
    def test_hydrogens_can_move_in_both_protein_chains_and_ligand(self):
        before = _snapshot()
        after = deepcopy(before)
        after["atoms"][(1, "HG")]["xyz"] = (1.0, 0.3, 0.4)
        after["atoms"][(3, "HG1")]["xyz"] = (1.0, 3.1, 0.0)
        after["atoms"][(4, "H2")]["xyz"] = (1.0, 0.0, 3.2)
        original_before, original_after = deepcopy(before), deepcopy(after)

        result = compare_snapshots(before, after)

        self.assertEqual(result["heavy_atom_count"], 4)
        self.assertEqual(result["hydrogen_atom_count"], 4)
        self.assertEqual(result["max_heavy_displacement_A"], 0.0)
        self.assertAlmostEqual(result["max_hydrogen_displacement_A"], 0.5)
        self.assertEqual(result["moved_hydrogen_count"], 3)
        self.assertEqual(set(result["moved_hydrogen_residues"]), {"H:35:SER", "L:35:THR", "X:1:E2L"})
        self.assertEqual(before, original_before)
        self.assertEqual(after, original_after)

    def test_unchanged_coordinates_have_no_moved_hydrogens(self):
        before = _snapshot()
        result = compare_snapshots(before, deepcopy(before))
        self.assertEqual(result["moved_hydrogen_count"], 0)
        self.assertEqual(result["moved_hydrogen_residues"], [])
        self.assertEqual(result["max_hydrogen_displacement_A"], 0.0)

    def test_two_moved_atoms_of_same_residue_are_counted_once_as_residue(self):
        before = _snapshot()
        before["atoms"][(1, "H")] = {
            "xyz": (0.0, -1.0, 0.0), "is_hydrogen": True, "residue": "H:35:SER",
        }
        after = deepcopy(before)
        after["atoms"][(1, "HG")]["xyz"] = (1.0, 0.1, 0.0)
        after["atoms"][(1, "H")]["xyz"] = (0.1, -1.0, 0.0)
        result = compare_snapshots(before, after)
        self.assertEqual(result["moved_hydrogen_count"], 2)
        self.assertEqual(result["moved_hydrogen_residues"], ["H:35:SER"])

    def test_any_heavy_atom_motion_above_tolerance_is_rejected(self):
        for key in ((1, "OG"), (2, "ND1"), (3, "OG1"), (4, "O2")):
            with self.subTest(atom=key):
                before = _snapshot()
                after = deepcopy(before)
                x, y, z = after["atoms"][key]["xyz"]
                after["atoms"][key]["xyz"] = (x + 1e-4, y, z)
                with self.assertRaises(ValueError):
                    compare_snapshots(before, after, tolerance=1e-6)

    def test_numerical_noise_below_heavy_atom_tolerance_is_allowed(self):
        before = _snapshot()
        after = deepcopy(before)
        after["atoms"][(1, "OG")]["xyz"] = (1e-7, 0.0, 0.0)
        result = compare_snapshots(before, after, tolerance=1e-6)
        self.assertAlmostEqual(result["max_heavy_displacement_A"], 1e-7)

    def test_histidine_tautomer_or_residue_variant_change_is_rejected(self):
        before = _snapshot()
        for replacement in ("HIS", "HIS_P", "HIS_D:protonated"):
            with self.subTest(replacement=replacement):
                after = deepcopy(before)
                after["residue_types"][1] = replacement
                with self.assertRaises(ValueError):
                    compare_snapshots(before, after)

    def test_residue_addition_or_deletion_is_rejected(self):
        before = _snapshot()
        for types in (before["residue_types"][:-1], before["residue_types"] + ["HOH"]):
            with self.subTest(types=types):
                after = deepcopy(before)
                after["residue_types"] = types
                with self.assertRaises(ValueError):
                    compare_snapshots(before, after)

    def test_missing_or_added_atoms_are_rejected(self):
        before = _snapshot()
        missing = deepcopy(before)
        del missing["atoms"][(4, "H2")]
        added = deepcopy(before)
        added["atoms"][(4, "H_NEW")] = {
            "xyz": (0.0, 0.0, 4.0), "is_hydrogen": True, "residue": "X:1:E2L",
        }
        for after in (missing, added):
            with self.subTest(atom_count=len(after["atoms"])):
                with self.assertRaises(ValueError):
                    compare_snapshots(before, after)

    def test_renaming_a_hydrogen_without_changing_atom_count_is_rejected(self):
        before = _snapshot()
        after = deepcopy(before)
        after["atoms"][(2, "HE2")] = after["atoms"].pop((2, "HD1"))
        with self.assertRaises(ValueError):
            compare_snapshots(before, after)

    def test_hydrogen_heavy_atom_class_changes_are_rejected(self):
        for key in ((4, "H2"), (4, "O2")):
            with self.subTest(atom=key):
                before = _snapshot()
                after = deepcopy(before)
                after["atoms"][key]["is_hydrogen"] = not before["atoms"][key]["is_hydrogen"]
                with self.assertRaises(ValueError):
                    compare_snapshots(before, after)


if __name__ == "__main__":
    unittest.main()
