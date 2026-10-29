"""Exact local-occluder SASA checks, including non-saturated controls."""
from dataclasses import replace
from pathlib import Path
import math
import tempfile
import unittest

import numpy as np

from fv_energy_analysis.c18_burial import (
    PROBE_RADIUS_A, _bio_structure, _burial_from_atoms, analyze_c18_burial,
)
from fv_energy_analysis.steric_sandwich import BONDI_RADII_A, StericAtom, read_steric_atoms


ROOT = Path(__file__).resolve().parents[1]
REAL_CASES = (
    ROOT / "results/2k_m7_hydrogen_optimized_rosettaligand_20261022/hydrogen_optimized_structures/parent_m7.pdb",
    ROOT / "results/2k_m7_hydrogen_optimized_full_rosettaligand_20261022/set_A_3mut_type2_contiguous_0_anchored_1_E1F_R2E_G3Y_estradiol_model_0_relaxed/prepared/mutant_input_rosetta.pdb",
)


def atom(name, xyz, *, chain="H", number=1, element="C", record="ATOM", residue="ALA"):
    return StericAtom(name, residue, chain, number, "", element,
                     np.array(xyz, dtype=float), record=record)


def ligand():
    # One self-occluding neighbor plus 18 remote ligand atoms. The deliberately
    # simple controls do not claim to represent chemically valid estradiol.
    return [atom("C10", (0, 0, 0), chain="X", record="HETATM", residue="E2L"),
            atom("C1", (-2.5, 0, 0), chain="X", record="HETATM", residue="E2L")] + [
        atom(f"Z{i}", (-20 - i * 8, 0, 0), chain="X", record="HETATM", residue="E2L")
        for i in range(18)
    ]


def full_sasa(lig, protein, n_points):
    from Bio.PDB.SASA import ShrakeRupley
    algorithm = ShrakeRupley(probe_radius=PROBE_RADIUS_A, n_points=n_points,
                            radii_dict=dict(BONDI_RADII_A))
    free_structure, free = _bio_structure(lig)
    bound_structure, bound = _bio_structure(list(lig) + list(protein))
    algorithm.compute(free_structure, level="A")
    algorithm.compute(bound_structure, level="A")
    index = next(i for i, item in enumerate(lig) if item.name == "C10")
    return float(free[index].sasa), float(bound[index].sasa)


def pdb_line(serial, item):
    x, y, z = item.xyz
    return (f"{item.record:<6}{serial:5d} {item.name:>4} {item.residue:>3} {item.chain:1}"
            f"{item.number:4d}    {x:8.3f}{y:8.3f}{z:8.3f}"
            f"{item.occupancy:6.2f}{0:6.2f}          {item.element:>2}\n")


class C18BurialTests(unittest.TestCase):
    def assert_full_equivalence(self, lig, protein, n_points):
        result = _burial_from_atoms(lig, protein, n_points=n_points)
        free, bound = full_sasa(lig, protein, n_points)
        self.assertEqual(result["isolated_full_ligand_SASA_A2"], free)
        self.assertEqual(result["complex_SASA_A2"], bound)
        self.assertAlmostEqual(result["buried_SASA_A2"], free - bound)
        return result

    def test_partial_exposure_local_equals_full_at_both_resolutions(self):
        protein = [atom("CA", (3.7, 0, 0)), atom("O", (0, 4.2, 0), element="O"),
                   atom("CB", (8, 0, 0), number=2), atom("CB", (30, 0, 0), number=3)]
        for n in (960, 1920):
            result = self.assert_full_equivalence(ligand(), protein, n)
            self.assertGreater(result["fraction_buried"], 0)
            self.assertLess(result["fraction_buried"], 1)
            self.assertEqual(result["protein_occluder_atom_count"], 2)
            self.assertTrue(result["protein_backbone_included"])

    def test_absent_or_distant_protein_gives_zero_burial(self):
        for protein in ([], [atom("CB", (6.201, 0, 0))]):
            result = self.assert_full_equivalence(ligand(), protein, 960)
            self.assertEqual(result["fraction_buried"], 0)
            self.assertEqual(result["protein_occluder_atom_count"], 0)

    def test_occluder_deletion_reduces_burial(self):
        protein = [atom("CB", (3.5, 0, 0)), atom("CB", (0, 0, 3.5), chain="L")]
        for n in (960, 1920):
            both = self.assert_full_equivalence(ligand(), protein, n)
            single = self.assert_full_equivalence(ligand(), protein[:1], n)
            absent = self.assert_full_equivalence(ligand(), [], n)
            self.assertGreater(both["fraction_buried"], single["fraction_buried"])
            self.assertGreater(single["fraction_buried"], absent["fraction_buried"])

    def test_geometric_boundary_keeps_possible_occluder(self):
        result = _burial_from_atoms(ligand(), [atom("CB", (6.2, 0, 0))], n_points=960)
        self.assertEqual(result["protein_occluder_atom_count"], 1)

    def test_no_coordinates_are_changed(self):
        lig, protein = ligand(), [atom("CA", (3.7, 0, 0))]
        originals = [a.xyz.copy() for a in lig + protein]
        result = _burial_from_atoms(lig, protein, n_points=960)
        for item, before in zip(lig + protein, originals):
            np.testing.assert_array_equal(item.xyz, before)
        self.assertFalse(result["coordinates_changed"])

    def test_invalid_inputs_fail_instead_of_emitting_misleading_fraction(self):
        for n in (0, -1, 3.2, True):
            with self.assertRaisesRegex(ValueError, "positive integer"):
                _burial_from_atoms(ligand(), [], n_points=n)
        for lig in ([], ligand()[1:], [replace(ligand()[0], element="O")]):
            with self.assertRaisesRegex(ValueError, "one carbon"):
                _burial_from_atoms(lig, [], n_points=960)
        with self.assertRaisesRegex(ValueError, "unique"):
            _burial_from_atoms(ligand() + [ligand()[0]], [], n_points=960)
        with self.assertRaisesRegex(ValueError, "Nonfinite"):
            _burial_from_atoms(ligand(), [atom("CB", (math.nan, 0, 0))], n_points=960)
        with self.assertRaisesRegex(ValueError, "Unsupported"):
            _burial_from_atoms(ligand(), [atom("ZN", (100, 0, 0), element="ZN")], n_points=960)
        # A coincident full-size ligand atom hides every target sphere point.
        zero_free = [ligand()[0], replace(ligand()[1], xyz=np.zeros(3), element="S")]
        with self.assertRaisesRegex(ValueError, "zero or invalid"):
            _burial_from_atoms(zero_free, [], n_points=960)

    def test_public_reader_selects_ligand_and_protein_scope(self):
        items = ligand() + [atom("CA", (3.5, 0, 0)),
                            atom("O", (0, 0, 0), chain="W", record="HETATM", residue="HOH", element="O"),
                            atom("CB", (0, 0, 0), chain="A"),
                            atom("H1", (0, 0, 0), element="H")]
        with tempfile.TemporaryDirectory() as folder:
            path = Path(folder) / "test.pdb"
            path.write_text("".join(pdb_line(i, a) for i, a in enumerate(items, 1)))
            result = analyze_c18_burial(path, n_points=960)
            self.assertEqual(result["ligand_heavy_atom_count"], 20)
            self.assertEqual(result["protein_heavy_atom_count"], 1)
            self.assertGreater(result["fraction_buried"], 0)
            path.write_text("".join(pdb_line(i, a) for i, a in enumerate(items[1:], 1)))
            with self.assertRaisesRegex(ValueError, "twenty"):
                analyze_c18_burial(path, n_points=960)

    @unittest.skipUnless(all(p.exists() for p in REAL_CASES), "Local audited m7/candidate structures absent")
    def test_audited_real_cases_equal_full_computation(self):
        expected = ((40.128479262159715, 40.191376565078464),
                    (35.47407884617254, 35.34828424033505))
        for case_index, path in enumerate(REAL_CASES):
            atoms = read_steric_atoms(path)
            lig = [a for a in atoms if a.record == "HETATM" and a.chain == "X" and a.residue == "E2L"]
            protein = [a for a in atoms if a.record == "ATOM" and a.chain in {"H", "L"}]
            for point_index, n in enumerate((960, 1920)):
                result = self.assert_full_equivalence(lig, protein, n)
                self.assertAlmostEqual(result["isolated_full_ligand_SASA_A2"], expected[case_index][point_index])
                self.assertEqual(result["fraction_buried"], 1)
                self.assertLess(result["measured_complex_heavy_atom_count"], 80)


if __name__ == "__main__":
    unittest.main()
