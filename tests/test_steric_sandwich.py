from dataclasses import replace
import math
from pathlib import Path
import tempfile
import unittest

import numpy as np

from fv_energy_analysis.steric_sandwich import (
    StericAtom, analyze_steric_sandwich, read_steric_atoms, steric_wall_coverage,
)


RING = tuple(f"C{i}" for i in range(1, 7))


def atom(name, xyz, *, residue="ALA", chain="H", number=1, element="C", record="ATOM"):
    return StericAtom(name, residue, chain, number, "", element, np.array(xyz, dtype=float), record=record)


def ligand():
    ring = [atom(name, (math.cos(i * math.pi / 3), math.sin(i * math.pi / 3), 0),
                 residue="E2L", chain="X", record="HETATM") for i, name in enumerate(RING)]
    return ring + [atom("C10", (0, 0, 1), residue="E2L", chain="X", record="HETATM"),
                   atom("C16", (0, 0, 0), residue="E2L", chain="X", record="HETATM")]


def walls(gap=0.5):
    # Ring carbons are one A from the axis; set exact nearest carbon gap.
    height = math.sqrt((3.4 + gap) ** 2 - 1)
    return [atom("CB", (0, 0, height), chain="H"), atom("CB", (0, 0, -height), chain="L")]


def measure(protein, **kwargs):
    return steric_wall_coverage(ligand(), protein, ring_atom_names=RING,
                               methyl_atom="C10", methyl_parent_atom="C16", **kwargs)


def pdb_line(serial, item, occupancy=1, alt=" ", element=None):
    x, y, z = item.xyz
    return (f"{item.record:<6}{serial:5d} {item.name:>4}{alt}{item.residue:>3} {item.chain:1}"
            f"{item.number:4d}{item.insertion:1}   {x:8.3f}{y:8.3f}{z:8.3f}"
            f"{occupancy:6.2f}{0:6.2f}          {(item.element if element is None else element):>2}\n")


class StericSandwichTests(unittest.TestCase):
    def test_two_distinct_opposite_walls_pass(self):
        result = measure(walls())
        self.assertTrue(result["steric_sandwich"])
        self.assertEqual(result["coverage"]["alpha_residue_count"], 1)
        self.assertEqual(result["coverage"]["opposing_residue_pair_count"], 1)
        self.assertEqual(len(result["contacts"]), 2)
        self.assertNotIn("rays", result)
        self.assertIn("not a published numerical rule", result["provenance"]["definition"])

    def test_empty_and_one_face_fail(self):
        for protein in ([], walls()[:1], walls()[1:]):
            self.assertFalse(measure(protein)["steric_sandwich"])

    def test_same_wrapping_residue_is_not_distinct_walls(self):
        protein = [replace(a, chain="H", number=2) for a in walls()]
        result = measure(protein)
        self.assertEqual(result["coverage"]["alpha_residue_count"], 1)
        self.assertEqual(result["coverage"]["opposing_residue_pair_count"], 0)
        self.assertFalse(result["steric_sandwich"])

    def test_exact_cutoff_and_just_outside(self):
        self.assertTrue(measure(walls(1.0))["steric_sandwich"])
        self.assertFalse(measure(walls(1.0001))["steric_sandwich"])

    def test_overlap_not_rewarded(self):
        result = measure(walls(-0.05))
        self.assertFalse(result["steric_sandwich"])
        self.assertEqual(result["coverage"]["overlapping_atom_pair_count"], 12)
        self.assertGreater(result["largest_ring_sidechain_vdw_overlap_A"], 0)

    def test_backbone_water_hydrogen_and_other_records_excluded(self):
        for change in ({"name": "CA"}, {"name": "H1", "element": "H"},
                       {"name": "O", "residue": "HOH", "record": "HETATM", "element": "O"}):
            result = measure([replace(a, **change) for a in walls()])
            self.assertFalse(result["steric_sandwich"])
            self.assertEqual(result["protein_sidechain_heavy_atom_count"], 0)

    def test_peripheral_opposite_atoms_fail(self):
        result = measure([atom("CB", (5, 0, 1)), atom("CB", (5, 0, -1), chain="L")])
        self.assertFalse(result["steric_sandwich"])
        self.assertEqual(result["coverage"]["opposing_residue_pair_count"], 0)

    def test_rigid_transform_invariance(self):
        angle = 0.723
        axis = np.array([1, 2, 3], dtype=float) / math.sqrt(14)
        cross = np.array([[0, -axis[2], axis[1]], [axis[2], 0, -axis[0]], [-axis[1], axis[0], 0]])
        rotation = np.eye(3) * math.cos(angle) + (1 - math.cos(angle)) * np.outer(axis, axis) + math.sin(angle) * cross
        transform = lambda atoms: [replace(a, xyz=a.xyz @ rotation + [10, -3, 9]) for a in atoms]
        original = measure(walls())
        changed = steric_wall_coverage(transform(ligand()), transform(walls()), ring_atom_names=RING,
                                      methyl_atom="C10", methyl_parent_atom="C16")
        self.assertEqual(original["coverage"], changed["coverage"])
        for old, new in zip(original["contacts"], changed["contacts"]):
            self.assertEqual(old["wall_atom"], new["wall_atom"])
            self.assertAlmostEqual(old["surface_gap_A"], new["surface_gap_A"], places=12)

    def test_atoms_audited_but_residue_pairs_deduplicated(self):
        protein = walls() + [replace(walls()[0], name="CG", xyz=walls()[0].xyz + [0.1, 0, 0])]
        result = measure(protein)
        self.assertEqual(len(result["contacts"]), 3)
        self.assertEqual(len(result["opposing_pairs"]), 2)
        self.assertEqual(result["coverage"]["opposing_residue_pair_count"], 1)

    def test_offset_contacts_crossing_footprint_pass(self):
        # Offset pair crosses outside the carbon-centre hexagon (x<=1),
        # but within the physical carbon footprint (x<=2.7).
        result = measure([atom("CB", (2.0, 0, 3.5)), atom("CB", (2.0, 0, -3.5), chain="L")])
        self.assertTrue(result["steric_sandwich"])
        crossing = result["opposing_pairs"][0]["plane_crossing_xyz_A"]
        self.assertGreater(crossing[0], 1)

    def test_sideways_rim_contacts_are_excluded(self):
        # Both atoms have allowable carbon gaps but spheres straddle the plane.
        result = measure([atom("CB", (4.5, 0, 0.5)), atom("CB", (4.5, 0, -0.5), chain="L")])
        self.assertEqual(result["contacts"], [])
        self.assertFalse(result["steric_sandwich"])

    def test_opposite_faces_without_crossing_footprint_fail(self):
        # Carbon separation is close, but both contacts lie beyond the footprint.
        result = measure([atom("CB", (3.8, 0, 2)), atom("CB", (3.8, 0, -2), chain="L")])
        self.assertEqual(len(result["contacts"]), 2)
        self.assertFalse(result["steric_sandwich"])

    def test_overlapping_nearest_atom_does_not_hide_valid_atom(self):
        protein = walls() + [atom("CG", (0, 0, 2.0))]
        result = measure(protein)
        self.assertTrue(result["steric_sandwich"])
        self.assertGreater(result["coverage"]["overlapping_atom_pair_count"], 0)
        self.assertTrue(all(c["wall_atom"] != "H:1:ALA:CG" for c in result["contacts"]))

    def test_numerical_boundaries_and_bad_coordinates(self):
        self.assertTrue(measure(walls(0.0))["steric_sandwich"])
        for cutoff in (-1, float("nan"), float("inf")):
            with self.assertRaises(ValueError):
                measure(walls(), maximum_surface_gap_A=cutoff)
        with self.assertRaises(ValueError):
            measure([replace(walls()[0], xyz=np.array([0., float("nan"), 3.]))])

    def test_bad_mapping_and_unsupported_radius_fail(self):
        for names in (RING[:5], ("C1",) * 6, (*RING[:5], "C99")):
            with self.assertRaises(ValueError):
                steric_wall_coverage(ligand(), walls(), ring_atom_names=names)
        with self.assertRaisesRegex(ValueError, "Bondi radius"):
            measure([replace(walls()[0], element="SE")])

    def test_pdb_scope_altloc_and_first_model(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "test.pdb"
            protein = walls()
            records = ligand() + protein + [replace(protein[0], chain="Z")]
            text = "MODEL        1\n" + "".join(pdb_line(i, a) for i, a in enumerate(records, 1))
            text += pdb_line(20, replace(protein[0], xyz=np.array([0, 0, 1.])), occupancy=0.1, alt="B")
            text += "ENDMDL\n" + pdb_line(21, atom("CB", [0, 0, 0], number=20)) + "ENDMDL\n"
            path.write_text(text)
            result = analyze_steric_sandwich(path, ligand_chain="X", ligand_resname="E2L", ligand_resseq=1,
                                            ring_atom_names=RING, protein_chains=("H", "L"),
                                            methyl_atom="C10", methyl_parent_atom="C16")
            self.assertTrue(result["steric_sandwich"])
            self.assertEqual(result["protein_sidechain_heavy_atom_count"], 2)
            self.assertEqual(len(read_steric_atoms(path)), len(records))


if __name__ == "__main__":
    unittest.main()
