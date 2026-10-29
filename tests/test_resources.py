"""Runtime resources must work from both source checkouts and installed wheels."""

from pathlib import Path
import unittest

from fv_energy_analysis.resources import molfile_to_params_path, resource_path


class ResourceTests(unittest.TestCase):
    def test_required_files_are_available(self):
        for relative in (
            "scripts/prepare_ligand_complexes.py",
            "scripts/pdb_to_arpeggio_cif.py",
            "scripts/discovery_studio_batch_interactions.pl",
            "docs/batch_output_columns.md",
            "docs/discovery_studio_columns.md",
            "docs/interaction_fingerprint_columns.md",
            "docs/e2_architecture.md",
        ):
            with self.subTest(resource=relative):
                self.assertTrue(resource_path(relative).is_file())

    def test_missing_file_is_explicit(self):
        with self.assertRaisesRegex(FileNotFoundError, "was not found"):
            resource_path("missing-resource-for-test.txt")

    def test_explicit_helper_path_is_respected(self):
        # Resolution only checks existence; this fixture is never executed.
        fixture = Path(__file__).resolve()
        self.assertEqual(molfile_to_params_path(fixture), fixture)


if __name__ == "__main__":
    unittest.main()
