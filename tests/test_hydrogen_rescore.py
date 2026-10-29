"""Pure safety and reporting checks for a separate hydrogen-only batch rescore.

These tests do not import PyRosetta, read real complexes, or run interaction
detectors.  All inputs and immutable-output artifacts are temporary fixtures.
"""
from copy import deepcopy
import csv
from hashlib import sha256
import json
from pathlib import Path
from tempfile import TemporaryDirectory
import unittest

from fv_energy_analysis.hydrogen_rescore import (
    make_summary_row,
    read_source,
    reuse_completed,
)


def _save_json(path: Path, value: dict) -> None:
    path.write_text(json.dumps(value), encoding="utf-8")


class HydrogenRescoreSourceTests(unittest.TestCase):
    def setUp(self):
        self.temporary = TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.source = Path(self.temporary.name) / "source"
        self.source.mkdir()
        self.manifest = {
            "status": "complete", "failure_count": 0, "completed_count": 2,
            "protocol": {
                "weights": "ligand", "legacy_pre_talaris": True,
                "ligand_chain": "X", "params_sha256": "a" * 64,
            },
        }
        self.rows = [
            {"candidate": "candidate_A.pdb", "substitutions": "E99F"},
            {"candidate": "candidate_B.pdb", "substitutions": "R100T;I102N"},
        ]
        self.write_source()

    def write_source(self):
        _save_json(self.source / "rescore_manifest.json", self.manifest)
        with (self.source / "rosettaligand_rescore.csv").open(
                "w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(handle, fieldnames=["candidate", "substitutions"])
            writer.writeheader()
            writer.writerows(self.rows)

    def test_reads_completed_batch_without_loading_case_pdbs(self):
        # Deliberately no result.json, parent PDB, or candidate PDB is present.
        manifest, rows = read_source(self.source, expected_count=2)
        self.assertEqual(manifest, self.manifest)
        self.assertEqual(rows, self.rows)

    def test_optional_expected_count_is_not_required(self):
        _, rows = read_source(self.source)
        self.assertEqual(len(rows), 2)

    def test_refuses_unfinished_or_partially_failed_source(self):
        for status, failures in (("running", 0), ("failed", 1),
                                 ("complete_with_failures", 1), ("complete", 1)):
            with self.subTest(status=status, failures=failures):
                self.manifest.update(status=status, failure_count=failures)
                self.write_source()
                with self.assertRaises(ValueError):
                    read_source(self.source)

    def test_refuses_manifest_and_csv_count_disagreement(self):
        self.manifest["completed_count"] = 3
        self.write_source()
        with self.assertRaises(ValueError):
            read_source(self.source)

    def test_refuses_unexpected_dataset_size(self):
        with self.assertRaises(ValueError):
            read_source(self.source, expected_count=1397)

    def test_refuses_duplicate_candidate_identifiers(self):
        self.rows[1]["candidate"] = self.rows[0]["candidate"]
        self.write_source()
        with self.assertRaises(ValueError):
            read_source(self.source)

    def test_requires_historical_rosettaligand_weights_and_legacy_mode(self):
        for weights, legacy in (("ref2015", False), ("ligand", False),
                                ("ref2015", True)):
            with self.subTest(weights=weights, legacy=legacy):
                self.manifest["protocol"].update(
                    weights=weights, legacy_pre_talaris=legacy)
                self.write_source()
                with self.assertRaises(ValueError):
                    read_source(self.source)


class HydrogenRescoreSummaryTests(unittest.TestCase):
    def setUp(self):
        # Both parent and mutant change on hydrogen optimization.  Comparing
        # mutant-after to parent-before would give the wrong relative effect.
        self.parent_before = {
            "bound_total": -100.0, "interface_delta_X": -16.95,
            "if_X_fa_atr": -20.0, "if_X_fa_rep": 1.21, "if_X_fa_elec": -0.5,
            "if_X_fa_sol": 5.0, "if_X_fa_pair": 0.0,
            "if_X_hbond_bb_sc": -0.5, "if_X_hbond_sc": -2.16,
        }
        self.parent_after = {
            "bound_total": -102.0, "interface_delta_X": -17.0,
            "if_X_fa_atr": -20.0, "if_X_fa_rep": 1.10, "if_X_fa_elec": -0.6,
            "if_X_fa_sol": 5.0, "if_X_fa_pair": 0.0,
            "if_X_hbond_bb_sc": -0.5, "if_X_hbond_sc": -2.0,
        }
        self.before = {
            "bound_total": -90.0, "interface_delta_X": -3.96,
            "if_X_fa_atr": -20.33, "if_X_fa_rep": 13.92, "if_X_fa_elec": -0.4,
            "if_X_fa_sol": 5.54, "if_X_fa_pair": 0.0,
            "if_X_hbond_bb_sc": -0.1, "if_X_hbond_sc": -2.59,
        }
        self.after = {
            "bound_total": -110.0, "interface_delta_X": -16.34,
            "if_X_fa_atr": -20.33, "if_X_fa_rep": 2.56, "if_X_fa_elec": -0.9,
            "if_X_fa_sol": 5.54, "if_X_fa_pair": 0.0,
            "if_X_hbond_bb_sc": -0.1, "if_X_hbond_sc": -3.11,
        }
        self.optimization = {
            "selected_trial": 1, "optimization_seconds": 0.07,
            "max_heavy_displacement_A": 0.0, "moved_hydrogen_count": 42,
        }

    def row(self):
        return make_summary_row(
            "candidate.pdb", "E99F;R100T;I102N", self.before, self.after,
            self.parent_before, self.parent_after, Path("optimized/candidate.pdb"),
            self.optimization,
        )

    def test_new_interface_delta_uses_optimized_parent_and_mutant(self):
        row = self.row()
        self.assertAlmostEqual(row["interface_parent_REU"], -17.0)
        self.assertAlmostEqual(row["interface_mutant_REU"], -16.34)
        self.assertAlmostEqual(row["interface_delta_REU"], 0.66)

    def test_before_values_retain_matched_original_comparison(self):
        row = self.row()
        self.assertAlmostEqual(row["before_interface_parent_REU"], -16.95)
        self.assertAlmostEqual(row["before_interface_mutant_REU"], -3.96)
        self.assertAlmostEqual(row["before_interface_delta_REU"], 12.99)
        self.assertAlmostEqual(row["hydrogen_optimization_interface_change_REU"], -12.38)

    def test_weighted_components_are_recalculated_from_new_metrics(self):
        row = self.row()
        self.assertAlmostEqual(row["parent_fa_rep_REU"], 1.10)
        self.assertAlmostEqual(row["mutant_fa_rep_REU"], 2.56)
        self.assertAlmostEqual(row["delta_fa_rep_REU"], 1.46)
        self.assertAlmostEqual(row["parent_fa_elec_REU"], -0.6)
        self.assertAlmostEqual(row["mutant_fa_elec_REU"], -0.9)
        self.assertAlmostEqual(row["delta_fa_elec_REU"], -0.3)

    def test_summary_does_not_mutate_saved_or_new_metrics(self):
        original = deepcopy((self.before, self.after, self.parent_before,
                             self.parent_after, self.optimization))
        self.row()
        self.assertEqual((self.before, self.after, self.parent_before,
                          self.parent_after, self.optimization), original)

    def test_stale_geometric_fingerprints_are_not_fabricated(self):
        row = self.row()
        # No Arpeggio, PLIP, or DS detector ran for these updated H positions.
        stale_prefixes = ("arpeggio_", "plip_", "ds_", "discovery_studio_")
        self.assertFalse(any(key.lower().startswith(stale_prefixes) for key in row))
        self.assertEqual(row["candidate"], "candidate.pdb")


class HydrogenRescoreReuseTests(unittest.TestCase):
    def setUp(self):
        self.temporary = TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.output = Path(self.temporary.name)
        self.record_path = self.output / "result.json"
        self.identity = {
            "pdb_sha256": "a" * 64,
            "protocol": {"weights": "ligand", "repeats": 1, "seed": 20261022},
        }
        self.artifact = self.output / "optimized.pdb"
        self.artifact.write_text("REMARK immutable optimized structure\nEND\n", encoding="utf-8")
        self.saved = {
            "status": "complete", "identity": self.identity,
            "artifact_sha256": {
                str(self.artifact.resolve()): sha256(self.artifact.read_bytes()).hexdigest(),
            },
        }
        _save_json(self.record_path, self.saved)

    def test_absent_result_is_not_a_reusable_completion(self):
        self.assertIsNone(reuse_completed(self.output / "absent.json", self.identity))

    def test_unchanged_completion_and_artifact_are_reused(self):
        before_record = self.record_path.read_bytes()
        before_artifact = self.artifact.read_bytes()
        self.assertEqual(reuse_completed(self.record_path, self.identity), self.saved)
        self.assertEqual(self.record_path.read_bytes(), before_record)
        self.assertEqual(self.artifact.read_bytes(), before_artifact)

    def test_modified_optimized_structure_is_rejected(self):
        self.artifact.write_text("REMARK modified structure\nEND\n", encoding="utf-8")
        with self.assertRaises(ValueError):
            reuse_completed(self.record_path, self.identity)

    def test_missing_optimized_structure_is_rejected(self):
        self.artifact.unlink()
        with self.assertRaises(ValueError):
            reuse_completed(self.record_path, self.identity)

    def test_modified_source_structure_or_protocol_is_rejected(self):
        for identity in (
            {**self.identity, "pdb_sha256": "b" * 64},
            {**self.identity, "protocol": {"weights": "ref2015", "repeats": 1, "seed": 20261022}},
            {**self.identity, "protocol": {"weights": "ligand", "repeats": 3, "seed": 20261022}},
            {**self.identity, "protocol": {"weights": "ligand", "repeats": 1, "seed": 42}},
        ):
            with self.subTest(identity=identity):
                with self.assertRaises(ValueError):
                    reuse_completed(self.record_path, identity)

    def test_unfinished_record_cannot_masquerade_as_completion(self):
        self.saved["status"] = "running"
        _save_json(self.record_path, self.saved)
        with self.assertRaises(ValueError):
            reuse_completed(self.record_path, self.identity)


if __name__ == "__main__":
    unittest.main()
