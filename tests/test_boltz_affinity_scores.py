"""Native Boltz affinity exports join only explicitly identified model-0 jobs."""

from __future__ import annotations

import json
from pathlib import Path
from tempfile import TemporaryDirectory
import unittest

from fv_energy_analysis.model_scores import (
    BOLTZ2_AFFINITY_FIELDS, BOLTZ2_FIELDS, GNINA_FIELDS, add_model_scores, score_columns,
)
from fv_energy_analysis.report import write_csv, write_json


JOB = "m7_NISE_CDRL3_AQYDQRPHT_estradiol"


def _native(job: str = JOB, **changes) -> dict[str, str]:
    row = {"condition": "with_msa", "job_name": job, "result_status": "success"}
    row.update({field: str((index + 1) / 10) for index, field in enumerate(BOLTZ2_AFFINITY_FIELDS)})
    row.update(changes)
    return row


def _wrapped(job: str = JOB, model: int = 0) -> str:
    return f"{job}__{job}_model_{model}_relaxed.pdb"


class NativeBoltzAffinityTests(unittest.TestCase):
    def setUp(self):
        self.directory = TemporaryDirectory()
        self.addCleanup(self.directory.cleanup)
        self.root = Path(self.directory.name)
        self.source = self.root / "boltz_affinity.csv"
        self.payload_path = self.root / "batch_summary_grouped_fingerprints.json"
        self.original = {
            "columns": [
                {"key": "candidate", "label": "Complex", "section": "Identity", "kind": "text"},
                {"key": "rosetta_interface_REU", "label": "Interface", "section": "RosettaLigand", "kind": "number"},
            ],
            "rows": [{"candidate": _wrapped(), "rosetta_interface_REU": -23.75, "network_delta": -4.25}],
            "parent": {"candidate": "Parent reference", "rosetta_interface_REU": -16.947652},
        }
        write_json(self.payload_path, self.original)
        write_json(self.root / "batch_manifest.json", {"parent": "e2pico_scfv_m7_estradiol_model_0_relaxed.pdb"})
        write_csv(self.source, [_native()])

    def payload(self):
        return json.loads(self.payload_path.read_text(encoding="utf-8"))

    def test_137_jobs_join_six_affinities_preserving_analysis_and_null_absent_fields(self):
        alphabet = "ACDEFGHIKLMNPQRSTVWY"
        jobs = [f"m7_NISE_CDRL3_AQYDQRP{alphabet[i // 20]}{alphabet[i % 20]}_estradiol" for i in range(137)]
        original_rows = [{"candidate": _wrapped(job), "rosetta_interface_REU": -20 - i / 100,
                          "network_delta": -i / 50} for i, job in enumerate(jobs)]
        self.original["rows"] = original_rows
        write_json(self.payload_path, self.original)
        write_csv(self.source, [_native(job) for job in jobs] + [{key: "" for key in _native()}])

        add_model_scores(self.root, self.source)

        payload = self.payload()
        self.assertEqual(payload["model_scores"]["matched_candidates"], 137)
        self.assertEqual(payload["model_scores"]["available_fields"], list(BOLTZ2_AFFINITY_FIELDS))
        self.assertEqual(payload["model_scores"]["available_sections"], ["Boltz-2"])
        self.assertEqual(payload["columns"][:2], self.original["columns"])
        for row, original in zip(payload["rows"], original_rows):
            self.assertEqual({key: row[key] for key in original}, original)
            for index, field in enumerate(BOLTZ2_AFFINITY_FIELDS):
                self.assertEqual(row[f"boltz2_{field}"], (index + 1) / 10)
            for field in BOLTZ2_FIELDS[:9]:
                self.assertIsNone(row[f"boltz2_{field}"])
            for field in GNINA_FIELDS:
                self.assertIsNone(row[f"gnina_{field}"])
        self.assertFalse(payload["model_scores"]["parent_scores_available"])
        self.assertTrue(all(payload["parent"][column["key"]] is None for column in score_columns()))
        provenance = payload["model_scores"]["mapping"][original_rows[0]["candidate"]]
        self.assertEqual(provenance, {"source_model": jobs[0], "source_row": 2})

    def test_optional_native_boltz_columns_are_retained_if_present(self):
        write_csv(self.source, [_native(confidence_score="0.91", complex_plddt="0.88")])
        add_model_scores(self.root, self.source)
        payload = self.payload()
        self.assertEqual(payload["rows"][0]["boltz2_confidence_score"], 0.91)
        self.assertEqual(payload["rows"][0]["boltz2_complex_plddt"], 0.88)
        self.assertIn("confidence_score", payload["model_scores"]["available_fields"])
        self.assertNotIn("ptm", payload["model_scores"]["available_fields"])

    def test_exact_native_job_identity_is_accepted_without_aliasing(self):
        self.original["rows"][0]["candidate"] = f"{JOB}.pdb"
        write_json(self.payload_path, self.original)
        add_model_scores(self.root, self.source)
        self.assertEqual(self.payload()["model_scores"]["matched_candidates"], 1)

    def test_duplicate_native_job_is_rejected_even_when_values_agree(self):
        write_csv(self.source, [_native(), _native()])
        with self.assertRaisesRegex(ValueError, "Duplicate score mapping"):
            add_model_scores(self.root, self.source)

    def test_only_full_identical_wrapper_for_model_zero_is_accepted(self):
        bad_names = (
            _wrapped(model=1),
            f"{JOB}__other_job_model_0_relaxed.pdb",
            f"other_job__{JOB}_model_0_relaxed.pdb",
            f"prefix_{JOB}__{JOB}_model_0_relaxed.pdb",
            f"{JOB}__{JOB}_model_0_extra_relaxed.pdb",
            f"{JOB}_model_0_relaxed.pdb",
        )
        for name in bad_names:
            with self.subTest(name=name):
                self.original["rows"][0]["candidate"] = name
                write_json(self.payload_path, self.original)
                before = self.payload_path.read_bytes()
                with self.assertRaisesRegex(ValueError, "No exact score mapping"):
                    add_model_scores(self.root, self.source)
                self.assertEqual(self.payload_path.read_bytes(), before)

    def test_non_success_native_status_is_rejected_before_writes(self):
        for status in ("error", "failed", ""):
            with self.subTest(status=status):
                write_csv(self.source, [_native(result_status=status)])
                before = self.payload_path.read_bytes()
                with self.assertRaisesRegex(ValueError, "Unsuccessful result_status"):
                    add_model_scores(self.root, self.source)
                self.assertEqual(self.payload_path.read_bytes(), before)

    def test_nonempty_row_without_job_identity_is_not_skipped(self):
        write_csv(self.source, [_native(job="")])
        with self.assertRaisesRegex(ValueError, "Missing or invalid job_name"):
            add_model_scores(self.root, self.source)

    def test_all_six_native_affinity_columns_are_required(self):
        row = _native()
        del row["affinity_pred_value2"]
        write_csv(self.source, [row])
        with self.assertRaisesRegex(ValueError, "missing required columns: affinity_pred_value2"):
            add_model_scores(self.root, self.source)

    def test_native_scores_are_not_allowed_to_be_nonfinite_or_nonnumeric(self):
        for value in ("nan", "inf", "not a score"):
            with self.subTest(value=value):
                write_csv(self.source, [_native(affinity_pred_value=value)])
                before = self.payload_path.read_bytes()
                with self.assertRaises(ValueError):
                    add_model_scores(self.root, self.source)
                self.assertEqual(self.payload_path.read_bytes(), before)

    def test_parent_can_only_match_an_explicit_native_identity(self):
        write_json(self.root / "batch_manifest.json", {"parent": _wrapped()})
        add_model_scores(self.root, self.source)
        self.assertTrue(self.payload()["model_scores"]["parent_scores_available"])
        self.assertEqual(self.payload()["parent"]["boltz2_affinity_pred_value"], 0.1)

    def test_native_alias_cannot_be_applied_to_merged_score_sources(self):
        row = {field: "0.1" for field in BOLTZ2_FIELDS + GNINA_FIELDS}
        row["input"] = JOB + ".cif"
        write_csv(self.source, [row])
        with self.assertRaisesRegex(ValueError, "No exact score mapping"):
            add_model_scores(self.root, self.source)


if __name__ == "__main__":
    unittest.main()
