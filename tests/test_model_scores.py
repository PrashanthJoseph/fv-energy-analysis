"""Saved-score joins must not change existing analysis or guess model matches."""

from __future__ import annotations

import json
from pathlib import Path
from tempfile import TemporaryDirectory
import unittest
from unittest.mock import patch

from fv_energy_analysis.model_scores import BOLTZ2_FIELDS, GNINA_FIELDS, _model_key, add_model_scores, score_columns
from fv_energy_analysis.report import write_csv, write_json


def _scores(identity: str = "mutant_model_0", **changes) -> dict[str, object]:
    row = {field: str((index + 1) / 100) for index, field in enumerate(BOLTZ2_FIELDS + GNINA_FIELDS)}
    row.update({"cif_file": f"nested/{identity}.cif", "input": f"nested/{identity}.cif",
                "chain_ptm_0": "0.9", "pair_iptm_0_1": "0.8", "design_type": "2"})
    row.update(changes)
    return row


class ModelScoreSchemaTests(unittest.TestCase):
    def test_twenty_five_individual_numeric_columns_in_two_sections(self):
        columns = score_columns()
        self.assertEqual(len(columns), 25)
        self.assertEqual(sum(c["section"] == "Boltz-2" for c in columns), 15)
        self.assertEqual(sum(c["section"] == "GNINA" for c in columns), 10)
        self.assertEqual(len({c["key"] for c in columns}), 25)
        self.assertTrue(all(c["kind"] == "number" for c in columns))
        self.assertEqual([c["label"] for c in columns if c["source_column"] == "complex_plddt"],
                         ["Complex pLDDT (0-1)"])

    def test_exact_normalization_accepts_both_path_styles_but_preserves_model_number(self):
        self.assertEqual(_model_key(r"D:\folder\mutant_model_0_relaxed.pdb"), "mutant_model_0")
        self.assertEqual(_model_key("/folder/mutant_model_0.mmcif"), "mutant_model_0")
        self.assertNotEqual(_model_key("mutant_model_0.cif"), _model_key("mutant_model_1.cif"))
        self.assertEqual(_model_key("a_relaxed_middle_model_0.cif"), "a_relaxed_middle_model_0")


class ModelScoreJoinTests(unittest.TestCase):
    def setUp(self):
        self.directory = TemporaryDirectory()
        self.addCleanup(self.directory.cleanup)
        self.root = Path(self.directory.name)
        self.source = self.root / "scores.csv"
        self.payload_path = self.root / "batch_summary_grouped_fingerprints.json"
        self.original = {
            "columns": [{"key": "candidate", "label": "Complex", "section": "Identity", "kind": "text"},
                        {"key": "rank", "label": "Rank", "section": "Identity", "kind": "integer"},
                        {"key": "hcdr3_plddt_97", "label": "H97 pLDDT", "section": "HCDR3 pLDDT", "kind": "number"}],
            "rows": [{"candidate": "mutant_model_0_relaxed.pdb", "rank": 17, "hcdr3_plddt_97": 91.25}],
            "parent": {"candidate": "Parent reference", "rank": None, "hcdr3_plddt_97": 85.0},
        }
        write_json(self.payload_path, self.original)
        write_json(self.root / "batch_manifest.json", {"parent": "parent_model_0_relaxed.pdb", "unrelated": True})
        write_csv(self.root / "batch_interaction_fingerprint_numeric.csv",
                  [{"candidate": "mutant_model_0_relaxed.pdb", "original_energy": -4.25}])
        write_csv(self.source, [_scores()])

    def payload(self):
        return json.loads(self.payload_path.read_text(encoding="utf-8"))

    def test_join_preserves_all_existing_values_columns_rank_and_unscaled_scores(self):
        output = add_model_scores(self.root, self.source)
        result = self.payload()
        self.assertTrue(output.is_file())
        self.assertEqual(result["columns"][:3], self.original["columns"])
        row = result["rows"][0]
        self.assertEqual({k: row[k] for k in self.original["rows"][0]}, self.original["rows"][0])
        self.assertEqual(row["boltz2_complex_plddt"], 0.06)
        self.assertEqual(row["gnina_rank1_CNNscore"], 0.17)
        self.assertEqual(len(result["columns"]), 28)
        manifest = json.loads((self.root / "batch_manifest.json").read_text(encoding="utf-8"))
        self.assertTrue(manifest["unrelated"])
        self.assertEqual(manifest["model_scores"]["source_csv"], str(self.source.resolve()))

    def test_unavailable_parent_scores_are_none_not_zero(self):
        add_model_scores(self.root, self.source)
        payload = self.payload()
        self.assertFalse(payload["model_scores"]["parent_scores_available"])
        self.assertTrue(all(payload["parent"][c["key"]] is None for c in score_columns()))

    def test_known_parent_can_be_matched_by_exact_model_name(self):
        write_csv(self.source, [_scores(), _scores("parent_model_0", confidence_score="0.8")])
        add_model_scores(self.root, self.source)
        self.assertTrue(self.payload()["model_scores"]["parent_scores_available"])
        self.assertEqual(self.payload()["parent"]["boltz2_confidence_score"], 0.8)

    def test_blank_score_remains_none_and_positive_affinity_is_not_clipped(self):
        write_csv(self.source, [_scores(ptm="", rank1_minimizedAffinity="110.2311")])
        add_model_scores(self.root, self.source)
        row = self.payload()["rows"][0]
        self.assertIsNone(row["boltz2_ptm"])
        self.assertEqual(row["gnina_rank1_minimizedAffinity"], 110.2311)

    def test_repeated_run_replaces_instead_of_duplicating_sections(self):
        add_model_scores(self.root, self.source)
        first = self.payload()
        add_model_scores(self.root, self.source)
        self.assertEqual(self.payload(), first)

    def test_model_one_is_never_used_for_model_zero(self):
        write_csv(self.source, [_scores("mutant_model_1")])
        original_bytes = self.payload_path.read_bytes()
        with self.assertRaisesRegex(ValueError, "No exact score mapping"):
            add_model_scores(self.root, self.source)
        self.assertEqual(self.payload_path.read_bytes(), original_bytes)

    def test_duplicate_mapping_is_rejected(self):
        write_csv(self.source, [_scores(), _scores()])
        with self.assertRaisesRegex(ValueError, "Duplicate score mapping"):
            add_model_scores(self.root, self.source)

    def test_missing_required_column_is_rejected_even_if_all_other_values_exist(self):
        row = _scores()
        del row["ptm"]
        write_csv(self.source, [row])
        with self.assertRaisesRegex(ValueError, "missing required columns: ptm"):
            add_model_scores(self.root, self.source)

    def test_all_matches_validate_before_any_report_write(self):
        payload = self.payload()
        payload["rows"].append({"candidate": "missing_model_0_relaxed.pdb", "rank": 1})
        write_json(self.payload_path, payload)
        original = {p.name: p.read_bytes() for p in self.root.iterdir() if p.is_file()}
        with self.assertRaisesRegex(ValueError, "No exact score mapping"):
            add_model_scores(self.root, self.source)
        self.assertEqual({p.name: p.read_bytes() for p in self.root.iterdir() if p.is_file()}, original)

    def test_nonfinite_and_nonnumeric_scores_are_rejected(self):
        original_bytes = self.payload_path.read_bytes()
        for value in ("nan", "inf", "-inf", "not a number"):
            with self.subTest(value=value):
                write_csv(self.source, [_scores(ptm=value)])
                with self.assertRaises(ValueError):
                    add_model_scores(self.root, self.source)
                self.assertEqual(self.payload_path.read_bytes(), original_bytes)

    def test_conflicting_source_identities_are_rejected(self):
        write_csv(self.source, [_scores(input="other_model_0.cif")])
        with self.assertRaisesRegex(ValueError, "conflicting model identity"):
            add_model_scores(self.root, self.source)

    def test_input_is_used_when_cif_file_is_blank(self):
        write_csv(self.source, [_scores(cif_file="", input=r"directory\mutant_model_0.cif")])
        add_model_scores(self.root, self.source)
        self.assertEqual(self.payload()["model_scores"]["matched_candidates"], 1)

    def test_detail_preserves_chain_pair_metrics_raw_provenance_and_unrelated_reports(self):
        result_dir = self.root / "mutant_result"
        result_dir.mkdir()
        write_json(result_dir / "automatic_summary.json", {"existing_network": {"formed": 3}})
        (result_dir / "automatic_summary.md").write_text("# Existing analysis\n", encoding="utf-8")
        write_csv(self.root / "batch_summary.csv", [{"candidate": "mutant_model_0_relaxed.pdb", "result_folder": str(result_dir)}])
        add_model_scores(self.root, self.source)
        add_model_scores(self.root, self.source)
        detail = json.loads((result_dir / "model_scores.json").read_text(encoding="utf-8"))
        self.assertEqual(detail["numeric_source_fields"]["chain_ptm_0"], 0.9)
        self.assertEqual(detail["numeric_source_fields"]["pair_iptm_0_1"], 0.8)
        self.assertEqual(detail["raw_source_row"]["chain_ptm_0"], "0.9")
        self.assertEqual(detail["source_row"], 2)
        summary = json.loads((result_dir / "automatic_summary.json").read_text(encoding="utf-8"))
        self.assertEqual(summary["existing_network"], {"formed": 3})
        markdown = (result_dir / "automatic_summary.md").read_text(encoding="utf-8")
        self.assertTrue(markdown.startswith("# Existing analysis"))
        self.assertEqual(markdown.count("<!-- model-scores:begin -->"), 1)

    def test_locked_csv_uses_clear_fallback_name(self):
        from fv_energy_analysis import model_scores

        original_write = model_scores.write_csv

        def locked(path, rows):
            if Path(path).name == "batch_summary_grouped_fingerprints.csv":
                raise PermissionError("open in spreadsheet application")
            return original_write(path, rows)

        with patch.object(model_scores, "write_csv", side_effect=locked):
            output = add_model_scores(self.root, self.source)
        self.assertEqual(output.name, "batch_summary_grouped_fingerprints_with_model_scores.csv")
        self.assertTrue(output.is_file())


if __name__ == "__main__":
    unittest.main()
