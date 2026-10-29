import csv
import json
from pathlib import Path
import tempfile
import unittest

from fv_energy_analysis.rescore_report import prepare_report, report_columns, SCORE_FIELDS


class RescoreReportTests(unittest.TestCase):
    def test_columns_preserve_tools_and_add_active_legacy_terms(self):
        columns = [
            {"key": "rosetta_hbond_sc_REU", "label": "hbond_sc only (REU)", "section": "Rosetta", "kind": "energy"},
            {"key": "ds_pi_stack", "label": "Pi-stacks (unique residues)", "section": "Discovery Studio", "kind": "integer"},
        ]
        new = report_columns(columns)
        self.assertEqual([c["key"] for c in new[:3]], ["rosetta_fa_pair_REU", "rosetta_hbond_bb_sc_REU", "rosetta_hbond_sc_REU"])
        self.assertTrue(all(c["section"] == "RosettaLigand" for c in new[:3]))
        self.assertEqual(new[-1], columns[-1])
        self.assertEqual(columns[0]["section"], "Rosetta")

    def _fixture(self, root):
        source, output = root / "source", root / "rescore"
        source.mkdir()
        output.mkdir()
        row = {"candidate": "mutant.pdb", "rank": 7, "mutations": "E99G", "ds_pi_stack": 1, "hcdr3_plddt_mean": 77.0}
        for key in SCORE_FIELDS:
            row[key] = 999
        payload = {"title": "Source", "rows": [row], "parent": {"candidate": "Parent reference", "ds_pi_stack": 0},
                   "columns": [{"key": "candidate", "label": "Complex", "section": "Identity", "kind": "text"},
                               {"key": "rosetta_hbond_sc_REU", "label": "hbond_sc only (REU)", "section": "Rosetta", "kind": "energy"}]}
        (source / "batch_summary_grouped_fingerprints.json").write_text(json.dumps(payload))
        self._csv(source / "batch_summary.csv", [{"candidate": "mutant.pdb", "candidate_path": "inputs/mutant.pdb", "status": "complete", "error": ""},
                                                  {"candidate": "bad.pdb", "candidate_path": "inputs/bad.pdb", "status": "failed", "error": "Wrong stereochemistry"}])
        (source / "ligand_validation_diagnostics.csv").write_text("candidate,differences\nbad.pdb,stereochemistry\n")
        (output / "rescore_manifest.json").write_text(json.dumps({"status": "complete", "failure_count": 0, "excluded_source_failures": 1}))
        score = {"candidate": "mutant.pdb", "weights": "ligand", "substitutions": "E99G", "rank_by_interface_delta": 1,
                 "interface_delta_REU": -1}
        for state in ("parent", "mutant"):
            for field in SCORE_FIELDS.values():
                score[f"{state}_{field}"] = 0
            score[f"{state}_ligand_hbonds"] = 3
            score[f"{state}_fa_atr_REU"] = -10 if state == "parent" else -11
            score[f"{state}_interface_REU"] = score[f"{state}_fa_atr_REU"]
        self._csv(output / "rosettaligand_rescore.csv", [score])
        network = {"mutant.pdb": {"scorefunction": "ligand", "network_added_pair_type_events": 2, "network_removed_pair_type_events": 1,
                                  "network_changed_unique_pairs_energy_delta_REU": -0.5, "comparison_quality_warning": False}}
        (output / "rescore_network_summary.json").write_text(json.dumps(network))
        return source, output, score

    @staticmethod
    def _csv(path, records):
        with path.open("w", newline="") as stream:
            writer = csv.DictWriter(stream, fieldnames=list(records[0]))
            writer.writeheader()
            writer.writerows(records)

    def test_report_replaces_scores_preserves_geometry_and_lists_exclusions(self):
        with tempfile.TemporaryDirectory() as tmp:
            source, output, _ = self._fixture(Path(tmp))
            original = (source / "batch_summary_grouped_fingerprints.json").read_bytes()
            summary = prepare_report(source, output)
            self.assertEqual(summary["rows"], 1)
            payload = json.loads((output / "batch_summary_grouped_fingerprints.json").read_text())
            row = payload["rows"][0]
            self.assertEqual((row["rank"], row["rosetta_interface_REU"], row["rosetta_interface_delta_REU"]), (1, -11, -1))
            self.assertEqual((row["ds_pi_stack"], row["hcdr3_plddt_mean"]), (1, 77))
            self.assertEqual(payload["parent"]["rosetta_interface_REU"], -10)
            self.assertEqual(row["network_changed_unique_pairs_energy_delta_REU"], -0.5)
            self.assertFalse(row["comparison_quality_warning"])
            with (output / "batch_summary.csv").open() as stream:
                records = list(csv.DictReader(stream))
            self.assertEqual(records[1]["status"], "excluded_source_validation")
            self.assertEqual((source / "batch_summary_grouped_fingerprints.json").read_bytes(), original)

    def test_report_rejects_missing_candidate(self):
        with tempfile.TemporaryDirectory() as tmp:
            source, output, score = self._fixture(Path(tmp))
            score["candidate"] = "another.pdb"
            self._csv(output / "rosettaligand_rescore.csv", [score])
            with self.assertRaisesRegex(ValueError, "match the source workbook"):
                prepare_report(source, output)

    def test_report_rejects_incomplete_decomposition(self):
        with tempfile.TemporaryDirectory() as tmp:
            source, output, score = self._fixture(Path(tmp))
            score["mutant_fa_pair_REU"] = -1
            self._csv(output / "rosettaligand_rescore.csv", [score])
            with self.assertRaisesRegex(ValueError, "does not sum"):
                prepare_report(source, output)


if __name__ == "__main__":
    unittest.main()
