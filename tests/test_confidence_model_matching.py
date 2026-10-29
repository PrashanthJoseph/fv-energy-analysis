"""Strict original-model matching for duplicated relaxed-export filenames."""

from __future__ import annotations

import importlib.util
from pathlib import Path
from tempfile import TemporaryDirectory
import unittest

from fv_energy_analysis.confidence import _matching_cifs, extract_hcdr3_plddt


class ConfidenceModelMatchingTests(unittest.TestCase):
    def test_plain_model_name_is_exact(self):
        source = Path("original/job_model_0.cif")
        self.assertEqual(_matching_cifs("job_model_0_relaxed.pdb", {"job_model_0": [source]}), [source])

    def test_strict_repeated_job_wrapper_is_unwrapped(self):
        source = Path("original/m7_CDRL3_AQYDRRPYT_estradiol_model_0.cif")
        candidate = "m7_CDRL3_AQYDRRPYT_estradiol__m7_CDRL3_AQYDRRPYT_estradiol_model_0_relaxed.pdb"
        self.assertEqual(_matching_cifs(candidate, {source.stem.casefold(): [source]}), [source])

    def test_unrelated_prefix_is_not_discarded(self):
        self.assertEqual(_matching_cifs("other__job_model_0_relaxed.pdb", {"job_model_0": [Path("a.cif")]}), [])

    def test_model_number_is_never_substituted(self):
        self.assertEqual(_matching_cifs("job__job_model_1_relaxed.pdb", {"job_model_0": [Path("a.cif")]}), [])

    def test_exact_and_unwrapped_equivalents_remain_ambiguous(self):
        direct, unwrapped = Path("a/job__job_model_0.cif"), Path("b/job_model_0.cif")
        sources = {"job__job_model_0": [direct], "job_model_0": [unwrapped]}
        self.assertEqual(_matching_cifs("job__job_model_0_relaxed.pdb", sources), [direct, unwrapped])

    def test_duplicate_index_entries_do_not_create_false_ambiguity(self):
        source = Path("job_model_0.cif")
        self.assertEqual(_matching_cifs("job_model_0_relaxed.pdb", {"job_model_0": [source, source]}), [source])


@unittest.skipUnless(importlib.util.find_spec("gemmi"), "optional gemmi is not installed")
class ConfidenceCompleteFvMatchingTests(unittest.TestCase):
    def setUp(self):
        self.directory = TemporaryDirectory()
        self.addCleanup(self.directory.cleanup)
        self.root = Path(self.directory.name)
        self.pdb = self.root / "canonical.pdb"
        lines = []
        for serial, number in enumerate(range(95, 114), 1):
            lines.append(f"ATOM  {serial:5d}  CA  ALA H{number:4d}    {serial:8.3f}{0.:8.3f}{0.:8.3f}{1.:6.2f}{0.:6.2f}           C  \n")
        for serial, name in ((20, "GLY"), (21, "SER")):
            lines.append(f"ATOM  {serial:5d}  CA  {name:3s} L{serial - 19:4d}    {serial:8.3f}{0.:8.3f}{0.:8.3f}{1.:6.2f}{0.:6.2f}           C  \n")
        self.pdb.write_text("".join(lines) + "END\n", encoding="utf-8")

    def _cif(self, light_names=("GLY", "SER"), duplicate_light=False, omit_light=False):
        path = self.root / "original.cif"
        lines = ["data_test", "loop_", "_ma_qa_metric.id", "_ma_qa_metric.type", "_ma_qa_metric.mode",
                 "1 pLDDT local", "#", "loop_", "_atom_site.label_atom_id", "_atom_site.label_comp_id",
                 "_atom_site.label_asym_id", "_atom_site.label_seq_id", "_atom_site.pdbx_PDB_model_num"]
        # Labels deliberately differ from canonical H/L.
        chains = [("A", ("ALA",) * 19)]
        if not omit_light:
            chains.append(("B", light_names))
        if duplicate_light:
            chains.append(("C", light_names))
        for chain, names in chains:
            lines.extend(f"CA {name} {chain} {i} 1" for i, name in enumerate(names, 1))
        lines += ["#", "loop_", "_ma_qa_metric_local.model_id", "_ma_qa_metric_local.label_asym_id",
                  "_ma_qa_metric_local.label_seq_id", "_ma_qa_metric_local.metric_id", "_ma_qa_metric_local.metric_value"]
        lines.extend(f"1 A {i} 1 80.0" for i in range(1, 20))
        path.write_text("\n".join(lines) + "\n#\n", encoding="utf-8")
        return path

    def test_complete_light_sequence_match_allows_transfer(self):
        result = extract_hcdr3_plddt(self._cif(), self.pdb)
        self.assertEqual(result["source_cif_chain"], "A")
        self.assertEqual(result["values"]["hcdr3_plddt_mean"], 80.)

    def test_wrong_light_variant_rejected_despite_identical_heavy_chain(self):
        with self.assertRaisesRegex(ValueError, "canonical L.*found 0"):
            extract_hcdr3_plddt(self._cif(light_names=("GLY", "THR")), self.pdb)

    def test_missing_light_chain_rejected(self):
        with self.assertRaisesRegex(ValueError, "canonical L.*found 0"):
            extract_hcdr3_plddt(self._cif(omit_light=True), self.pdb)

    def test_ambiguous_light_chain_rejected(self):
        with self.assertRaisesRegex(ValueError, "canonical L.*found 2"):
            extract_hcdr3_plddt(self._cif(duplicate_light=True), self.pdb)


if __name__ == "__main__":
    unittest.main()
