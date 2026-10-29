"""Original ModelCIF confidence transfer; no scoring or user inputs are run."""

from __future__ import annotations

import importlib.util
from pathlib import Path
from tempfile import TemporaryDirectory
import unittest
from unittest.mock import patch

from fv_energy_analysis.confidence import RESIDUES, confidence_columns, extract_hcdr3_plddt


NAMES = (
    "ALA", "GLY", "SER", "VAL", "TYR", "PHE", "ARG", "ASP", "ASN", "GLN",
    "TRP", "ILE", "LEU", "MET", "LYS", "HIS", "CYS", "THR", "GLU",
)


def _pdb(path: Path, *, chain: str = "H", numbers=None, names=NAMES, insertion: str = "") -> Path:
    numbers = range(95, 95 + len(names)) if numbers is None else numbers
    lines = []
    for serial, (number, name) in enumerate(zip(numbers, names), 1):
        code = insertion if number == 97 else ""
        # Deliberately zero relaxed B factors: they must never supply confidence.
        lines.append(
            f"ATOM  {serial:5d}  CA  {name:3s} {chain}{number:4d}{code:1s}   "
            f"{float(serial):8.3f}{0.0:8.3f}{0.0:8.3f}{1.0:6.2f}{0.0:6.2f}           C  \n"
        )
    path.write_text("".join(lines) + "END\n", encoding="utf-8")
    return path


def _cif(path: Path, *, declared: bool = True, ambiguous: bool = False,
         omitted_score: int | None = None, changed_score: tuple[int, str] | None = None,
         duplicate_score: bool = False, source_names=NAMES) -> Path:
    lines = ["data_confidence", "#", "loop_", "_ma_qa_metric.id", "_ma_qa_metric.type", "_ma_qa_metric.mode",
             "1 pLDDT local" if declared else "1 other local", "#", "loop_",
             "_atom_site.label_atom_id", "_atom_site.label_comp_id", "_atom_site.label_asym_id",
             "_atom_site.label_seq_id", "_atom_site.pdbx_PDB_model_num", "_atom_site.B_iso_or_equiv"]
    # CIF H is light and CIF L is heavy, opposite to the normalized analysis PDB.
    chains = [("H", ("PRO",) * len(NAMES)), ("L", source_names)]
    if ambiguous:
        chains.append(("X", source_names))
    for chain, names in chains:
        lines.extend(f"CA {name} {chain} {i} 1 2.0" for i, name in enumerate(names, 1))
    lines += ["#", "loop_", "_ma_qa_metric_local.model_id", "_ma_qa_metric_local.label_asym_id",
              "_ma_qa_metric_local.label_seq_id", "_ma_qa_metric_local.metric_id", "_ma_qa_metric_local.metric_value"]
    for chain, names in chains:
        for i, _ in enumerate(names, 1):
            if i == omitted_score and chain == "L":
                continue
            value = str(60.0 + i) if chain == "L" else "10.0"
            if changed_score and i == changed_score[0] and chain == "L":
                value = changed_score[1]
            lines.append(f"1 {chain} {i} 1 {value}")
    if duplicate_score:
        lines.append("1 L 3 1 63.0")
    lines.append("#")
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path


class ConfidenceSchemaTests(unittest.TestCase):
    def test_separate_section_has_mean_min_and_seventeen_residues(self):
        columns = confidence_columns()
        self.assertEqual(len(columns), 19)
        self.assertEqual({column["section"] for column in columns}, {"HCDR3 pLDDT"})
        self.assertEqual([column["key"] for column in columns[2:]],
                         [f"hcdr3_plddt_{number}" for number in range(97, 114)])

    def test_missing_optional_dependency_is_actionable(self):
        with patch.dict("sys.modules", {"gemmi": None}):
            with self.assertRaisesRegex(RuntimeError, r"fv-energy-analysis\[reporting\]"):
                extract_hcdr3_plddt(Path("unused.cif"), Path("unused.pdb"))


@unittest.skipUnless(importlib.util.find_spec("gemmi"), "optional gemmi is not installed")
class ConfidenceTransferTests(unittest.TestCase):
    def setUp(self):
        self.directory = TemporaryDirectory()
        self.addCleanup(self.directory.cleanup)
        self.root = Path(self.directory.name)
        self.pdb = _pdb(self.root / "normalized.pdb")
        self.cif = self.root / "original.cif"

    def test_swapped_chain_labels_map_by_complete_sequence_and_position(self):
        result = extract_hcdr3_plddt(_cif(self.cif), self.pdb)
        self.assertEqual(result["source_cif_chain"], "L")
        self.assertEqual(result["canonical_chain"], "H")
        self.assertEqual(result["score_source"], "_ma_qa_metric_local.metric_value")
        self.assertEqual([row["residue_number"] for row in result["residues"]], list(RESIDUES))
        self.assertEqual([row["source_seq_id"] for row in result["residues"]], list(range(3, 20)))
        self.assertEqual([row["plddt"] for row in result["residues"]], list(range(63, 80)))
        self.assertEqual(result["values"]["hcdr3_plddt_mean"], 71.0)
        self.assertEqual(result["values"]["hcdr3_plddt_min"], 63.0)
        self.assertEqual(result["values"]["hcdr3_plddt_113"], 79.0)

    def test_local_plddt_must_be_explicitly_declared(self):
        with self.assertRaisesRegex(ValueError, "explicitly declared local pLDDT"):
            extract_hcdr3_plddt(_cif(self.cif, declared=False), self.pdb)

    def test_one_missing_interval_score_is_not_silently_filled_from_b_factor(self):
        with self.assertRaisesRegex(ValueError, "Missing or ambiguous confidence"):
            extract_hcdr3_plddt(_cif(self.cif, omitted_score=3), self.pdb)

    def test_invalid_scores_are_rejected(self):
        for value in ("-0.01", "100.01", "nan", "inf", "?"):
            with self.subTest(value=value), self.assertRaises(ValueError):
                extract_hcdr3_plddt(_cif(self.cif, changed_score=(3, value)), self.pdb)

    def test_duplicate_local_score_is_ambiguous(self):
        with self.assertRaisesRegex(ValueError, "Duplicate pLDDT"):
            extract_hcdr3_plddt(_cif(self.cif, duplicate_score=True), self.pdb)

    def test_identical_source_chains_are_not_arbitrarily_selected(self):
        with self.assertRaisesRegex(ValueError, "full-chain sequence match.*found 2"):
            extract_hcdr3_plddt(_cif(self.cif, ambiguous=True), self.pdb)

    def test_sequence_mismatch_outside_interval_still_rejects_mapping(self):
        names = ("VAL",) + NAMES[1:]
        with self.assertRaisesRegex(ValueError, "full-chain sequence match.*found 0"):
            extract_hcdr3_plddt(_cif(self.cif, source_names=names), self.pdb)

    def test_missing_canonical_heavy_chain_is_rejected(self):
        with self.assertRaisesRegex(ValueError, "Canonical chain H is absent"):
            extract_hcdr3_plddt(_cif(self.cif), _pdb(self.pdb, chain="L"))

    def test_incomplete_target_interval_is_rejected(self):
        with self.assertRaisesRegex(ValueError, "Incomplete canonical H97"):
            extract_hcdr3_plddt(_cif(self.cif), _pdb(self.pdb, numbers=range(94, 113)))

    def test_inserted_target_number_is_ambiguous(self):
        with self.assertRaisesRegex(ValueError, "Missing or ambiguous confidence"):
            extract_hcdr3_plddt(_cif(self.cif), _pdb(self.pdb, insertion="A"))


if __name__ == "__main__":
    unittest.main()
