"""Fresh-folder selection and immutable OptH resume/parent-reuse boundaries."""
from copy import deepcopy
import json
from pathlib import Path
from tempfile import TemporaryDirectory
import unittest

from fv_energy_analysis.hydrogen_batch import (
    digest, list_inputs, reuse_case, validate_reusable_parent,
)
import fv_energy_analysis.hydrogen_batch as module


class FreshHydrogenBatchTests(unittest.TestCase):
    def setUp(self):
        self.temporary = TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.root = Path(self.temporary.name)
        self.parent = self.root/"parent.pdb"
        self.parent.write_text("REMARK original parent\nEND\n")
        self.output = self.root/"optimized.pdb"
        self.output.write_text("REMARK optimized parent\nEND\n")
        self.protocol = {"source_protocol": {"weights": "ligand", "legacy_pre_talaris": True,
                         "params_sha256": "p", "weights_sha256": "w", "pyrosetta_version": "version",
                         "ligand_chain": "X"}, "passes_per_complex": 1, "seed": 20261022,
                         "selection": "lowest bound energy of original and one optimized trial"}
        self.record = {"status": "complete", "identity": {"protocol": deepcopy(self.protocol),
                       "pdb_sha256": digest(self.parent), "source_sha256": "raw", "input_ligand_chain": "d"},
                       "optimized_pdb": str(self.output), "artifact_sha256": {str(self.output): digest(self.output)}}
        self.record["identity"]["protocol"]["implementation_sha256"] = {
            name: digest(Path(module.__file__).with_name(name))
            for name in ("hydrogen_benchmark.py", "pyrosetta_score.py", "hbonds.py")}

    def test_selection_can_expand_from_smoke_without_renaming_candidates(self):
        folder = self.root/"inputs"
        folder.mkdir()
        for name in ("A_relaxed.pdb", "B_relaxed.pdb"):
            (folder/name).write_text("END\n")
        self.assertEqual([p.name for p in list_inputs(folder,self.parent,pattern="A*.pdb")], ["A_relaxed.pdb"])
        self.assertEqual([p.name for p in list_inputs(folder,self.parent)], ["A_relaxed.pdb","B_relaxed.pdb"])

    def test_recursive_duplicate_candidate_names_are_rejected(self):
        for directory in (self.root/"inputs/a", self.root/"inputs/b"):
            directory.mkdir(parents=True)
            (directory/"same.pdb").write_text("END\n")
        with self.assertRaisesRegex(ValueError,"unique"):
            list_inputs(self.root/"inputs",self.parent,recursive=True)

    def test_parent_reuse_requires_identical_prepared_parent(self):
        validate_reusable_parent(self.record,self.protocol,self.parent)
        self.parent.write_text("REMARK changed parent\nEND\n")
        with self.assertRaisesRegex(ValueError,"different inputs"):
            validate_reusable_parent(self.record,self.protocol,self.parent)

    def test_parent_reuse_rejects_changed_seed_weights_or_optimization_count(self):
        for change in ("seed","weights","passes"):
            protocol = deepcopy(self.protocol)
            if change == "seed":
                protocol["seed"] += 1
            elif change == "weights":
                protocol["source_protocol"]["weights"] = "ref2015"
            else:
                self.record["identity"]["protocol"]["passes_per_complex"] = 5
            with self.subTest(change=change), self.assertRaisesRegex(ValueError,"different inputs"):
                validate_reusable_parent(self.record,protocol,self.parent)

    def test_parent_reuse_rejects_tampered_optimized_structure(self):
        self.output.write_text("REMARK changed H coordinates\nEND\n")
        with self.assertRaisesRegex(ValueError,"artifact changed"):
            validate_reusable_parent(self.record,self.protocol,self.parent)

    def test_candidate_resume_checks_settings_and_persisted_coordinates(self):
        path = self.root/"result.json"
        path.write_text(json.dumps(self.record))
        expected = {key:self.record["identity"][key] for key in ("protocol","source_sha256","input_ligand_chain")}
        self.assertEqual(reuse_case(path,expected),self.record)
        with self.assertRaisesRegex(ValueError,"incompatible"):
            reuse_case(path,{**expected,"source_sha256":"different"})
        self.output.write_text("REMARK altered\nEND\n")
        with self.assertRaisesRegex(ValueError,"modified"):
            reuse_case(path,expected)


if __name__ == "__main__":
    unittest.main()
