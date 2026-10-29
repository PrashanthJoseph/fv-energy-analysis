import csv
from pathlib import Path
import tempfile
import unittest

from fv_energy_analysis.rescore_network import _geometry, comparison_warning


class RescoreNetworkTests(unittest.TestCase):
    def test_quality_warning_uses_current_totals_and_repulsion(self):
        parent = {"bound_total": -800.0, "if_X_fa_rep": 1.0}
        self.assertFalse(comparison_warning(parent, {"bound_total": -802.0, "if_X_fa_rep": 2.0}, "X"))
        self.assertTrue(comparison_warning(parent, {"bound_total": -650.0, "if_X_fa_rep": 2.0}, "X"))
        self.assertTrue(comparison_warning(parent, {"bound_total": -802.0, "if_X_fa_rep": 12.0}, "X"))

    def test_geometry_filters_to_mutation_centres(self):
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / "geometry.csv"
            with path.open("w", newline="") as handle:
                writer = csv.DictWriter(handle, fieldnames=["pair_id", "residue1", "residue2", "interaction_types"])
                writer.writeheader()
                writer.writerows([
                    {"pair_id": "H:99--X:1", "residue1": "H:99", "residue2": "X:1", "interaction_types": "hbond;proximal"},
                    {"pair_id": "H:1--H:2", "residue1": "H:1", "residue2": "H:2", "interaction_types": "vdw"},
                ])
            rows = _geometry(path, {"H:99"})
            self.assertEqual([row["pair_id"] for row in rows], ["H:99--X:1"])

    def test_geometry_rejects_old_scoring_labels(self):
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / "geometry.csv"
            path.write_text("pair_id,residue1,residue2,interaction_types\nH:99--X:1,H:99,X:1,rosetta_hbond\n")
            with self.assertRaisesRegex(ValueError, "score-dependent"):
                _geometry(path, {"H:99"})


if __name__ == "__main__":
    unittest.main()
