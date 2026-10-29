"""Refresh only the E3 classification and the selected new hydroxyl proton."""
from hashlib import sha256
from pathlib import Path

import numpy as np
import pytest

from fv_energy_analysis.e3_io import replace_hydroxyl_hydrogen
from fv_energy_analysis.e3_refresh import refresh_one


def test_hydrogen_persistence_preserves_other_bytes():
    source = ("REMARK original\n"
              "HETATM    1  H25 E3L X   1       1.000   2.000   3.000  1.00  0.00           H\n"
              "ATOM      2  CA  ALA H   1       4.000   5.000   6.000  1.00  0.00           C\nEND\n").encode()
    changed = replace_hydroxyl_hydrogen(source, [1.5, 2.5, 3.5])
    for old, new in zip(source.splitlines(), changed.splitlines()):
        if old[12:16].strip() != b"H25":
            assert old == new
        else:
            assert old[:30] == new[:30] and old[54:] == new[54:]
    with pytest.raises(ValueError, match="finite"):
        replace_hydroxyl_hydrogen(source, [float("nan"), 0, 0])


ROOT = Path(__file__).resolve().parents[1]
STEM = "set_A_3mut_type4_contiguous_0_anchored_0_E1Y_R2T_P5D_estradiol_model_0_relaxed"
SAVED = ROOT / "results/full_2436_e3_steric_20261027/cases" / STEM / "e3_discrimination.json"


@pytest.mark.skipif(not SAVED.is_file(), reason="Local regression dataset absent")
def test_backbone_trp_overlap_no_longer_receives_false_clearance(tmp_path):
    import json
    from fv_energy_analysis.paths import native_path
    previous = json.loads(SAVED.read_text())
    probe = native_path(previous["probe_pdb"])
    before = probe.read_bytes()
    result = refresh_one({"candidate": STEM + ".pdb", "artifact": str(SAVED),
                          "output": str(tmp_path / "corrected"), "identity": {"test": True},
                          "overlap_warning_A": .4, "resume": False})
    assert result["status"] == "complete", result.get("error")
    assert result["summary"]["e3_steric_warning"] is True
    assert result["summary"]["e3_added_o_non_hbond_overlap_residue_count"] == 1
    contact = next(c for c in result["contacts"] if c["protein_atom"] == "L:35:TRP:N")
    assert contact["hbond_geometry"]["status"] == "not_supported"
    trials = contact["hbond_geometry"]["trials"]
    assert any(t["angle_DHA_deg"] == pytest.approx(83.1665752) for t in trials)
    assert all(t["passes_geometry"] is False for t in trials)
    assert probe.read_bytes() == before
    assert result["source_e3_probe_sha256"] == sha256(before).hexdigest()
    after = Path(result["probe_pdb"]).read_bytes()
    for old, new in zip(before.splitlines(), after.splitlines()):
        if old.startswith(b"HETATM") and old[12:16].strip() == b"H25":
            assert old[:30] == new[:30] and old[54:] == new[54:]
        else:
            assert old == new
    line = next(l for l in after.decode().splitlines() if l.startswith("HETATM") and l[12:16].strip() == "H25")
    xyz = [float(line[i:i+8]) for i in (30, 38, 46)]
    np.testing.assert_allclose(xyz, result["selected_hydrogen_xyz_A"], atol=1e-9)
