"""Artifact reuse must not silently mix old chemistry or hydrogen poses."""
import json
from hashlib import sha256

import pytest

from fv_energy_analysis.architecture_batch import atom_snapshot, verify_opth


def atom(x=1.0):
    return f"HETATM    1  C10 E2L X   1    {x:8.3f}{2.0:8.3f}{3.0:8.3f}  1.00 10.00           C\n"


def record(tmp_path):
    pdb = tmp_path / "optimized.pdb"
    pdb.write_text(atom())
    data = {"status": "complete", "optimized_pdb": str(pdb),
            "artifact_sha256": {str(pdb): sha256(pdb.read_bytes()).hexdigest()},
            "identity": {"pdb_sha256": "pre-optimization digest is deliberately not used",
                         "protocol": {"seed": 42, "passes_per_complex": 1,
                                      "source_protocol": {"weights": "ligand", "legacy_pre_talaris": True,
                                                          "params_sha256": "params"}}}}
    saved = tmp_path / "result.json"
    saved.write_text(json.dumps(data))
    return pdb, saved, data


def test_verify_matches_named_atoms_not_file_format(tmp_path):
    pdb, saved, _ = record(tmp_path)
    prepared = tmp_path / "prepared.pdb"
    prepared.write_text("REMARK preparation metadata\n" + pdb.read_text() + "END\n")
    result = verify_opth(prepared, saved, "params")
    assert result["all_named_atoms_and_coordinates_identical"]
    assert result["replicates"] == 1


def test_refuse_changed_coordinate_even_if_structure_hash_valid(tmp_path):
    _, saved, _ = record(tmp_path)
    prepared = tmp_path / "prepared.pdb"
    prepared.write_text(atom(1.001))
    with pytest.raises(ValueError, match="differs"):
        verify_opth(prepared, saved, "params")


def test_refuse_changed_params_or_tampered_optimized_file(tmp_path):
    pdb, saved, _ = record(tmp_path)
    with pytest.raises(ValueError, match="protocol"):
        verify_opth(pdb, saved, "different params")
    pdb.write_text(atom(2))
    with pytest.raises(ValueError, match="hash"):
        verify_opth(pdb, saved, "params")


def test_duplicate_named_atoms_rejected(tmp_path):
    pdb = tmp_path / "duplicate.pdb"
    pdb.write_text(atom() + atom(2))
    with pytest.raises(ValueError, match="Duplicate"):
        atom_snapshot(pdb)


def test_resume_reuses_only_hash_matching_completed_architecture(monkeypatch, tmp_path):
    from fv_energy_analysis import architecture_batch as module
    from fv_energy_analysis.ds_architecture import VERSION
    from fv_energy_analysis.steric_sandwich import METHOD_VERSION
    pdb, _, _ = record(tmp_path)
    ds = tmp_path / "ds.csv"
    ds.write_text("unchanged saved DS records")
    (tmp_path / "final").mkdir()
    (tmp_path / "final/summary.json").write_text(json.dumps({"parent": str(pdb)}))
    calibration = {"production_rule": "test identity"}
    cached = {"candidate": "case.pdb", "status": "complete", "result": {
        "e2_architecture_definition_version": VERSION,
        "e2_steric_sandwich": {"method": METHOD_VERSION},
        "e2_c18_calibration": calibration, "provenance": {
            "prepared_pdb_sha256": module.digest(pdb), "ds_csv_sha256": module.digest(ds),
            "optimized_pdb": str(pdb), "optimized_pdb_sha256": module.digest(pdb)}}}
    cache = tmp_path / "cached.json"
    cache.write_text(json.dumps(cached))
    task = {"candidate": "case.pdb", "source_result": str(tmp_path), "role": "parent",
            "ds_csv": str(ds), "cache_path": str(cache), "resume": True,
            "c18_calibration": calibration, "params_hash": "params", "opth_record": "not reached"}
    def refuse(*args):
        raise ValueError("cache miss must revalidate instead of silently reusing")
    monkeypatch.setattr(module, "verify_opth", refuse)
    assert module._one(task)["reused"] is True
    ds.write_text("changed DS record")
    missed = module._one(task)
    assert missed["status"] == "failed"
    assert "cache miss" in missed["error"]


@pytest.mark.parametrize("changed_key", ["architecture", "steric"])
def test_resume_revalidates_old_architecture_definition(monkeypatch, tmp_path, changed_key):
    from fv_energy_analysis import architecture_batch as module
    from fv_energy_analysis.ds_architecture import VERSION
    from fv_energy_analysis.steric_sandwich import METHOD_VERSION
    pdb, _, _ = record(tmp_path)
    ds = tmp_path / "ds.csv"
    ds.write_text("unchanged saved DS records")
    (tmp_path / "final").mkdir()
    (tmp_path / "final/summary.json").write_text(json.dumps({"parent": str(pdb)}))
    calibration = {"production_rule": "test identity"}
    previous = {
        "e2_architecture_definition_version": VERSION,
        "e2_steric_sandwich": {"method": METHOD_VERSION},
        "e2_c18_calibration": calibration,
        "provenance": {"prepared_pdb_sha256": module.digest(pdb),
                       "ds_csv_sha256": module.digest(ds), "optimized_pdb": str(pdb),
                       "optimized_pdb_sha256": module.digest(pdb)},
    }
    if changed_key == "architecture":
        previous["e2_architecture_definition_version"] = "historical-definition"
    else:
        previous["e2_steric_sandwich"]["method"] = METHOD_VERSION + "-obsolete"
    cache = tmp_path / "cached.json"
    cache.write_text(json.dumps({"candidate": "case.pdb", "status": "complete", "result": previous}))
    task = {"candidate": "case.pdb", "source_result": str(tmp_path), "role": "parent",
            "ds_csv": str(ds), "cache_path": str(cache), "resume": True,
            "c18_calibration": calibration, "params_hash": "params", "opth_record": "not reached"}
    def refuse(*args):
        raise ValueError("changed definition must revalidate")
    monkeypatch.setattr(module, "verify_opth", refuse)
    missed = module._one(task)
    assert missed["status"] == "failed"
    assert "changed definition" in missed["error"]


def test_report_describes_enclosure_without_ray_counts(tmp_path):
    from fv_energy_analysis.architecture_batch import render_report
    report = render_report("case.pdb", {}, {}, tmp_path)
    assert "two-sided A-ring enclosure" in report
    assert "0–1.0 A" in report
    assert "projected union" in report
    assert "C19-exclusion test" in report
    assert "carbon-normal paths" not in report
    assert "overlapping rays" not in report


def test_public_architecture_never_reads_arpeggio(monkeypatch, tmp_path):
    from fv_energy_analysis import ds_architecture
    from fv_energy_analysis.estrogen_architecture import analyze_estrogen_architecture
    calls = []
    monkeypatch.setattr(ds_architecture, "analyze_ds_architecture", lambda *a, **k: calls.append((a, k)) or {"ok": True})
    invalid_arpeggio = tmp_path / "invalid_arpeggio.json"
    invalid_arpeggio.write_text("THIS IS NOT JSON")
    result = analyze_estrogen_architecture(
        tmp_path / "not_read_here.pdb", invalid_arpeggio,
        o3_atom="O1", o17_atom="O2", c18_methyl_atom="C10",
        a_ring_atoms=("C1", "C2", "C3", "C4", "C5", "C6"),
        saved_features={"e2_o3_hbond": True},
    )
    assert result == {"ok": True}
    assert len(calls[0][0]) == 2
    assert "saved_features" not in calls[0][1]
    assert "arpeggio_json" not in calls[0][1]


def test_unavailable_ds_not_promoted_to_complete(monkeypatch, tmp_path):
    from fv_energy_analysis import estrogen_architecture
    from fv_energy_analysis.fingerprint_reports import _saved_architecture
    from fv_energy_analysis.atom_roles import save_atom_roles
    (tmp_path / "final").mkdir()
    pdb, params, template = [tmp_path / name for name in ("a.pdb", "E2L.params", "E2L_0001.pdb")]
    for path in (pdb, params, template):
        path.write_text("placeholder")
    roles = {"o3_atom": "O1", "o17_atom": "O2", "c18_methyl_atom": "C10", "a_ring_atoms": list("ABCDEF")}
    save_atom_roles(params, template, {"valid": True, "estrogen_atom_roles": roles})
    (tmp_path / "workflow_manifest.json").write_text(json.dumps({"params": str(params)}))
    (tmp_path / "final/summary.json").write_text(json.dumps({"mutant": str(pdb), "fv_chains": ["H", "L"], "ligand_resname": "E2L"}))
    monkeypatch.setattr(estrogen_architecture, "analyze_estrogen_architecture", lambda *a, **k:
                        {"e2_architecture_status": "unavailable: DS missing", "e2_o3_hbond": None})
    result = _saved_architecture(tmp_path, "mutant", tmp_path / "missing.csv")
    assert result["e2_architecture_status"].startswith("unavailable")
    assert result["e2_o3_hbond"] is None
