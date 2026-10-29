"""Architecture report wiring; no detector or scoring execution."""
import json

from fv_energy_analysis.fingerprint_reports import ARCHITECTURE_COLUMNS, _grouped_columns, _saved_architecture
from fv_energy_analysis.spreadsheet import _sections


def test_architecture_is_a_separate_contiguous_section():
    columns = _grouped_columns()
    groups = _sections(columns)
    assert sum(name == "E2 architecture" for name, _, _ in groups) == 1
    assert [c["key"] for c in columns if c["section"] == "E2 architecture"] == [c[0] for c in ARCHITECTURE_COLUMNS]
    labels = [c[1] for c in ARCHITECTURE_COLUMNS]
    assert not any("ray" in label.lower() or "/ 6" in label or "/ 12" in label for label in labels)
    assert "Opposing residue pairs" in labels
    assert "A-ring atom-pair overlaps" in labels


def test_missing_recorded_inputs_are_unavailable(tmp_path):
    result = _saved_architecture(tmp_path, "mutant", tmp_path / "missing.csv")
    assert result["e2_architecture_status"].startswith("unavailable")
    assert "e2_strict_screen_pass" not in result


def test_saved_architecture_uses_recorded_paths_and_ds(monkeypatch, tmp_path):
    from fv_energy_analysis import estrogen_architecture
    final = tmp_path / "final"
    final.mkdir()
    pdb = tmp_path / "ordinary.pdb"
    pdb.touch()
    params = tmp_path / "E2L.params"
    params.write_text("NAME E2L\n")
    template = tmp_path / "E2L_0001.pdb"
    template.touch()
    from fv_energy_analysis.atom_roles import save_atom_roles
    roles = {"o3_atom": "O1", "o17_atom": "O2", "c18_methyl_atom": "C18", "a_ring_atoms": ["C1", "C2", "C3", "C4", "C5", "C6"]}
    save_atom_roles(params, template, {"valid": True, "estrogen_atom_roles": roles})
    ds = tmp_path / "discovery_studio_interactions.csv"
    ds.touch()
    (final / "summary.json").write_text(json.dumps({"mutant": str(pdb), "fv_chains": ["H", "L"]}))
    (tmp_path / "workflow_manifest.json").write_text(json.dumps({"params": str(params)}))
    calls = []
    def analyze(*args, **kwargs):
        calls.append((args, kwargs))
        return {"e2_a_ring_pi_stack_count": 2}
    monkeypatch.setattr(estrogen_architecture, "analyze_estrogen_architecture", analyze)
    result = _saved_architecture(tmp_path, "mutant", ds)
    assert result["e2_architecture_status"] == "complete"
    assert calls[0][0][0] == pdb
    assert calls[0][1]["discovery_studio_csv"] == ds
    assert "params_path" not in calls[0][1]
    assert calls[0][1]["o3_atom"] == "O1"
    assert calls[0][1]["fv_chains"] == {"H", "L"}
