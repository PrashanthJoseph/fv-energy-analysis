"""Native detector names, ligand selection and missing-data regressions."""
import json

import pytest

from fv_energy_analysis.fingerprint_counts import (
    _arpeggio_counts, _arpeggio_raw, _delta, _ds_counts, _ligand_resname,
    _plip_counts, _raw_count_rows, _rosetta_counts, _rosetta_energies,
)
from fv_energy_analysis.fingerprint_reports import ESSENTIAL_COUNTS, _grouped_columns, _markdown, _raw_markdown
from fv_energy_analysis.report import write_json


def test_arpeggio_keeps_native_names_and_uses_requested_ligand(tmp_path):
    def record(ligand, contacts):
        return {"bgn": {"label_comp_id": ligand, "label_comp_type": "B"},
                "end": {"label_comp_id": "ARG", "label_comp_type": "P"},
                "type": "atom-atom", "contact": contacts}
    path = tmp_path / "contacts.json"
    path.write_text(json.dumps([record("EST", ["ionic", "hbond", "CATIONPI"]),
                                record("E2L", ["ionic", "ionic"])]))
    counts = _arpeggio_counts(path, "EST")
    assert counts["ionic"] == counts["hbond"] == counts["CATIONPI"] == 1
    assert counts["salt_bridge"] == counts["attractive_charge"] == "NA"
    assert counts["conventional_hbond"] == counts["pi_cation"] == "NA"
    assert _arpeggio_raw(path, "EST") == {"ionic": 1, "hbond": 1, "CATIONPI": 1}
    assert all(label == key for label, key in ESSENTIAL_COUNTS["Arpeggio"])


@pytest.mark.parametrize("summary", [{}, {"ligand_hbonds": {}}, {"ligand_hbonds": {"parent": None}}])
def test_missing_rosetta_counts_are_not_zero(summary):
    counts = _rosetta_counts(summary, "parent")
    assert counts["hbond_total"] is None
    assert _delta(counts["hbond_total"], 2) is None
    assert _rosetta_counts({"ligand_hbonds": {"parent": []}}, "parent")["hbond_total"] == 0


@pytest.mark.parametrize("data", [None, {}, {"interaction_counts": None},
                                  {"interaction_counts": {"hydrogen_bonds": 0}, "binding_site_count": 0}])
def test_missing_plip_counts_or_unmatched_ligand_are_not_zero(data):
    counts = _plip_counts({"plip_ligand_interactions": {"parent": data}}, "parent")
    assert counts["hbond_total"] is None
    assert counts["total"] is None


def test_recorded_plip_zero_is_zero_but_unrecorded_class_is_missing():
    counts = _plip_counts({"plip_ligand_interactions": {"parent": {
        "interaction_counts": {"hydrogen_bonds": 0}, "total_interaction_count": 0,
        "binding_site_count": 1}}}, "parent")
    assert counts["hbond_total"] == counts["total"] == 0
    assert counts["salt_bridge"] is None


def test_ligand_name_comes_from_metadata_or_params_not_default(tmp_path):
    assert _ligand_resname({"ligand_resname": "ABC"}, tmp_path) == "ABC"
    assert _ligand_resname({"plip_ligand_interactions": {"parent": {"ligand_resname": "XYZ"}}}, tmp_path) == "XYZ"
    params = tmp_path / "ligand.params"
    params.write_text("NAME long_internal_residue_name\nIO_STRING EST Z\n")
    (tmp_path / "workflow_manifest.json").write_text(json.dumps({"params": str(params)}))
    assert _ligand_resname({}, tmp_path) == "EST"
    with pytest.raises(ValueError, match="No recorded ligand"):
        _ligand_resname({}, tmp_path / "missing")


def test_ds_pi_stacks_count_unique_residue_partners(tmp_path):
    path = tmp_path / "ds.csv"
    path.write_text("interaction_type,protein_residue,residue_pair_id\n"
                    "pi_pi_stacked,H:99:TRP,H:99:TRP<->X:1:E2L\n"
                    "pi_pi_stacked,H:99:TRP,H:99:TRP<->X:1:E2L\n"
                    "pi_pi_stacked,L:91:TYR,L:91:TYR<->X:1:E2L\n")
    counts = _ds_counts(path)
    assert counts["pi_stack"] == 2
    assert counts["total"] == 3  # Native record total is deliberately unchanged.


def test_ds_conventional_hbond_column_is_not_combined(tmp_path):
    path = tmp_path / "ds.csv"
    path.write_text("interaction_type\nconventional_hbond\ncarbon_hbond\nwater_bridge_hbond\nwater_hbond\n")
    counts = _ds_counts(path)
    assert counts["conventional_hbond"] == 1
    assert counts["hbond_total"] == 4
    assert ("Conventional H-bonds", "conventional_hbond") in ESSENTIAL_COUNTS["Discovery Studio"]
    assert not any(key == "hbond_total" for _, key in ESSENTIAL_COUNTS["Discovery Studio"])


def test_hbond_sc_energy_label_and_missing_energy():
    columns = {column["key"]: column for column in _grouped_columns()}
    assert columns["rosetta_hbond_sc_REU"]["label"] == "hbond_sc only (REU)"
    assert _rosetta_energies({})["parent_hbond_sc_REU"] is None


def test_missing_counts_survive_raw_reports_and_strict_json(tmp_path):
    for role in ("parent", "mutant"):
        directory = tmp_path / f"arpeggio_{role}"
        directory.mkdir()
        (directory / f"{role}_arpeggio.json").write_text("[]")
    parent_ds = tmp_path / "parent_ds.csv"
    parent_ds.write_text("interaction_type\n")
    (tmp_path / "discovery_studio_interactions.csv").write_text("interaction_type\n")
    summary = {"ligand_resname": "EST", "ligand_hbonds": {"parent": []},
               "plip_ligand_interactions": {"mutant": {"interaction_counts": {"hydrogen_bonds": 2}}}}
    rows = _raw_count_rows("candidate", tmp_path, summary, parent_ds)
    plip = next(row for row in rows if row["tool"] == "PLIP" and row["detector_interaction_type"] == "hydrogen_bonds")
    assert plip["parent_count"] is None
    assert plip["mutant_count"] == 2
    assert plip["delta_count"] is None
    assert "NaN" in _raw_markdown(rows)
    output = tmp_path / "counts.json"
    write_json(output, rows)
    assert "null" in output.read_text()
    assert "NaN" not in output.read_text()
