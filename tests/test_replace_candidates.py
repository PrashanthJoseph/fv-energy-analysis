"""Recovery replaces exactly the rejected identities without rescoring old rows."""
from copy import deepcopy
import csv
from hashlib import sha256
import json

import pytest

from fv_energy_analysis import replace_candidates as module
from fv_energy_analysis.ds_architecture import VERSION
from fv_energy_analysis.e3_reports import e3_columns, flatten_e3_result
from fv_energy_analysis.report import write_json
from fv_energy_analysis.steric_sandwich import METHOD_VERSION


SCIENCE = [
    {"key": "rank", "label": "Rank", "section": "Identity", "kind": "integer"},
    {"key": "candidate", "label": "Complex", "section": "Identity", "kind": "text"},
    {"key": "rosetta_interface_REU", "label": "Interface", "section": "RosettaLigand", "kind": "energy"},
    {"key": "rosetta_interface_delta_REU", "label": "Delta", "section": "RosettaLigand", "kind": "delta_energy"},
    {"key": "e2_geometric_sandwich", "label": "Enclosure", "section": "E2 architecture", "kind": "boolean"},
]
STATUS = [{"key": "status", "label": "Analysis status", "section": "Analysis status", "kind": "text"},
          {"key": "error", "label": "Failure reason", "section": "Analysis status", "kind": "text"}]


def e3(candidate):
    return {"schema_version": 2, "status": "complete", "candidate": candidate,
            "summary": {"e3_analysis_scope": "fixed_pocket_steric_only",
                "e3_added_o_closest_distance_A": 3.5, "e3_added_o_closest_residue": "H:1:ALA",
                "e3_added_o_deepest_overlap_A": 0., "e3_added_o_overlap_residue_count": 0,
                "e3_added_o_non_hbond_overlap_residue_count": 0,
                "e3_added_o_hbond_supported_overlap_residue_count": 0,
                "e3_added_o_unassessed_overlap_residue_count": 0,
                "e3_added_o_overlap_residues": "none", "e3_steric_assessment": "no_large_overlap_in_fixed_pose",
                "e3_steric_warning": False, "e3_oh_h_best_sampled_protein_overlap_A": 0.},
            "protocol": {"method": "geometry-v2", "scope": "fixed_pocket_steric_only",
                "warning_overlap_A": .4, "radii_A": {"C": 1.7, "O": 1.52},
                "hbond_geometry_criteria": {"angle_DHA_degrees": "> 100"},
                "geometry_combination_policy": "one shared H orientation"}}


def row(name, delta, rank=None, failed=False):
    result = {"rank": rank, "candidate": name, "rosetta_interface_REU": None if failed else -20+delta,
              "rosetta_interface_delta_REU": None if failed else delta,
              "e2_geometric_sandwich": None if failed else True,
              "status": "failed" if failed else "complete", "error": "stereo mismatch" if failed else None}
    result.update(flatten_e3_result(None, missing_status="excluded_input_rejection") if failed
                  else flatten_e3_result(e3(name)))
    return result


def test_identity_recovery_preserves_old_values_and_recalculates_global_rank():
    base = [row("old.pdb", -1, 1), row("fixed-b.pdb", 0, failed=True), row("fixed-c.pdb", 0, failed=True)]
    previous = deepcopy(base)
    new = [row("fixed-b.pdb", -3, 1), row("fixed-c.pdb", 0, 2)]
    merged, replaced = module.merge_rows(base, new, {r["candidate"]: e3(r["candidate"]) for r in new})
    assert replaced == {"fixed-b.pdb", "fixed-c.pdb"}
    assert [r["candidate"] for r in merged] == ["fixed-b.pdb", "old.pdb", "fixed-c.pdb"]
    assert [r["rank"] for r in merged] == [1, 2, 3]
    assert base == previous
    assert {k: v for k, v in merged[1].items() if k != "rank"} == {k: v for k, v in base[0].items() if k != "rank"}
    assert all(r["status"] == "complete" and r["error"] is None for r in merged)


@pytest.mark.parametrize("names", [["wrong.pdb"], ["old.pdb"], ["fixed.pdb", "fixed.pdb"], []])
def test_missing_unexpected_duplicate_or_already_completed_replacement_is_rejected(names):
    base = [row("old.pdb", -1, 1), row("fixed.pdb", 0, failed=True)]
    new = [row(name, -2) for name in names]
    with pytest.raises(ValueError, match="Duplicate|inventory differs"):
        module.merge_rows(base, new, {name: e3(name) for name in names})


def test_missing_or_failed_e3_cannot_become_a_completed_row():
    base, new = [row("fixed.pdb", 0, failed=True)], [row("fixed.pdb", -1)]
    with pytest.raises(ValueError, match="E3 inventory"):
        module.merge_rows(base, new, {})
    with pytest.raises(ValueError, match="E3 analysis is not complete"):
        module.merge_rows(base, new, {"fixed.pdb": {"candidate": "fixed.pdb", "status": "failed"}})


def test_parent_imported_annotations_may_be_reused_but_scientific_mismatch_is_rejected():
    columns = [*SCIENCE, {"key": "plddt", "label": "pLDDT", "section": "HCDR3 pLDDT", "kind": "number"}]
    base = {"rosetta_interface_REU": -20, "rosetta_interface_delta_REU": 0, "e2_geometric_sandwich": True, "plddt": 80}
    other = {**base, "plddt": None}
    module._compare_parent(base, other, columns)
    other["rosetta_interface_REU"] = -19
    with pytest.raises(ValueError, match="parent baseline differs"):
        module._compare_parent(base, other, columns)


def pdb_atom(serial, name, element, residue="E3L", chain="X", number=1, record="HETATM"):
    return (f"{record:<6}{serial:5d} {name:>4} {residue:>3} {chain}{number:4d}    "
            f"{serial:8.3f}{0:8.3f}{0:8.3f}{1:6.2f}{0:6.2f}          {element:>2}\n")


def probe_and_input(root, candidate):
    root.mkdir(parents=True, exist_ok=True)
    prepared = root / "prepared.pdb"
    prepared.write_text(pdb_atom(1, "CA", "C", "ALA", "H", 1, "ATOM"))
    probe = root / "e3_probe.pdb"
    ligand = [(f"C{i}", "C") for i in range(1, 19)] + [(f"O{i}", "O") for i in range(1, 4)] + [(f"H{i}", "H") for i in range(1, 25)]
    probe.write_text(prepared.read_text()+"".join(pdb_atom(i+2, name, element) for i, (name, element) in enumerate(ligand)))
    result = e3(candidate)
    result.update(probe_pdb=str(probe), probe_pdb_sha256=module.digest(probe),
                  input_pdb=str(prepared), input_pdb_sha256=module.digest(prepared),
                  chemistry={"inchikey": "PROQIPRRNZUXQM-ZXXIGWHRSA-N", "formula": "C18H24O3", "formal_charge": 0, "ligand_hydrogen_count": 24})
    return result, prepared


def test_e3_must_match_the_scored_structure_and_saved_protonated_probe(tmp_path):
    result, prepared = probe_and_input(tmp_path, "case.pdb")
    context = module._geometry_context(result)
    module._validate_e3(result, "case.pdb", prepared, context)
    changed = tmp_path / "different.pdb"
    changed.write_text(prepared.read_text().replace("   1.000", "   2.000"))
    with pytest.raises(ValueError, match="differs from the scored"):
        module._validate_e3(result, "case.pdb", changed, context)
    result["probe_pdb_sha256"] = "bad"
    with pytest.raises(ValueError, match="hash differs"):
        module._validate_e3(result, "case.pdb", prepared, context)


@pytest.fixture
def complete_sources(tmp_path, monkeypatch):
    base, replacement, new_e3 = [tmp_path / name for name in ("base", "replacement", "new_e3")]
    for folder in (base, replacement, new_e3):
        folder.mkdir()
    old_rows = [row("old.pdb", -1, 1), row("fixed.pdb", 0, failed=True)]
    columns = [*SCIENCE, *e3_columns(), *STATUS]
    parent = {**row("Parent reference", 0), "rank": None, "status": "reference"}
    payload = {"columns": columns, "rows": old_rows, "parent": parent,
               "scoring_profile": "rosettaligand_legacy", "scoring": {"weights": "ligand", "legacy_pre_talaris": True},
               "e3_discrimination_protocol": {"hydrogen_bond_geometry_assessed": True},
               "e2_recognition_protocol": {"definition_version": VERSION, "steric_method": METHOD_VERSION, "interaction_source": "Discovery Studio"}}
    base_csv = module._write_table(base, payload)
    base_e3 = base / "e3_manifest.json"
    write_json(base_e3, {"parent": e3("parent_m7")})
    context = {"weights": "ligand", "parent_hash": "shared"}
    manifest = {"status": "complete_with_failures", "candidate_count": 2, "completed_count": 1, "failure_count": 1,
                "scientific_protocol": context, "csv": str(base_csv), "csv_sha256": module.digest(base_csv),
                "e3_discrimination": str(base_e3), "sources": [{"source": str(base), "candidate_count": 2, "completed_count": 1, "failure_count": 1}]}
    write_json(base / "combined_manifest.json", manifest)
    with (base / "candidate_sources.csv").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=["candidate", "source", "source_rank"])
        writer.writeheader()
        writer.writerows({"candidate": r["candidate"], "source": str(base), "source_rank": r["rank"]} for r in old_rows)
    case = replacement / "case"
    case.mkdir()
    replacement_e3_result, prepared = probe_and_input(new_e3 / "cases/fixed", "fixed.pdb")
    write_json(new_e3 / "cases/fixed/e3_discrimination.json", replacement_e3_result)
    (case / "final").mkdir()
    write_json(case / "final/summary.json", {"mutant": str(prepared), "parent": str(prepared)})
    new_row = row("fixed.pdb", -3, 1)
    new_payload = {**deepcopy(payload), "columns": SCIENCE, "rows": [{c["key"]: new_row[c["key"]] for c in SCIENCE}],
                   "parent": {c["key"]: parent[c["key"]] for c in SCIENCE}}
    new_payload.pop("e3_discrimination_protocol")
    module._write_table(replacement, new_payload)
    write_json(replacement / "batch_manifest.json", {"status": "complete", "candidate_count": 1, "completed_count": 1, "failure_count": 0,
               "candidates": [{"candidate": "fixed.pdb", "status": "complete", "result_folder": str(case)}]})
    write_json(new_e3 / "e3_batch_manifest.json", {"status": "complete", "status_counts": {"complete": 1}, "parent": None})
    monkeypatch.setattr(module, "_scientific_context", lambda *args: context)
    return base, replacement, new_e3, payload


def test_replacement_writes_complete_bound_snapshot_and_updates_provenance(complete_sources, tmp_path):
    base, replacement, new_e3, original = complete_sources
    before = (base / f"{module.TABLE}.json").read_bytes()
    output = tmp_path / "final"
    csv_path = module.replace_failed_candidates(base, replacement, output, new_e3, 2)
    result = module.read_json(csv_path.with_suffix(".json"))
    assert [r["candidate"] for r in result["rows"]] == ["fixed.pdb", "old.pdb"]
    assert [r["rank"] for r in result["rows"]] == [1, 2]
    assert result["rows"][1]["e3_steric_warning"] is False
    assert result["parent"] == original["parent"]
    manifest = module.read_json(output / "combined_manifest.json")
    assert manifest["status"] == "complete" and manifest["failure_count"] == 0
    assert sum(s["candidate_count"] for s in manifest["sources"]) == 2
    assert manifest["replacement"]["replaced_candidates"] == ["fixed.pdb"]
    source = result["e3_discrimination_protocol"]["source_combined_batch"]
    snapshot = module.read_json(output / "pre_e3_combined" / f"{module.TABLE}.json")
    assert source == str(output / "pre_e3_combined")
    assert not any(c["section"] == "E3 discrimination" for c in snapshot["columns"])
    assert snapshot["rows"][0]["status"] == "complete"
    with (output / "candidate_sources.csv").open(encoding="utf-8-sig") as stream:
        sources = {r["candidate"]: r["source"] for r in csv.DictReader(stream)}
    assert sources == {"fixed.pdb": str(replacement), "old.pdb": str(base)}
    assert (base / f"{module.TABLE}.json").read_bytes() == before


def test_replacement_rejects_protocol_mismatch_before_creating_output(complete_sources, tmp_path, monkeypatch):
    base, replacement, new_e3, _ = complete_sources
    monkeypatch.setattr(module, "_scientific_context", lambda *args: {"different": True})
    output = tmp_path / "bad"
    with pytest.raises(ValueError, match="scientific protocol"):
        module.replace_failed_candidates(base, replacement, output, new_e3)
    assert not output.exists()
