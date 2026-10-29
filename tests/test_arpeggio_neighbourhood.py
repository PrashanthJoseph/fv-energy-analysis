import json
import os
from pathlib import Path
from unittest.mock import Mock

import pytest

from fv_energy_analysis.arpeggio import (
    filter_arpeggio_contacts, mutation_neighbourhood_selection,
    run_arpeggio, validate_arpeggio_coverage,
)


def write_structure(path, mutant=False):
    residues = [(99, "PHE" if mutant else "GLU", 0.0),
                (10, "ALA", 12.0 if mutant else 5.0),
                (11, "SER", 5.0 if mutant else 12.0), (12, "GLY", 25.0)]
    lines = [f"ATOM  {i:5d}  CA  {name:3s} H{number:4d}    {x:8.3f}{0.:8.3f}{0.:8.3f}  1.00 20.00           C\n"
             for i, (number, name, x) in enumerate(residues, 1)]
    lines.append(f"HETATM{5:5d}  C1  E2L X{1:4d}    {90.:8.3f}{0.:8.3f}{0.:8.3f}  1.00 20.00           C\n")
    path.write_text("".join(lines) + "END\n")


def test_neighbours_union_both_states_and_keep_ligand(tmp_path):
    parent, mutant = tmp_path / "parent.pdb", tmp_path / "mutant.pdb"
    write_structure(parent)
    write_structure(mutant, True)
    plan = mutation_neighbourhood_selection(parent, mutant, ["H"], "X")
    assert plan["mutated_residues"] == ["H:99"]
    assert plan["selected_residues"] == ["H:10", "H:11", "H:99", "X:1"]
    assert plan["selectors"] == ["/H/10/", "/H/11/", "/H/99/", "/X/1/"]
    assert plan["full_structure_retained"]
    # No mutation: ligand-only query, never an empty/full-structure fallback.
    assert mutation_neighbourhood_selection(parent, parent, ["H"], "X")["selectors"] == ["/X/1/"]


@pytest.mark.parametrize("radius", [0, -1, float("nan"), float("inf")])
def test_invalid_neighbourhood_radius_rejected(tmp_path, radius):
    with pytest.raises(ValueError, match="radius"):
        mutation_neighbourhood_selection(tmp_path / "p", tmp_path / "m", ["H"], "X", radius)


def test_missing_or_mismatched_chain_rejected(tmp_path):
    parent, mutant = tmp_path / "p.pdb", tmp_path / "m.pdb"
    write_structure(parent)
    write_structure(mutant, True)
    with pytest.raises(ValueError, match="matching residue"):
        mutation_neighbourhood_selection(parent, mutant, ["L"], "X")
    mutant.write_text(mutant.read_text().replace("H  12", "H  13"))
    with pytest.raises(ValueError, match="matching residue"):
        mutation_neighbourhood_selection(parent, mutant, ["H"], "X")


def test_query_file_keeps_intact_input_and_fixed_seed(tmp_path):
    runner = Mock()
    cif = tmp_path / "full_structure.cif"
    cif.write_text("unchanged input")
    directory = tmp_path / "contacts"
    run_arpeggio(runner, "pdbe-arpeggio", directory, cif, selectors=["/H/99/", "/X/1/"])
    args, kwargs = runner.run.call_args
    command = args[0]
    assert command[command.index("-sf") + 1] == str(directory / "selection.txt")
    assert (directory / "selection.txt").read_text() == "/H/99/\n/X/1/\n"
    assert command[-1] == str(cif) and cif.read_text() == "unchanged input"
    assert kwargs["env_overrides"] == {"PYTHONHASHSEED": "0"}
    with pytest.raises(ValueError, match="empty"):
        run_arpeggio(runner, "pdbe-arpeggio", directory, cif, selectors=[])


def record(a, b):
    return {"bgn": {"auth_asym_id": "H", "auth_seq_id": a},
            "end": {"auth_asym_id": "H", "auth_seq_id": b}, "contact": ["hbond"]}


def test_scope_filter_is_symmetric_and_does_not_mutate_linked_cache(tmp_path):
    source, linked = tmp_path / "cache.json", tmp_path / "child.json"
    source.write_text(json.dumps([record(99, 10), record(11, 99), record(1, 2)]))
    original = source.read_bytes()
    os.link(source, linked)
    filter_arpeggio_contacts(linked, linked, ["H:99"])
    assert source.read_bytes() == original
    assert json.loads(linked.read_text()) == [record(99, 10), record(11, 99)]


def test_local_reuse_requires_coverage_full_outputs_remain_usable(tmp_path):
    manifest = tmp_path / "workflow_manifest.json"
    manifest.write_text(json.dumps({"arpeggio_selection": {"scope": "neighbourhood", "selected_residues": ["H:99", "X:1"]}}))
    validate_arpeggio_coverage(tmp_path, {"selected_residues": ["X:1"]})
    with pytest.raises(ValueError, match="does not cover"):
        validate_arpeggio_coverage(tmp_path, {"selected_residues": ["H:100"]})
    with pytest.raises(ValueError, match="does not cover"):
        validate_arpeggio_coverage(tmp_path, None)
    manifest.write_text(json.dumps({"arpeggio_selection": {"scope": "full"}}))
    validate_arpeggio_coverage(tmp_path, {"selected_residues": ["H:100"]})
