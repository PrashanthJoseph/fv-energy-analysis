"""Combine saved cohorts without hidden omissions or scientific recomputation."""
import csv
from hashlib import sha256
import json
from pathlib import Path

import pytest

from fv_energy_analysis.combine_batches import combine
from fv_energy_analysis.ds_architecture import VERSION
from fv_energy_analysis.report import write_csv, write_json
from fv_energy_analysis.steric_sandwich import METHOD_VERSION


COLUMNS = [
    {"key": "rank", "label": "Rank", "section": "Identity", "kind": "integer"},
    {"key": "candidate", "label": "Complex", "section": "Identity", "kind": "text"},
    {"key": "rosetta_interface_REU", "label": "Interface", "section": "RosettaLigand", "kind": "energy"},
    {"key": "rosetta_interface_delta_REU", "label": "Delta", "section": "RosettaLigand", "kind": "delta_energy"},
    {"key": "e2_geometric_sandwich", "label": "Enclosure", "section": "E2 architecture", "kind": "boolean"},
]


def read(path):
    return json.loads(path.read_text())


def cohort(tmp_path, label, candidates):
    shared = tmp_path / "shared"
    shared.mkdir(exist_ok=True)
    parent = shared / "hydrogen_optimized_structures" / "parent_m7.pdb"
    parent.parent.mkdir(exist_ok=True)
    parent.write_text("unchanged parent coordinates\n")
    params, template = shared / "E2L.params", shared / "E2L_0001.pdb"
    params.write_text("same validated ligand definition\n")
    template.write_text("same ligand topology template\n")
    params_hash = sha256(params.read_bytes()).hexdigest()
    parent_record = shared / "parent" / "result.json"
    parent_record.parent.mkdir(exist_ok=True)
    write_json(parent_record, {"status": "complete", "optimized_pdb": str(parent),
        "artifact_sha256": {str(parent): sha256(parent.read_bytes()).hexdigest()},
        "identity": {"protocol": {"passes_per_complex": 1, "seed": 42, "source_protocol": {
            "weights": "ligand", "legacy_pre_talaris": True, "params_sha256": params_hash,
            "pyrosetta_version": "test-build", "weights_sha256": "same-weights"}}}})
    source = tmp_path / label
    source.mkdir()
    workflow = source / "case"
    workflow.mkdir()
    write_json(workflow / "workflow_manifest.json", {
        "parent": str(parent), "params": str(params), "ligand_template_pdb": str(template),
        "scoring_protocol": {"weights": "ligand", "params_sha256": params_hash}})
    rows = [{"rank": i, "candidate": name, "rosetta_interface_REU": -20 + delta,
             "rosetta_interface_delta_REU": delta, "e2_geometric_sandwich": True}
            for i, (name, delta) in enumerate(candidates, 1)]
    payload = {"columns": COLUMNS, "rows": rows,
        "parent": {"rank": None, "candidate": "Parent reference", "rosetta_interface_REU": -20,
                   "rosetta_interface_delta_REU": 0, "e2_geometric_sandwich": True},
        "scoring_profile": "rosettaligand_legacy", "scoring": {"weights": "ligand", "legacy_pre_talaris": True},
        "e2_recognition_protocol": {"definition_version": VERSION, "steric_method": METHOD_VERSION,
                                     "interaction_source": "Discovery Studio"}}
    write_json(source / "batch_summary_grouped_fingerprints.json", payload)
    write_csv(source / "batch_summary_grouped_fingerprints.csv", rows)
    write_json(source / "batch_manifest.json", {
        "status": "complete", "candidate_count": len(rows), "completed_count": len(rows), "failure_count": 0,
        "parent": str(parent), "candidates": [{"candidate": row["candidate"], "status": "complete",
                                                "result_folder": str(workflow)} for row in rows],
        "discovery_studio": {"status": "complete", "add_hydrogens": False, "build_info": "same-build"}})
    return source


def change_payload(source, change):
    path = source / "batch_summary_grouped_fingerprints.json"
    payload = read(path)
    change(payload)
    write_json(path, payload)
    write_csv(source / "batch_summary_grouped_fingerprints.csv", payload["rows"])


def test_combines_disjoint_cohorts_recomputes_ranks_and_preserves_metrics(tmp_path):
    a = cohort(tmp_path, "a", [("a.pdb", 2), ("z.pdb", -3)])
    b = cohort(tmp_path, "b", [("b.pdb", -3)])
    original = (a / "batch_summary_grouped_fingerprints.json").read_bytes()
    expected = tmp_path / "expected.json"
    expected.write_text(json.dumps(["z.pdb", "a.pdb", "b.pdb"]))
    csv_path = combine([a, b], tmp_path / "combined", 3, expected)
    result = read(csv_path.with_suffix(".json"))
    assert [row["candidate"] for row in result["rows"]] == ["b.pdb", "z.pdb", "a.pdb"]
    assert [row["rank"] for row in result["rows"]] == [1, 2, 3]
    assert [row["rosetta_interface_delta_REU"] for row in result["rows"]] == [-3, -3, 2]
    assert (a / "batch_summary_grouped_fingerprints.json").read_bytes() == original
    manifest = read(csv_path.parent / "combined_manifest.json")
    assert manifest["status"] == "complete"
    assert manifest["candidate_count"] == 3
    assert manifest["scoring_rerun"] is False and manifest["detectors_rerun"] == []
    assert manifest["scientific_protocol"]["opth_passes"] == 1


def test_same_basename_cannot_silently_duplicate_across_folders_or_case(tmp_path):
    a = cohort(tmp_path, "a", [("folder/A.pdb", -1)])
    b = cohort(tmp_path, "b", [("D:\\other\\a.pdb", -1)])
    with pytest.raises(ValueError, match="Duplicate candidate"):
        combine([a, b], tmp_path / "combined")
    assert not (tmp_path / "combined").exists()


@pytest.mark.parametrize("mismatch", ["schema", "parent", "architecture", "score"])
def test_refuses_incompatible_scientific_tables(tmp_path, mismatch):
    a = cohort(tmp_path, "a", [("a.pdb", -1)])
    b = cohort(tmp_path, "b", [("b.pdb", -2)])
    def change(payload):
        if mismatch == "schema":
            payload["columns"][0]["label"] = "Changed schema"
        elif mismatch == "parent":
            payload["parent"]["rosetta_interface_REU"] = -21
        elif mismatch == "architecture":
            payload["e2_recognition_protocol"]["steric_method"] = "retired-rays"
        else:
            payload["scoring"]["weights"] = "ref2015"
    change_payload(b, change)
    with pytest.raises(ValueError):
        combine([a, b], tmp_path / "combined")


def test_source_csv_and_json_must_agree(tmp_path):
    a = cohort(tmp_path, "a", [("a.pdb", -1)])
    with (a / "batch_summary_grouped_fingerprints.csv").open(newline="") as handle:
        rows = list(csv.DictReader(handle))
    rows[0]["rosetta_interface_delta_REU"] = "-99"
    write_csv(a / "batch_summary_grouped_fingerprints.csv", rows)
    with pytest.raises(ValueError, match="CSV/JSON value differs"):
        combine([a], tmp_path / "combined")


def test_reconciles_expected_count_and_exact_identity_inventory(tmp_path):
    a = cohort(tmp_path, "a", [("a.pdb", -1)])
    with pytest.raises(ValueError, match="Expected 2"):
        combine([a], tmp_path / "combined", 2)
    expected = tmp_path / "expected.json"
    expected.write_text('["different.pdb"]')
    with pytest.raises(ValueError, match="inventory differs"):
        combine([a], tmp_path / "combined", 1, expected)


def test_resolves_architecture_only_refresh_to_original_scoring_provenance(tmp_path):
    original = cohort(tmp_path, "original", [("a.pdb", -1)])
    refreshed = tmp_path / "refreshed"
    refreshed.mkdir()
    payload = read(original / "batch_summary_grouped_fingerprints.json")
    payload["architecture_protocol"] = {"version": VERSION, "steric_method": METHOD_VERSION,
        "interaction_source": "Discovery Studio", "source_batch": str(original)}
    payload.pop("e2_recognition_protocol")
    write_json(refreshed / "batch_summary_grouped_fingerprints.json", payload)
    write_csv(refreshed / "batch_summary_grouped_fingerprints.csv", payload["rows"])
    write_json(refreshed / "batch_manifest.json", {"status": "complete", "candidate_count": 1,
        "completed_count": 1, "failure_count": 0, "source_batch": str(original)})
    result = combine([refreshed], tmp_path / "combined", 1)
    assert result.is_file()


def test_parent_coordinate_hash_mismatch_rejected(tmp_path):
    a = cohort(tmp_path, "a", [("a.pdb", -1)])
    record = tmp_path / "shared" / "parent" / "result.json"
    value = read(record)
    value["artifact_sha256"][value["optimized_pdb"]] = "tampered"
    write_json(record, value)
    with pytest.raises(ValueError, match="Parent OptH protocol/artifact mismatch"):
        combine([a], tmp_path / "combined")


def test_failures_only_included_explicitly_without_inventing_metrics(tmp_path):
    a = cohort(tmp_path, "a", [("good.pdb", -1)])
    manifest_path = a / "batch_manifest.json"
    manifest = read(manifest_path)
    manifest.update(status="complete_with_failures", failure_count=1, candidate_count=2)
    manifest["candidates"].append({"candidate": "bad.pdb", "status": "failed"})
    write_json(manifest_path, manifest)
    write_csv(a / "batch_summary.csv", [
        {"candidate": "good.pdb", "status": "complete"},
        {"candidate": "bad.pdb", "status": "failed", "error": "Ligand stereochemistry mismatch"}])
    with pytest.raises(ValueError, match="complete without failures"):
        combine([a], tmp_path / "strict")
    path = combine([a], tmp_path / "explicit", 2, include_failures=True)
    payload = read(path.with_suffix(".json"))
    assert payload["rows"][1] == {"candidate": "bad.pdb", "status": "failed",
                                  "error": "Ligand stereochemistry mismatch", "rank": None}
    assert payload["rows"][0]["status"] == "complete"
    assert payload["columns"][-2]["key"] == "status"
    result = read(path.parent / "combined_manifest.json")
    assert result["status"] == "complete_with_failures"
    assert result["completed_count"] == 1 and result["failure_count"] == 1
    with path.open(newline="") as handle:
        row = list(csv.DictReader(handle))[1]
    assert row["rosetta_interface_delta_REU"] == ""


def all_rejected_cohort(tmp_path, reference):
    source = tmp_path / "rejected"
    source.mkdir()
    reference_manifest = read(reference / "batch_manifest.json")
    workflow = read(reference / "case" / "workflow_manifest.json")
    parent = Path(reference_manifest["parent"])
    record = read(parent.parent.parent / "parent" / "result.json")
    opt = record["identity"]["protocol"]
    opt["source_protocol"]["template_sha256"] = sha256(Path(workflow["ligand_template_pdb"]).read_bytes()).hexdigest()
    cases = [{"candidate": name, "status": "failed", "error": "Ligand stereochemistry mismatch"}
             for name in ("bad-a.pdb", "bad-b.pdb")]
    opth = source / "hydrogen_batch_manifest.json"
    write_json(opth, {"status": "complete_with_failures", "candidate_count": 2, "completed_count": 0,
        "failure_count": 2, "protocol": opt, "cases": cases, "parent_optimized_pdb": str(parent)})
    manifest = {"status": "complete_with_failures", "candidate_count": 2, "completed_count": 0,
        "failure_count": 2, "preparation_only": True, "reference_batch": str(reference),
        "source_opth_manifest": str(opth), "parent": str(parent), "params": workflow["params"],
        "ligand_template_pdb": workflow["ligand_template_pdb"], "candidates": cases}
    write_json(source / "batch_manifest.json", manifest)
    payload = read(reference / "batch_summary_grouped_fingerprints.json")
    payload["rows"] = []
    write_json(source / "batch_summary_grouped_fingerprints.json", payload)
    with (source / "batch_summary_grouped_fingerprints.csv").open("w", newline="") as handle:
        csv.DictWriter(handle, fieldnames=[column["key"] for column in payload["columns"]]).writeheader()
    write_csv(source / "batch_summary.csv", cases)
    return source


def test_all_rejected_cohort_preserves_every_failed_candidate(tmp_path):
    reference = cohort(tmp_path, "reference", [("good.pdb", -1)])
    rejected = all_rejected_cohort(tmp_path, reference)
    path = combine([reference, rejected], tmp_path / "combined", 3, include_failures=True)
    result = read(path.parent / "combined_manifest.json")
    assert result["status"] == "complete_with_failures"
    assert result["candidate_count"] == 3
    assert result["completed_count"] == 1 and result["failure_count"] == 2
    payload = read(path.with_suffix(".json"))
    assert [row["candidate"] for row in payload["rows"]] == ["good.pdb", "bad-a.pdb", "bad-b.pdb"]
    assert all(row["status"] == "failed" and row["rank"] is None for row in payload["rows"][1:])
    assert all("rosetta_interface_REU" not in row for row in payload["rows"][1:])
    assert result["sources"][1]["preparation_only"] is True


@pytest.mark.parametrize("mismatch", ["parent", "params", "seed", "inventory"])
def test_all_rejected_cohort_cannot_hide_incompatible_inputs(tmp_path, mismatch):
    reference = cohort(tmp_path, "reference", [("good.pdb", -1)])
    rejected = all_rejected_cohort(tmp_path, reference)
    manifest_path = rejected / "batch_manifest.json"
    manifest = read(manifest_path)
    if mismatch == "params":
        other = tmp_path / "different.params"
        other.write_text("different ligand parameters")
        manifest["params"] = str(other)
        write_json(manifest_path, manifest)
    elif mismatch == "parent":
        original_parent = Path(manifest["parent"])
        other_parent = tmp_path / "different_opt" / "hydrogen_optimized_structures" / "parent_m7.pdb"
        other_parent.parent.mkdir(parents=True)
        other_parent.write_text("different parent coordinates")
        record = read(original_parent.parent.parent / "parent" / "result.json")
        record["optimized_pdb"] = str(other_parent)
        record["artifact_sha256"] = {str(other_parent): sha256(other_parent.read_bytes()).hexdigest()}
        other_record = other_parent.parent.parent / "parent" / "result.json"
        other_record.parent.mkdir()
        write_json(other_record, record)
        manifest["parent"] = str(other_parent)
        write_json(manifest_path, manifest)
    else:
        opth_path = Path(manifest["source_opth_manifest"])
        opth = read(opth_path)
        if mismatch == "seed":
            opth["protocol"]["seed"] = 99
        else:
            opth["cases"][1]["candidate"] = "omitted.pdb"
        write_json(opth_path, opth)
    with pytest.raises(ValueError):
        combine([reference, rejected], tmp_path / "combined", 3, include_failures=True)


def generalized_scoring(source, weights="ref2015", legacy=False, *, opth_weights=None,
                        opth_legacy=None, record_opth=True):
    """Saved test evidence for an explicitly selected non-production score."""
    workflow_path = source / "case/workflow_manifest.json"
    workflow = read(workflow_path)
    weightfile = Path(weights)
    weight_hash = sha256(weightfile.read_bytes()).hexdigest() if weightfile.is_file() else "weights-" + weights
    workflow["scoring_protocol"].update(weights=weights, legacy_pre_talaris=legacy,
                                       weights_sha256=weight_hash, pyrosetta_version="test-build")
    write_json(workflow_path, workflow)
    payload_path = source / "batch_summary_grouped_fingerprints.json"
    payload = read(payload_path)
    payload["scoring"] = {"weights": weights, "legacy_pre_talaris": legacy, "weights_sha256": weight_hash}
    if not legacy:
        payload.pop("scoring_profile", None)
    write_json(payload_path, payload)
    parent = Path(workflow["parent"])
    record_path = parent.parent.parent / "parent/result.json"
    if not record_opth:
        record_path.unlink()
        return
    record = read(record_path)
    selected = weights if opth_weights is None else opth_weights
    selected_legacy = legacy if opth_legacy is None else opth_legacy
    selected_file = Path(selected)
    selected_hash = sha256(selected_file.read_bytes()).hexdigest() if selected_file.is_file() else "weights-" + selected
    record["identity"]["protocol"]["source_protocol"].update(
        weights=selected, legacy_pre_talaris=selected_legacy, weights_sha256=selected_hash,
        template_sha256=sha256(Path(workflow["ligand_template_pdb"]).read_bytes()).hexdigest())
    write_json(record_path, record)


def test_ref2015_requires_opt_in_then_combines_matched_scientific_evidence(tmp_path):
    a = cohort(tmp_path, "a", [("a.pdb", -1)])
    b = cohort(tmp_path, "b", [("b.pdb", -2)])
    generalized_scoring(a)
    generalized_scoring(b)
    with pytest.raises(ValueError, match="Expected legacy RosettaLigand"):
        combine([a, b], tmp_path / "strict")
    path = combine([a, b], tmp_path / "general", allow_other_scoring=True)
    context = read(path.parent / "combined_manifest.json")["scientific_protocol"]
    assert context["weights"] == "ref2015"
    assert context["legacy_pre_talaris"] is False
    assert context["opth_weights"] == "ref2015" and context["opth_passes"] == 1


def test_generalized_mode_records_distinct_opth_and_rescoring_weights(tmp_path):
    a = cohort(tmp_path, "a", [("a.pdb", -1)])
    generalized_scoring(a, "ref2015", False, opth_weights="ligand", opth_legacy=True)
    path = combine([a], tmp_path / "combined", allow_other_scoring=True)
    context = read(path.parent / "combined_manifest.json")["scientific_protocol"]
    assert context["weights"] == "ref2015" and context["legacy_pre_talaris"] is False
    assert context["opth_weights"] == "ligand" and context["opth_legacy_pre_talaris"] is True
    assert context["weights_sha256"] != context["opth_weights_sha256"]


def test_custom_legacy_weights_require_unchanged_file_bytes(tmp_path):
    a = cohort(tmp_path, "a", [("a.pdb", -1)])
    weights = tmp_path / "my_ligand_weights.wts"
    weights.write_text("fa_atr 0.8\n")
    generalized_scoring(a, str(weights), True)
    path = combine([a], tmp_path / "combined", allow_other_scoring=True)
    assert read(path.parent / "combined_manifest.json")["scientific_protocol"]["weights_sha256"] == sha256(weights.read_bytes()).hexdigest()
    weights.write_text("fa_atr 0.5\n")
    with pytest.raises(ValueError, match="Custom scoring weights changed"):
        combine([a], tmp_path / "changed", allow_other_scoring=True)


@pytest.mark.parametrize("mismatch", ["weights", "legacy", "weight_hash", "engine_version"])
def test_generalized_mode_rejects_mixed_or_inconsistent_scoring(tmp_path, mismatch):
    a = cohort(tmp_path, "a", [("a.pdb", -1)])
    b = cohort(tmp_path, "b", [("b.pdb", -2)])
    generalized_scoring(a)
    generalized_scoring(b)
    path = b / "case/workflow_manifest.json"
    workflow = read(path)
    if mismatch == "weights":
        workflow["scoring_protocol"]["weights"] = "beta_nov16"
    elif mismatch == "legacy":
        workflow["scoring_protocol"]["legacy_pre_talaris"] = True
    elif mismatch == "weight_hash":
        workflow["scoring_protocol"]["weights_sha256"] = "different-weights"
    else:
        workflow["scoring_protocol"]["pyrosetta_version"] = "different-build"
    write_json(path, workflow)
    with pytest.raises(ValueError):
        combine([a, b], tmp_path / "combined", allow_other_scoring=True)


def test_generalized_mode_does_not_claim_an_unrecorded_opth_run(tmp_path):
    a = cohort(tmp_path, "a", [("a.pdb", -1)])
    generalized_scoring(a, record_opth=False)
    path = combine([a], tmp_path / "combined", allow_other_scoring=True)
    context = read(path.parent / "combined_manifest.json")["scientific_protocol"]
    assert context["opth_provenance"] == "not_recorded"
    assert context["opth_passes"] is None and context["opth_weights"] is None
    assert "parent_optimized_pdb_sha256" not in context
    assert context["parent_pdb_sha256"]


def test_generalized_mode_accepts_worker_nested_weight_identity_without_opth(tmp_path):
    """The full runner records weights_identity, not a flat weights hash."""
    source = cohort(tmp_path, "source", [("candidate.pdb", -1)])
    generalized_scoring(source, "ligand", True, record_opth=False)
    path = source / "case/workflow_manifest.json"
    workflow = read(path)
    protocol = workflow["scoring_protocol"]
    expected_hash = protocol.pop("weights_sha256")
    protocol["weights_identity"] = {"weights": "ligand", "legacy_pre_talaris": True,
                                    "weights_sha256": expected_hash}
    write_json(path, workflow)
    result = combine([source], tmp_path / "combined", allow_other_scoring=True)
    context = read(result.parent / "combined_manifest.json")["scientific_protocol"]
    assert context["weights_sha256"] == expected_hash
    assert context["pyrosetta_version"] == "test-build"
    assert context["opth_provenance"] == "not_recorded"


@pytest.mark.parametrize("mismatch", ["weights", "legacy_pre_talaris", "weights_sha256"])
def test_generalized_mode_rejects_conflicting_nested_weight_identity(tmp_path, mismatch):
    source = cohort(tmp_path, "source", [("candidate.pdb", -1)])
    generalized_scoring(source, record_opth=False)
    path = source / "case/workflow_manifest.json"
    workflow = read(path)
    protocol = workflow["scoring_protocol"]
    identity = {key: protocol[key] for key in ("weights", "legacy_pre_talaris", "weights_sha256")}
    identity[mismatch] = {"weights": "ligand", "legacy_pre_talaris": True,
                          "weights_sha256": "different-weight-file"}[mismatch]
    protocol["weights_identity"] = identity
    write_json(path, workflow)
    with pytest.raises(ValueError, match="nested identity differ"):
        combine([source], tmp_path / "combined", allow_other_scoring=True)
    assert not (tmp_path / "combined").exists()


def test_generalized_mode_does_not_invent_missing_nested_weight_provenance(tmp_path):
    source = cohort(tmp_path, "source", [("candidate.pdb", -1)])
    generalized_scoring(source, record_opth=False)
    path = source / "case/workflow_manifest.json"
    workflow = read(path)
    protocol = workflow["scoring_protocol"]
    protocol.pop("weights_sha256")
    protocol["weights_identity"] = {"weights": "ref2015", "legacy_pre_talaris": False,
                                    "weights_sha256": None}
    write_json(path, workflow)
    with pytest.raises(ValueError, match="weight hash and PyRosetta version are required"):
        combine([source], tmp_path / "combined", allow_other_scoring=True)


@pytest.mark.parametrize("different_opth", ["weight_set", "unrecorded"])
def test_generalized_mode_cannot_mix_opth_provenance_across_cohorts(tmp_path, different_opth):
    a = cohort(tmp_path, "a", [("a.pdb", -1)])
    second = tmp_path / "second"
    second.mkdir()
    b = cohort(second, "b", [("b.pdb", -2)])
    generalized_scoring(a)
    if different_opth == "weight_set":
        generalized_scoring(b, opth_weights="ligand", opth_legacy=True)
    else:
        generalized_scoring(b, record_opth=False)
    with pytest.raises(ValueError, match="Scientific protocol/parent/parameters differ"):
        combine([a, b], tmp_path / "combined", allow_other_scoring=True)


def test_generalized_mode_validates_every_completed_candidate(tmp_path):
    a = cohort(tmp_path, "a", [("a.pdb", -1), ("b.pdb", -2)])
    generalized_scoring(a)
    other = a / "different_case"
    other.mkdir()
    workflow = read(a / "case/workflow_manifest.json")
    workflow["scoring_protocol"]["weights_sha256"] = "another-scoring-file"
    write_json(other / "workflow_manifest.json", workflow)
    manifest = read(a / "batch_manifest.json")
    manifest["candidates"][1]["result_folder"] = str(other)
    write_json(a / "batch_manifest.json", manifest)
    with pytest.raises(ValueError, match="Candidate scoring/parent/ligand provenance differs"):
        combine([a], tmp_path / "combined", allow_other_scoring=True)
