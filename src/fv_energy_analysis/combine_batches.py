"""Combine compatible completed cohorts without repeating scientific analysis.

All source rows are retained once, with a newly calculated global interface
rank. Source tables, protocols, parent coordinates and parameters are validated
before a fresh CSV/JSON and their provenance manifest are written.
"""
from __future__ import annotations

import argparse
from copy import deepcopy
import csv
from datetime import datetime, timezone
from hashlib import sha256
import json
import math
from pathlib import Path

from .ds_architecture import VERSION
from .paths import native_path
from .report import write_csv, write_json
from .steric_sandwich import METHOD_VERSION


def digest(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def read_json(path: Path):
    return json.loads(path.read_text(encoding="utf-8-sig"))


def candidate_name(value) -> str:
    name = str(value).replace("\\", "/").rsplit("/", 1)[-1]
    if not name or name in {".", "..", "None"}:
        raise ValueError(f"Invalid candidate identity: {value!r}")
    return name


def _typed(value, column):
    if value is None or value == "":
        return None
    kind = column["kind"]
    if kind == "text":
        return str(value)
    if kind == "boolean":
        if isinstance(value, bool):
            return value
        if str(value).lower() in {"true", "false"}:
            return str(value).lower() == "true"
        raise ValueError(f"Invalid boolean for {column['key']}: {value!r}")
    if kind not in {"integer", "number", "energy", "delta_energy"}:
        raise ValueError(f"Unsupported column kind: {kind}")
    if isinstance(value, bool):
        raise ValueError(f"Boolean in numeric column {column['key']}")
    number = float(value)
    if not math.isfinite(number) or kind == "integer" and not number.is_integer():
        raise ValueError(f"Invalid number for {column['key']}: {value!r}")
    return int(number) if kind == "integer" else number


def _validate_table(source, payload, manifest, include_failures=False, allow_other_scoring=False):
    allowed = {"complete", "complete_with_failures"} if include_failures else {"complete"}
    if (manifest.get("status") not in allowed
            or manifest.get("failure_count", 0) and not include_failures
            or manifest.get("architecture_failure_count", 0)):
        raise ValueError(f"Source must be complete without failures: {source}")
    rows, columns = payload["rows"], payload["columns"]
    empty_failure_cohort = (include_failures and manifest.get("preparation_only") is True
                            and manifest.get("failure_count", 0) > 0)
    if (not rows and not empty_failure_cohort or len(rows) != manifest["completed_count"]
            or len(rows) + int(manifest.get("failure_count", 0)) != manifest["candidate_count"]):
        raise ValueError(f"Source row/count mismatch: {source}")
    keys = [column["key"] for column in columns]
    if len(set(keys)) != len(keys):
        raise ValueError(f"Duplicate source column keys: {source}")
    with (source / "batch_summary_grouped_fingerprints.csv").open(encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle)
        csv_rows = list(reader)
        if (set(reader.fieldnames or []) != set(keys) or len(reader.fieldnames or []) != len(keys)
                or len(csv_rows) != len(rows) or any(None in row for row in csv_rows)):
            raise ValueError(f"Source CSV/JSON schema or row count differs: {source}")
    for index, (row, csv_row) in enumerate(zip(rows, csv_rows), 1):
        for column in columns:
            if _typed(row.get(column["key"]), column) != _typed(csv_row.get(column["key"]), column):
                raise ValueError(f"Source CSV/JSON value differs: {source}, row {index}, {column['key']}")
    if allow_other_scoring:
        scoring = payload.get("scoring", {})
        if (not isinstance(scoring.get("weights"), str) or not scoring["weights"].strip()
                or not isinstance(scoring.get("legacy_pre_talaris"), bool)):
            raise ValueError(f"Explicit scoring weights and actual legacy setting are required: {source}")
    elif (payload.get("scoring_profile") != "rosettaligand_legacy"
          or payload.get("scoring", {}).get("weights") != "ligand"
          or payload.get("scoring", {}).get("legacy_pre_talaris") is not True):
        raise ValueError(f"Expected legacy RosettaLigand scores: {source}")
    protocol = payload.get("e2_recognition_protocol") or payload.get("architecture_protocol", {})
    definition = protocol.get("definition_version", protocol.get("version"))
    if (definition != VERSION or protocol.get("steric_method") != METHOD_VERSION
            or protocol.get("interaction_source") != "Discovery Studio"):
        raise ValueError(f"Source has an outdated E2 architecture definition: {source}")


def _failure_rows(source, manifest):
    if not manifest.get("failure_count", 0):
        return []
    with (source / "batch_summary.csv").open(encoding="utf-8-sig", newline="") as handle:
        failed = [row for row in csv.DictReader(handle) if row.get("status") != "complete"]
    if len(failed) != int(manifest["failure_count"]):
        raise ValueError(f"Source failure count differs from detailed rows: {source}")
    names = [candidate_name(row["candidate"]) for row in failed]
    declared = [candidate_name(case["candidate"]) for case in manifest.get("candidates", [])
                if case.get("status") != "complete"]
    if len(set(names)) != len(names) or set(names) != set(declared):
        raise ValueError(f"Source failure identities differ from manifest: {source}")
    return [{"candidate": candidate_name(row["candidate"]), "status": "failed",
             "error": row.get("error") or "Scientific analysis did not complete"} for row in failed]


def _parent_context(parent: Path, params: Path, template: Path, allow_other_scoring=False) -> dict:
    parent_hash = digest(parent)
    params_hash = digest(params)
    record_path = parent.parent.parent / "parent" / "result.json"
    if allow_other_scoring and not record_path.is_file():
        return {"weights": None, "legacy_pre_talaris": None, "weights_sha256": None,
                "params_sha256": params_hash, "template_sha256": digest(template),
                "parent_pdb_sha256": parent_hash, "opth_passes": None, "opth_seed": None,
                "pyrosetta_version": None, "opth_provenance": "not_recorded",
                "opth_weights": None, "opth_legacy_pre_talaris": None,
                "opth_weights_sha256": None, "opth_pyrosetta_version": None}
    record = read_json(record_path)
    opt = record["identity"]["protocol"]
    opt_score = opt["source_protocol"]
    if (record.get("status") != "complete"
            or record.get("artifact_sha256", {}).get(record["optimized_pdb"]) != parent_hash
            or digest(native_path(record["optimized_pdb"])) != parent_hash
            or opt.get("passes_per_complex") != 1
            or (not allow_other_scoring and opt_score.get("weights") != "ligand")
            or (not allow_other_scoring and opt_score.get("legacy_pre_talaris") is not True)
            or opt_score.get("params_sha256") != params_hash):
        raise ValueError(f"Parent OptH protocol/artifact mismatch: {record_path}")
    if allow_other_scoring and (
            not isinstance(opt_score.get("weights"), str) or not opt_score["weights"].strip()
            or not isinstance(opt_score.get("legacy_pre_talaris"), bool)
            or not opt_score.get("weights_sha256") or not opt_score.get("pyrosetta_version")
            or opt_score.get("template_sha256") != digest(template)):
        raise ValueError(f"Parent OptH weight/template provenance is incomplete: {record_path}")
    context = {"weights": opt_score["weights"], "legacy_pre_talaris": opt_score["legacy_pre_talaris"],
            "params_sha256": params_hash, "template_sha256": digest(template),
            "parent_optimized_pdb_sha256": parent_hash, "opth_passes": 1,
            "opth_seed": opt["seed"], "pyrosetta_version": opt_score["pyrosetta_version"],
            "weights_sha256": opt_score["weights_sha256"]}
    if allow_other_scoring:
        context.update(parent_pdb_sha256=parent_hash, opth_provenance="recorded_one_pass",
                       opth_weights=opt_score["weights"],
                       opth_legacy_pre_talaris=opt_score["legacy_pre_talaris"],
                       opth_weights_sha256=opt_score["weights_sha256"],
                       opth_pyrosetta_version=opt_score["pyrosetta_version"])
    return context


def _score_context(protocol, preparation, workflow_path):
    """Separate the actual rescore definition from any preceding OptH score."""
    from .pyrosetta_score import resolve_legacy_pre_talaris
    weights = protocol.get("weights")
    if not isinstance(weights, str) or not weights.strip():
        raise ValueError(f"Scoring weights are missing: {workflow_path}")
    identity = protocol.get("weights_identity", {})
    if not isinstance(identity, dict):
        raise ValueError(f"Scoring weights identity is invalid: {workflow_path}")
    if identity.get("weights", weights) != weights:
        raise ValueError(f"Scoring weights and nested identity differ: {workflow_path}")
    declared_legacy = protocol.get("legacy_pre_talaris")
    nested_legacy = identity.get("legacy_pre_talaris")
    legacy = resolve_legacy_pre_talaris(
        weights, declared_legacy if declared_legacy is not None else nested_legacy)
    if nested_legacy is not None and resolve_legacy_pre_talaris(weights, nested_legacy) != legacy:
        raise ValueError(f"Scoring legacy setting and nested identity differ: {workflow_path}")
    declared_hash, nested_hash = protocol.get("weights_sha256"), identity.get("weights_sha256")
    if declared_hash and nested_hash and declared_hash != nested_hash:
        raise ValueError(f"Scoring weight hash and nested identity differ: {workflow_path}")
    same_opt = (weights == preparation.get("opth_weights")
                and legacy == preparation.get("opth_legacy_pre_talaris"))
    weights_hash = nested_hash or declared_hash or (preparation.get("opth_weights_sha256") if same_opt else None)
    version = protocol.get("pyrosetta_version") or preparation.get("opth_pyrosetta_version")
    if not weights_hash or not version:
        raise ValueError(f"Scoring weight hash and PyRosetta version are required: {workflow_path}")
    if (protocol.get("params_sha256") != preparation["params_sha256"]
            or protocol.get("template_sha256", preparation["template_sha256"]) != preparation["template_sha256"]):
        raise ValueError(f"Scoring protocol and ligand parameters/template differ: {workflow_path}")
    # A custom file must still match the bytes declared by the scientific run.
    candidate = native_path(weights)
    if candidate.is_file() or "/" in weights or "\\" in weights:
        if not candidate.is_file() or digest(candidate) != weights_hash:
            raise ValueError(f"Custom scoring weights changed or are unavailable: {workflow_path}")
    return {**preparation, "weights": weights, "legacy_pre_talaris": legacy,
            "weights_sha256": weights_hash, "pyrosetta_version": version}


def _scientific_context(source: Path, manifest: dict, seen=None, allow_other_scoring=False) -> dict:
    """Resolve refreshes, including all-rejected OptH cohorts, to their inputs."""
    seen = set() if seen is None else seen
    if source in seen:
        raise ValueError("Cycle in source-batch provenance")
    seen.add(source)
    while manifest.get("source_batch"):
        source = native_path(manifest["source_batch"]).resolve()
        if source in seen:
            raise ValueError("Cycle in source-batch provenance")
        seen.add(source)
        manifest = read_json(source / "batch_manifest.json")
    parent = native_path(manifest["parent"])
    if manifest.get("preparation_only") is True:
        reference = native_path(manifest["reference_batch"]).resolve()
        reference_manifest = read_json(reference / "batch_manifest.json")
        if reference_manifest.get("status") != "complete" or reference_manifest.get("failure_count", 0):
            raise ValueError("An all-rejected cohort requires a complete scientific reference batch")
        context = _scientific_context(reference, reference_manifest, seen, allow_other_scoring)
        current = _parent_context(parent, native_path(manifest["params"]), native_path(manifest["ligand_template_pdb"]), allow_other_scoring)
        if allow_other_scoring:
            current.update({key: context[key] for key in ("weights", "legacy_pre_talaris", "weights_sha256", "pyrosetta_version")})
        if current != {key: value for key, value in context.items() if key != "discovery_studio_build"}:
            raise ValueError("Rejected cohort parent/parameters differ from its scientific reference")
        opth = read_json(native_path(manifest["source_opth_manifest"]))
        count = int(manifest["candidate_count"])
        names = [candidate_name(case["candidate"]) for case in opth.get("cases", [])]
        declared = [candidate_name(case["candidate"]) for case in manifest["candidates"]]
        if (opth.get("status") != "complete_with_failures" or opth.get("completed_count") != 0
                or opth.get("candidate_count") != count or opth.get("failure_count") != count
                or len(names) != count or len(set(names)) != count or set(names) != set(declared)
                or any(case.get("status") != "failed" for case in opth["cases"])):
            raise ValueError("Rejected cohort does not match the complete OptH failure inventory")
        protocol = opth["protocol"]
        score = protocol["source_protocol"]
        expected_opt = {key: current.get("opth_" + key) if allow_other_scoring and key in
                        {"weights", "legacy_pre_talaris", "weights_sha256", "pyrosetta_version"} else current[key]
                        for key in ("weights", "legacy_pre_talaris", "params_sha256", "template_sha256",
                                    "weights_sha256", "pyrosetta_version")}
        if (protocol.get("passes_per_complex") != current["opth_passes"]
                or protocol.get("seed") != current["opth_seed"]
                or any(score.get(key) != expected_opt[key] for key in
                       ("weights", "legacy_pre_talaris", "params_sha256", "template_sha256",
                        "weights_sha256", "pyrosetta_version"))
                or digest(native_path(opth["parent_optimized_pdb"])) != current["parent_optimized_pdb_sha256"]):
            raise ValueError("Rejected cohort OptH protocol/parent differs from its reference")
        return context
    first = next((case for case in manifest["candidates"] if case["status"] == "complete"), None)
    if first is None:
        raise ValueError(f"Scientific source has no completed candidates: {source}")
    workflow_path = native_path(first["result_folder"]) / "workflow_manifest.json"
    workflow = read_json(workflow_path)
    protocol = workflow.get("scoring_protocol", {})
    params, template = native_path(workflow["params"]), native_path(workflow["ligand_template_pdb"])
    context = _parent_context(parent, params, template, allow_other_scoring)
    if allow_other_scoring:
        context = _score_context(protocol, context, workflow_path)
        hashes = {parent.resolve(): context["parent_pdb_sha256"],
                  params.resolve(): context["params_sha256"], template.resolve(): context["template_sha256"]}
        def artifact_hash(value):
            path = native_path(value).resolve()
            if path not in hashes:
                hashes[path] = digest(path)
            return hashes[path]
        for case in manifest["candidates"]:
            if case["status"] != "complete":
                continue
            path = native_path(case["result_folder"]) / "workflow_manifest.json"
            current = read_json(path)
            if (artifact_hash(current["parent"]) != context["parent_pdb_sha256"]
                    or artifact_hash(current["params"]) != context["params_sha256"]
                    or artifact_hash(current["ligand_template_pdb"]) != context["template_sha256"]
                    or _score_context(current.get("scoring_protocol", {}), context, path) != context):
                raise ValueError(f"Candidate scoring/parent/ligand provenance differs: {path}")
    elif protocol.get("weights") != "ligand" or protocol.get("params_sha256") != context["params_sha256"]:
        raise ValueError(f"Scoring protocol and ligand parameters differ: {workflow_path}")
    if digest(native_path(workflow["parent"])) != context.get("parent_pdb_sha256", context.get("parent_optimized_pdb_sha256")):
        raise ValueError(f"Workflow parent differs from batch parent: {workflow_path}")
    ds = manifest.get("discovery_studio", {})
    if ds.get("status") != "complete" or ds.get("add_hydrogens") is not False:
        raise ValueError(f"Completed DS with AddHydrogens=False required: {source}")
    return {**context, "discovery_studio_build": ds.get("build_info")}


def _expected_names(path: Path | None) -> set[str] | None:
    if path is None:
        return None
    values = read_json(path)
    if isinstance(values, dict):
        values = values.get("expected_candidates", values.get("candidates"))
    if not isinstance(values, list) or not all(isinstance(value, str) for value in values):
        raise ValueError("Expected-candidates JSON must contain an array of candidate strings")
    names = [candidate_name(value) for value in values]
    if len(set(name.casefold() for name in names)) != len(names):
        raise ValueError("Duplicate identifiers in expected-candidates list")
    return set(names)


def combine(sources: list[Path], output: Path, expected_count=None, expected_candidates=None,
            include_failures=False, allow_other_scoring=False) -> Path:
    if not sources:
        raise ValueError("At least one source is required")
    sources = [native_path(source).resolve() for source in sources]
    output = native_path(output).resolve()
    if len(set(sources)) != len(sources):
        raise ValueError("Duplicate source folders")
    if any(output == source or source in output.parents for source in sources):
        raise ValueError("Use a separate combined output folder outside all sources")
    if output.exists() and any(output.iterdir()):
        raise FileExistsError(f"Combined output must be fresh: {output}")
    expected = _expected_names(native_path(expected_candidates) if expected_candidates else None)
    rows, failed_rows, origins, source_records = [], [], {}, []
    baseline, schema, context, parent = None, None, None, None
    for source in sources:
        payload_path, manifest_path = source / "batch_summary_grouped_fingerprints.json", source / "batch_manifest.json"
        payload, manifest = read_json(payload_path), read_json(manifest_path)
        _validate_table(source, payload, manifest, include_failures, allow_other_scoring)
        current_schema = [tuple(column[key] for key in ("key", "label", "section", "kind")) for column in payload["columns"]]
        current_parent = {key: value for key, value in payload["parent"].items() if key not in {"rank", "candidate"}}
        current_context = _scientific_context(source, manifest, allow_other_scoring=allow_other_scoring)
        if allow_other_scoring:
            declared = payload["scoring"]
            if (declared["weights"] != current_context["weights"]
                    or declared["legacy_pre_talaris"] is not current_context["legacy_pre_talaris"]
                    or declared.get("weights_sha256", current_context["weights_sha256"]) != current_context["weights_sha256"]):
                raise ValueError(f"Report scoring declaration differs from scientific artifacts: {source}")
        if baseline is None:
            baseline, schema, context, parent = payload, current_schema, current_context, current_parent
        elif current_schema != schema:
            raise ValueError(f"Source column schema differs: {source}")
        elif current_parent != parent:
            differing = sorted(key for key in set(parent) | set(current_parent) if parent.get(key) != current_parent.get(key))
            raise ValueError(f"Parent baseline differs in {differing}: {source}")
        elif current_context != context:
            raise ValueError(f"Scientific protocol/parent/parameters differ: {source}")
        failures = _failure_rows(source, manifest) if include_failures else []
        for original, failed in [*((row, False) for row in payload["rows"]),
                                 *((row, True) for row in failures)]:
            row = dict(original)
            name = candidate_name(row["candidate"])
            identity = name.casefold()
            if identity in origins:
                raise ValueError(f"Duplicate candidate {name}: {origins[identity]['source']} and {source}")
            row["candidate"] = name
            origins[identity] = {"candidate": name, "source": str(source), "source_rank": row.get("rank")}
            if failed:
                failed_rows.append(row)
            else:
                if include_failures:
                    row.update(status="complete", error=None)
                rows.append(row)
        source_records.append({"source": str(source), "candidate_count": len(payload["rows"]) + len(failures),
                               "completed_count": len(payload["rows"]), "failure_count": len(failures),
                               "preparation_only": manifest.get("preparation_only", False),
                               "manifest_sha256": digest(manifest_path), "payload_sha256": digest(payload_path),
                               "csv_sha256": digest(source / "batch_summary_grouped_fingerprints.csv")})
    names = {row["candidate"] for row in [*rows, *failed_rows]}
    total = len(rows) + len(failed_rows)
    if expected_count is not None and total != expected_count:
        raise ValueError(f"Expected {expected_count} candidates, found {total}")
    if expected is not None and names != expected:
        raise ValueError(f"Combined candidate inventory differs; missing={sorted(expected-names)}, unexpected={sorted(names-expected)}")
    for row in rows:
        for key in ("rosetta_interface_delta_REU", "rosetta_interface_REU"):
            if isinstance(row.get(key), bool) or not math.isfinite(float(row.get(key))):
                raise ValueError(f"Missing/nonfinite interface energy for {row['candidate']}: {key}")
    rows.sort(key=lambda row: (float(row["rosetta_interface_delta_REU"]), row["candidate"]))
    for rank, row in enumerate(rows, 1):
        row["rank"] = rank
    completed_count = len(rows)
    for row in failed_rows:
        row["rank"] = None
    rows.extend(sorted(failed_rows, key=lambda row: row["candidate"]))
    combined = deepcopy(baseline)
    combined.pop("architecture_protocol", None)
    combined.update(title="Combined Fv energy analysis" if allow_other_scoring else "Combined m7 RosettaLigand analysis", rows=rows,
                    summary_subtitle=f"{completed_count} completed, {len(failed_rows)} failed; global rank by mutant-minus-parent interface energy.",
                    note="One row per unique candidate. All source metrics are retained; ranks were recalculated across cohorts.",
                    e2_recognition_protocol={"definition_version": VERSION, "steric_method": METHOD_VERSION,
                                             "interaction_source": "Discovery Studio", "maximum_surface_gap_A": 1.0,
                                             "c19_obstruction_included": False},
                    combined_sources=source_records)
    if include_failures:
        combined["columns"].extend([
            {"key": "status", "label": "Analysis status", "section": "Analysis status", "kind": "text", "width": 22},
            {"key": "error", "label": "Failure reason", "section": "Analysis status", "kind": "text", "width": 80},
        ])
        combined["parent"].update(status="reference", error=None)
    # Score mappings are source provenance, not a partial mapping inherited
    # from whichever cohort happened to be read first.
    combined.pop("model_scores", None)
    combined.pop("hcdr3_confidence", None)
    output.mkdir(parents=True, exist_ok=True)
    csv_path = output / "batch_summary_grouped_fingerprints.csv"
    keys = [column["key"] for column in combined["columns"]]
    write_csv(csv_path, [{key: row.get(key) for key in keys} for row in rows])
    write_json(output / "batch_summary_grouped_fingerprints.json", combined)
    write_csv(output / "candidate_sources.csv", [origins[row["candidate"].casefold()] for row in rows])
    write_json(output / "combined_manifest.json", {
        "status": "complete_with_failures" if failed_rows else "complete", "completed_utc": datetime.now(timezone.utc).isoformat(),
        "candidate_count": len(rows), "completed_count": completed_count, "failure_count": len(failed_rows),
        "expected_count": expected_count, "expected_candidates_sha256": digest(native_path(expected_candidates)) if expected_candidates else None,
        "columns": len(keys), "rank_metric": "rosetta_interface_delta_REU ascending; candidate identity breaks ties",
        "scientific_protocol": context, "sources": source_records, "scoring_rerun": False, "detectors_rerun": [],
        "metrics_changed": ["rank"], "include_failures": include_failures,
        **({"allow_other_scoring": True} if allow_other_scoring else {}),
        "csv": str(csv_path), "csv_sha256": digest(csv_path),
    })
    return csv_path


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, action="append", required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--expected-count", type=int)
    parser.add_argument("--expected-candidates", type=Path)
    parser.add_argument("--include-failures", action="store_true",
                        help="Keep failed candidates with explicit status/reason and blank scientific metrics")
    parser.add_argument("--allow-other-scoring", action="store_true",
                        help="Allow other matched scoring/OptH protocols; verify their actual weight hashes and legacy settings")
    args = parser.parse_args(argv)
    print(combine(args.source, args.output, args.expected_count, args.expected_candidates,
                  args.include_failures, args.allow_other_scoring), flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
