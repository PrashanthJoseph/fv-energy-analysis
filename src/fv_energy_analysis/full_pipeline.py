"""Run the packaged folder pipeline; child stages block and reuse saved results.

Scientific stages run in the active Linux environment or an explicitly selected
WSL environment. Native Discovery Studio runs in Windows; Linux callers must
explicitly skip it, and their status records that those features are unavailable.
"""

from __future__ import annotations

import argparse
import csv
from hashlib import sha256
import json
import os
from pathlib import Path, PureWindowsPath
import re
import shlex
import subprocess
import sys
import traceback
from datetime import datetime, timezone


PATH_KINDS = {"input_folder": "dir", "parent": "file", "params": "file", "ligand_template_pdb": "file"}
OPTIONAL_PATH_KINDS = {
    "workspace": "dir", "science_python": "file",
    "mutant": "file",
    "reuse_mutants_from": "dir", "ligand_source": "file",
    "cif_folder": "dir", "parent_cif": "file", "scores_csv": "file",
    "hydrogen_parent": "file", "reuse_parent_result": "file",
    "c18_calibration": "file", "combined_expected_candidates": "file",
    "failure_reference_batch": "dir",
    "reuse_parent_cache": "dir",
}


def utc_now():
    return datetime.now(timezone.utc).isoformat()


def configure_workspace_imports(config):
    """Support source checkouts and installed wheels without temporary helpers."""
    source = str(Path(__file__).resolve().parents[1])
    if source not in sys.path:
        sys.path.insert(0, source)


def wsl_path(value):
    value = str(value)
    if value.startswith("/"):
        return value
    path = PureWindowsPath(value)
    if not re.fullmatch(r"[A-Za-z]:", path.drive):
        raise ValueError(f"An absolute Windows drive path is required: {value}")
    return f"/mnt/{path.drive[0].lower()}/" + "/".join(path.parts[1:])


def read_config(config_path):
    config_path = Path(config_path).resolve()
    from .pipeline_options import normalize_options
    config = normalize_options(json.loads(config_path.read_text(encoding="utf-8-sig")))
    def resolved(value):
        path = Path(value).expanduser()
        return (path if path.is_absolute() else config_path.parent / path).resolve()
    if config.get("mutant"):
        if config.get("input_folder"):
            raise ValueError("Supply mutant for one complex or input_folder for a batch, not both")
        selected = resolved(config["mutant"])
        config.update(mutant=str(selected), input_folder=str(selected.parent),
                      input_glob=selected.name, recursive=False)
    for key, kind in PATH_KINDS.items():
        if not config.get(key):
            raise ValueError(f"Missing required configuration field: {key}")
        path = resolved(config[key])
        if not (path.is_dir() if kind == "dir" else path.is_file()):
            raise FileNotFoundError(f"Configuration {key} is not an existing {kind}: {path}")
        config[key] = str(path)
    for key, kind in OPTIONAL_PATH_KINDS.items():
        if config.get(key):
            path = resolved(config[key])
            if not (path.is_dir() if kind == "dir" else path.is_file()):
                raise FileNotFoundError(f"Configuration {key} is not an existing {kind}: {path}")
            config[key] = str(path)
    for key in ("output",):
        if not config.get(key):
            raise ValueError(f"Missing required configuration field: {key}")
    config["output"] = str(resolved(config["output"]))
    config.setdefault("workspace", str(config_path.parent))
    backend = config.setdefault("science_backend", "wsl" if config.get("wsl_python") else "native")
    if backend not in {"native", "wsl"}:
        raise ValueError("science_backend must be native or wsl")
    if backend == "wsl":
        if not config.get("distro") or not str(config.get("wsl_python", "")).startswith("/"):
            raise ValueError("WSL requires distro and an absolute Linux wsl_python interpreter path")
    elif os.name == "nt" and not config.get("postprocessing_only"):
        raise ValueError("Scientific stages require Linux/PyRosetta; select science_backend=wsl on Windows")
    if not config.get("skip_discovery_studio") and os.name != "nt":
        raise ValueError("Native Discovery Studio requires a Windows controller; use skip_discovery_studio=true for an explicitly partial Linux run")
    if config.get("skip_discovery_studio") and config.get("e3_discrimination"):
        raise ValueError("The full E3 report requires completed DS architecture; run E3 separately for a DS-free diagnostic")
    config["jobs"] = int(config.get("jobs", 2))
    if config["jobs"] < 1:
        raise ValueError("jobs must be a positive integer")
    config["weights"] = str(config.get("weights", "ligand")).strip()
    if not config["weights"]:
        raise ValueError("weights must name a Rosetta scoring function")
    from .pyrosetta_score import scoring_weights_identity
    for key, override in (("weights", "legacy_pre_talaris"), ("hydrogen_weights", "hydrogen_legacy_pre_talaris")):
        value = str(config.get(key, config["weights"]))
        candidate = resolved(value)
        if candidate.is_file():
            value = str(candidate)
            if config.get(override, config.get("legacy_pre_talaris")) is None:
                raise ValueError(f"Set {override}=true or false explicitly for custom Rosetta weights")
        identity = scoring_weights_identity(value, config.get(override, config.get("legacy_pre_talaris")))
        config[key] = str(identity["weights"])
        config[key + "_identity"] = identity
    config.setdefault("include_empty_metadata_sections", True)
    if config.get("hydrogen_output") or (config.get("optimize_hydrogens", True) and not config.get("postprocessing_only")):
        config.setdefault("hydrogen_parent", config["parent"])
        config["hydrogen_output"] = str(resolved(config.get("hydrogen_output", str(Path(config["output"]).with_name(Path(config["output"]).name + "_opth")))))
        matches = Path(config["input_folder"]).rglob(config.get("input_glob", "*.pdb")) if config.get("recursive") else Path(config["input_folder"]).glob(config.get("input_glob", "*.pdb"))
        if Path(config["hydrogen_parent"]).suffix.lower() in {".cif", ".mmcif"} or any(p.is_file() and p.suffix.lower() in {".cif", ".mmcif"} for p in matches):
            raise ValueError("The OptH folder stage currently accepts PDB complexes. Convert mmCIF first, or disable optimize_hydrogens to use the scientific workflow's mmCIF conversion.")
    if config.get("discovery_studio_root"):
        config["discovery_studio_root"] = str(resolved(config["discovery_studio_root"]))
    if not config.get("skip_discovery_studio") and not config.get("reuse_discovery_studio"):
        from .discovery_studio import DEFAULT_DS_ROOT
        ds_root = Path(config.get("discovery_studio_root", DEFAULT_DS_ROOT))
        if not (ds_root / "bin/perl.bat").is_file():
            raise FileNotFoundError(f"Discovery Studio is unavailable: {ds_root / 'bin/perl.bat'}")
    if config.get("parent_cif") and not config.get("cif_folder"):
        raise ValueError("parent_cif requires cif_folder")
    if config.get("reuse_mutants_from") and Path(config["output"]) == Path(config["reuse_mutants_from"]):
        raise ValueError("Output must be separate from the saved batch used for reuse")
    if config.get("fresh_start"):
        if config.get("reuse_mutants_from"):
            raise ValueError("A fresh start must not reuse a saved mutant batch")
        if Path(config["output"]).exists():
            raise FileExistsError("Fresh output already exists; refusing to mix runs")
    for key in ("e3_output", "e3_source_output", "combined_output", "expected_xlsx", "report_output"):
        if config.get(key):
            config[key] = str(resolved(config[key]))
    if config.get("combine_sources"):
        config["combine_sources"] = [str(resolved(value)) for value in config["combine_sources"]]
        if config.get("e3_discrimination"):
            config["combined_include_failures"] = True
    for key in ("params", "ligand_template_pdb", "ligand_source"):
        expected = config.get(f"{key}_sha256")
        if expected and sha256(Path(config[key]).read_bytes()).hexdigest() != expected.lower():
            raise ValueError(f"Validated scientific input has changed: {key}")
    from .report_options import validate_report_selection
    validate_report_selection(config)
    return config


def science_backend(config):
    return config.get("science_backend", "wsl" if config.get("wsl_python") else "native")


def science_path(config, value):
    return wsl_path(value) if science_backend(config) == "wsl" else str(value)


def science_python(config):
    return config["wsl_python"] if science_backend(config) == "wsl" else config.get("science_python", sys.executable)


def scientific_command(config, arguments):
    return wsl_command(config, arguments) if science_backend(config) == "wsl" else arguments


def wsl_command(config, arguments):
    workspace = wsl_path(config["workspace"])
    python = config["wsl_python"]
    source = workspace + "/src"
    source_export = f"PYTHONPATH={shlex.quote(source)} " if (Path(config["workspace"]) / "src").is_dir() else ""
    command = (
        f"cd {shlex.quote(workspace)} && "
        f"export {source_export}"
        f"PATH={shlex.quote(str(Path(python).parent).replace(chr(92), '/') + ':/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin')} "
        "OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 && "
        f"exec {shlex.join(arguments)}"
    )
    return ["wsl.exe", "-d", config["distro"], "--", "/bin/bash", "-c", command]


def hydrogen_command(config):
    arguments = [science_python(config), "-u", "-m", "fv_energy_analysis.hydrogen_batch"]
    for flag, key in (("--parent", "hydrogen_parent"), ("--input-folder", "input_folder"),
                      ("--output", "hydrogen_output"), ("--params", "params"),
                      ("--template", "ligand_template_pdb"), ("--reuse-parent-result", "reuse_parent_result")):
        if config.get(key):
            arguments.extend([flag, science_path(config, config[key])])
    arguments.extend(["--jobs", str(config["jobs"]), "--seed", str(config.get("hydrogen_seed", 20261022)),
                      "--weights", science_weights(config, config.get("hydrogen_weights", config.get("weights", "ligand"))),
                      "--ligand-chain", config.get("hydrogen_ligand_chain", config.get("mutant_ligand_chain", config.get("ligand_chain", "d"))),
                      "--parent-ligand-chain", config.get("hydrogen_parent_ligand_chain", config.get("parent_ligand_chain", config.get("ligand_chain", "d"))),
                      "--glob", config.get("input_glob", "*.pdb")])
    legacy = config.get("hydrogen_legacy_pre_talaris", config.get("legacy_pre_talaris"))
    if legacy is not None:
        arguments.append("--legacy-pre-talaris" if legacy else "--no-legacy-pre-talaris")
    if config.get("recursive"):
        arguments.append("--recursive")
    if config.get("expected_count") is not None:
        arguments.extend(["--expected-count", str(config["expected_count"])])
    return scientific_command(config, arguments)


def batch_command(config):
    python = science_python(config)
    arguments = [python, "-u", "-m", "fv_energy_analysis.batch_workflow"]
    for flag, key in (("--parent", "parent"), ("--input-folder", "input_folder"),
                      ("--output", "output"),
                      ("--params", "params"), ("--ligand-template-pdb", "ligand_template_pdb")):
        arguments.extend([flag, science_path(config, config[key])])
    for flag, key in (("--reuse-mutants-from", "reuse_mutants_from"), ("--ligand-source", "ligand_source"),
                      ("--reuse-parent-cache", "reuse_parent_cache")):
        if config.get(key):
            arguments.extend([flag, science_path(config, config[key])])
    for flag, value in (
        ("--mutant-chain-map", config.get("mutant_chain_map")),
        ("--parent-chain-map", config.get("parent_chain_map")),
        ("--ligand-chain", config.get("ligand_chain")),
        ("--parent-ligand-chain", config.get("parent_ligand_chain")),
        ("--mutant-ligand-chain", config.get("mutant_ligand_chain")),
        ("--rosetta-ligand-chain", config.get("rosetta_ligand_chain")),
    ):
        if value:
            arguments.extend([flag, str(value)])
    arguments.extend(["--jobs", str(config["jobs"]), "--glob", config.get("input_glob", "*.pdb"),
                      "--weights", science_weights(config, config.get("weights", "ligand")), "--arpeggio-scope", config.get("arpeggio_scope", "neighbourhood"),
                      "--arpeggio-neighbourhood-radius", str(config.get("arpeggio_neighbourhood_radius", 6.0))])
    if config.get("recursive"):
        arguments.append("--recursive")
    for key, flag in (("scheme", "--scheme"), ("cdr_definition", "--cdr-definition"),
                      ("fit_region", "--fit-region"), ("rmsd_atoms", "--rmsd-atoms"),
                      ("rotamer_rmsd_cutoff", "--rotamer-rmsd-cutoff"),
                      ("hcdr3_tip_sequence", "--hcdr3-tip-sequence"),
                      ("network_energy_presence_cutoff", "--network-energy-presence-cutoff"),
                      ("network_energy_change_cutoff", "--network-energy-change-cutoff"),
                      ("fv_chains", "--fv-chains"), ("ligand_resname", "--ligand-resname"),
                      ("ligand_smiles", "--ligand-smiles")):
        if key in config:
            value = config[key]
            arguments.extend([flag, ",".join(value) if isinstance(value, list) else str(value)])
    legacy = config.get("legacy_pre_talaris")
    if legacy is None:
        legacy = config.get("weights_identity", {}).get("legacy_pre_talaris")
    if legacy is not None:
        arguments.append("--legacy-pre-talaris" if legacy else "--no-legacy-pre-talaris")
    for key, flag in (("arpeggio", "--skip-arpeggio"), ("plip", "--skip-plip"),
                      ("network", "--skip-network-analysis"), ("rotamers", "--skip-rotamer-analysis"),
                      ("e2_architecture", "--skip-e2-architecture")):
        if not config.get("analyses", {}).get(key, True):
            arguments.append(flag)
    if config.get("stop_on_error"):
        arguments.append("--stop-on-error")
    return scientific_command(config, arguments)


def science_weights(config, value):
    """Translate a custom file path, retaining built-in Rosetta names as names."""
    return science_path(config, value) if Path(value).is_file() else str(value)


def fingerprint_command(config):
    """Append optional model metadata only when it belongs to this dataset."""
    arguments = [sys.executable, "-u", "-m", "fv_energy_analysis.fingerprints",
                 "--batch-results", config["output"]]
    for flag, key in (("--cif-folder", "cif_folder"), ("--parent-cif", "parent_cif"),
                      ("--scores-csv", "scores_csv")):
        if config.get(key):
            arguments.extend([flag, config[key]])
    if config.get("include_empty_metadata_sections"):
        arguments.append("--include-empty-metadata-sections")
    for key, flag in (("discovery_studio", "--skip-discovery-studio"),
                      ("e2_architecture", "--skip-e2-architecture"),
                      ("hcdr3_plddt", "--skip-hcdr3-plddt"),
                      ("boltz2_scores", "--skip-boltz2-scores"),
                      ("gnina_scores", "--skip-gnina-scores")):
        if not config.get("analyses", {}).get(key, True):
            arguments.append(flag)
    return arguments


def optimized_analysis_config(config):
    """OptH normalizes ligand names/chain, but retains protein chain labels."""
    folder = Path(config["hydrogen_output"]) / "hydrogen_optimized_structures"
    return {**config, "input_folder": str(folder), "parent": str(folder / "parent_m7.pdb"),
            "ligand_chain": "X", "parent_ligand_chain": "X", "mutant_ligand_chain": "X",
            "rosetta_ligand_chain": "X"}


def preflight_commands(config):
    """Check dependencies in the interpreters that will actually run each stage."""
    if not config.get("postprocessing_only"):
        science = [science_python(config), "-u", "-m", "fv_energy_analysis.doctor", "--profile", "analysis"]
        for tool in ("arpeggio", "plip"):
            if not config.get("analyses", {}).get(tool, True):
                science.extend(["--skip-tool", tool])
        yield "science_dependencies", scientific_command(config, science)
    arguments = [sys.executable, "-u", "-m", "fv_energy_analysis.doctor", "--profile", "reporting"]
    if not config.get("skip_discovery_studio") and config.get("discovery_studio_root"):
        arguments.extend(["--discovery-studio-root", config["discovery_studio_root"]])
    yield "reporting_dependencies", arguments


def merge_preparation_failures(config, hydrogen):
    """Keep rejected raw inputs in the final master counts and failure sheet."""
    failures = [case for case in hydrogen["cases"] if case["status"] != "complete"]
    if not failures:
        return
    from fv_energy_analysis.report import write_csv, write_json
    from fv_energy_analysis.paths import native_path
    output = Path(config["output"])
    manifest_path = output / "batch_manifest.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    summary_path = output / "batch_summary.csv"
    with summary_path.open(encoding="utf-8-sig", newline="") as stream:
        rows = list(csv.DictReader(stream))
    names = {Path(row["candidate"]).name for row in rows}
    expected_names = {case["candidate"] for case in failures}
    already_present = names & expected_names
    if already_present:
        existing = {Path(row["candidate"]).name: row for row in rows}
        if (already_present != expected_names
                or manifest.get("preparation_failure_count") != len(failures)
                or any(existing[name]["status"] != "failed" for name in expected_names)):
            raise ValueError("Preparation-failure inventory conflicts with the completed analysis")
        for case in failures:
            expected_error = "OptH preparation: " + case["error_type"] + ": " + case["error"]
            if existing[case["candidate"]].get("error") != expected_error:
                raise ValueError("Recorded preparation failure differs from the saved OptH diagnostic")
        return
    for case in failures:
        if case["candidate"] in names:
            raise ValueError("Candidate present in both OptH failures and full analysis")
        row = {"candidate": case["candidate"], "input_path": str(native_path(hydrogen["input_folder"]) / case["candidate"]),
               "result_folder": case["diagnostic_folder"], "status": "failed",
               "error": "OptH preparation: " + case["error_type"] + ": " + case["error"]}
        rows.append(row)
        manifest["candidates"].append({**row, "stage": "hydrogen_optimization"})
    manifest.update(candidate_count=hydrogen["candidate_count"],
                    failure_count=int(manifest.get("failure_count", 0)) + len(failures),
                    preparation_failure_count=len(failures), status="complete_with_failures")
    write_csv(summary_path, rows)
    write_json(output / "preparation_failures.json", failures)
    write_json(manifest_path, manifest)
    with (output / "batch_summary.md").open("a", encoding="utf-8") as stream:
        stream.write(f"\n\nOptH preparation rejected {len(failures)} additional raw inputs. "
                     "They remain in the master CSV and workbook failure sheet with per-case diagnostics.\n")


def finalization_commands(config):
    """Refresh cheap saved geometry, then merge disjoint compatible cohorts."""
    refresh = config.get("architecture_refresh")
    if refresh:
        arguments = [sys.executable, "-u", "-m", "fv_energy_analysis.architecture_batch"]
        for flag, value in (("--source", refresh["source"]), ("--output", refresh["output"]),
                            ("--opth-results", refresh["opth_results"]), ("--params", config["params"]),
                            ("--template", config["ligand_template_pdb"]), ("--c18-calibration", config["c18_calibration"])):
            arguments.extend([flag, str(value)])
        arguments.extend(["--jobs", str(config["jobs"]), "--expected-count", str(refresh["expected_count"])])
        if refresh.get("resume"):
            arguments.append("--resume")
        yield "refresh_previous_architecture", arguments
    if config.get("combine_sources"):
        arguments = [sys.executable, "-u", "-m", "fv_energy_analysis.combine_batches",
                     "--output", config["combined_output"], "--expected-count", str(config["combined_expected_count"])]
        for source in config["combine_sources"]:
            arguments.extend(["--source", source])
        if config.get("combined_expected_candidates"):
            arguments.extend(["--expected-candidates", config["combined_expected_candidates"]])
        if config.get("combined_include_failures"):
            arguments.append("--include-failures")
        if needs_general_scoring(config):
            arguments.append("--allow-other-scoring")
        yield "combine_full_csv", arguments
    if config.get("e3_discrimination"):
        source = e3_source_folder(config)
        if not config.get("combine_sources"):
            source_arguments = [sys.executable, "-u", "-m", "fv_energy_analysis.combine_batches",
                "--source", config["output"], "--output", str(source), "--include-failures"]
            if needs_general_scoring(config):
                source_arguments.append("--allow-other-scoring")
            yield "e3_source_inventory", source_arguments
        arguments = [science_python(config), "-u", "-m", "fv_energy_analysis.e3_discrimination",
                     "--source", science_path(config, source), "--output", science_path(config, final_report_folder(config)),
                     "--params", science_path(config, config["params"]),
                     "--template", science_path(config, config["ligand_template_pdb"]),
                     "--jobs", str(config["jobs"]), "--resume"]
        if "e3_overlap_warning" in config:
            arguments.extend(["--overlap-warning", str(config["e3_overlap_warning"])])
        yield "e3_steric_discrimination", scientific_command(config, arguments)


def e3_source_folder(config):
    output = Path(config["output"])
    return Path(config.get("combined_output", config.get("e3_source_output", output.with_name(output.name + "_combined"))))


def needs_general_scoring(config):
    """Request the merger's explicit generic mode for non-production protocols."""
    return (config.get("weights", "ligand") != "ligand"
            or config.get("legacy_pre_talaris") is False
            or config.get("hydrogen_weights", config.get("weights", "ligand")) != "ligand"
            or config.get("hydrogen_legacy_pre_talaris") is False
            or config.get("optimize_hydrogens") is False)


def final_report_folder(config):
    output = Path(config["output"])
    if config.get("e3_discrimination"):
        return Path(config.get("e3_output", output.with_name(output.name + "_e3")))
    return Path(config.get("combined_output", output)) if config.get("combine_sources") else output


def discovery_studio_command(config, *, reuse=False):
    arguments = [sys.executable, "-u", "-m", "fv_energy_analysis.discovery_studio",
                 "--batch-results", config["output"], "--skip-fingerprint-reports"]
    if config.get("discovery_studio_root"):
        arguments.extend(["--discovery-studio-root", config["discovery_studio_root"]])
    if reuse or config.get("reuse_discovery_studio"):
        arguments.append("--reports-only")
    return arguments


def discovery_studio_identity(config):
    """Bind reusable DS results to the exact saved PDBs, helper and installation."""
    from .paths import native_path
    from .resources import resource_path
    from .discovery_studio import DEFAULT_DS_ROOT
    batch = Path(config["output"])
    manifest = json.loads((batch / "batch_manifest.json").read_text(encoding="utf-8"))
    coordinates = {}
    for item in manifest.get("candidates", []):
        if item["status"] != "complete":
            continue
        folder = native_path(item["result_folder"])
        summary = json.loads((folder / "final/summary.json").read_text(encoding="utf-8"))
        for role in ("parent", "mutant"):
            pdb = native_path(summary[role])
            if role == "parent" and "parent" in coordinates:
                continue
            key = "parent" if role == "parent" else Path(item["candidate"]).name
            coordinates[key] = sha256(pdb.read_bytes()).hexdigest()
    ds_root = Path(config.get("discovery_studio_root", DEFAULT_DS_ROOT))
    buildinfo = ds_root / "buildinfo.txt"
    return {"pdb_sha256": coordinates, "installation_root": str(ds_root),
            "build_info": sha256(buildinfo.read_bytes()).hexdigest() if buildinfo.is_file() else None,
            "helper_sha256": sha256(resource_path("scripts/discovery_studio_batch_interactions.pl").read_bytes()).hexdigest()}


def combination_identity(command):
    """Validate source tables and bind a report-only merge to their exact metrics."""
    from .combine_batches import _scientific_context, _validate_table
    sources = [Path(command[index + 1]) for index, value in enumerate(command) if value == "--source"]
    records = []
    for source in sources:
        manifest = json.loads((source / "batch_manifest.json").read_text(encoding="utf-8"))
        payload_path = source / "batch_summary_grouped_fingerprints.json"
        payload = json.loads(payload_path.read_text(encoding="utf-8"))
        generic = {"allow_other_scoring": True} if "--allow-other-scoring" in command else {}
        _validate_table(source, payload, manifest, "--include-failures" in command, **generic)
        record = {"source": str(source.resolve()), "payload_sha256": sha256(payload_path.read_bytes()).hexdigest(),
                  "csv_sha256": sha256((source / "batch_summary_grouped_fingerprints.csv").read_bytes()).hexdigest(),
                  "scientific_protocol": _scientific_context(source, manifest, **generic)}
        if manifest.get("failure_count"):
            record["failures_sha256"] = sha256((source / "batch_summary.csv").read_bytes()).hexdigest()
        records.append(record)
    expected_file = command[command.index("--expected-candidates") + 1] if "--expected-candidates" in command else None
    return {"sources": records, "include_failures": "--include-failures" in command,
            "expected_count": command[command.index("--expected-count") + 1] if "--expected-count" in command else None,
            "expected_candidates_sha256": sha256(Path(expected_file).read_bytes()).hexdigest() if expected_file else None,
            "implementation_sha256": sha256(Path(__file__).with_name("combine_batches.py").read_bytes()).hexdigest()}


def combination_artifacts(output):
    """Exclude workbook-only metadata from the immutable merge identity."""
    manifest = json.loads((output / "combined_manifest.json").read_text(encoding="utf-8"))
    if manifest.get("status") not in {"complete", "complete_with_failures"}:
        raise ValueError("Saved combined report is incomplete")
    keys = ("status", "candidate_count", "completed_count", "failure_count", "columns", "scientific_protocol", "include_failures")
    return {"files": {name: sha256((output / name).read_bytes()).hexdigest() for name in
                      ("batch_summary_grouped_fingerprints.json", "batch_summary_grouped_fingerprints.csv", "candidate_sources.csv")},
            "manifest": {key: manifest.get(key) for key in keys}}


def reuse_combination(output, identity):
    if not output.exists() or not any(output.iterdir()):
        return False
    record = output / "full_pipeline_combination_identity.json"
    if not record.is_file():
        raise FileExistsError(f"Existing combined output has no verified controller provenance; choose a fresh output: {output}")
    saved = json.loads(record.read_text(encoding="utf-8"))
    if saved.get("identity") != identity or saved.get("artifacts") != combination_artifacts(output):
        raise ValueError(f"Combined input/protocol or saved report has changed; choose a fresh output: {output}")
    return True


def initialize_rejected_cohort(config, hydrogen):
    """Retain every rejected input when no mutant passes ligand preparation."""
    from copy import deepcopy
    from fv_energy_analysis.paths import native_path
    from fv_energy_analysis.report import write_json
    if not config.get("failure_reference_batch"):
        raise ValueError("All inputs were rejected during OptH. See its per-case diagnostics; failure_reference_batch is required to export a parent-referenced rejection-only cohort.")
    reference = Path(config["failure_reference_batch"])
    reference_manifest = json.loads((reference / "batch_manifest.json").read_text(encoding="utf-8"))
    parent = native_path(hydrogen["parent_optimized_pdb"])
    if sha256(parent.read_bytes()).digest() != sha256(native_path(reference_manifest["parent"]).read_bytes()).digest():
        raise ValueError("Rejected-cohort parent differs from the reference batch")
    output = Path(config["output"])
    payload = deepcopy(json.loads((reference / "batch_summary_grouped_fingerprints.json").read_text(encoding="utf-8")))
    payload["rows"] = []
    payload["title"] = "Ligand preparation rejections"
    payload["summary_subtitle"] = f"All {hydrogen['candidate_count']} inputs were rejected; no mutant scientific scores were calculated."
    payload["source_note"] = "Parent baseline reused from " + str(reference) + ". Rejected mutant metrics remain blank."
    write_json(output / "batch_summary_grouped_fingerprints.json", payload)
    with (output / "batch_summary_grouped_fingerprints.csv").open("w", encoding="utf-8", newline="") as stream:
        csv.DictWriter(stream, fieldnames=[c["key"] for c in payload["columns"]]).writeheader()
    with (output / "batch_summary.csv").open("w", encoding="utf-8", newline="") as stream:
        csv.DictWriter(stream, fieldnames=["candidate", "status", "error"]).writeheader()
    manifest = {"status": "complete", "candidate_count": 0, "completed_count": 0, "failure_count": 0,
                "candidates": [], "preparation_only": True, "parent": str(parent),
                "params": config["params"], "ligand_template_pdb": config["ligand_template_pdb"],
                "reference_batch": str(reference),
                "source_opth_manifest": str(Path(config["hydrogen_output"]) / "hydrogen_batch_manifest.json")}
    write_json(output / "batch_manifest.json", manifest)
    merge_preparation_failures(config, hydrogen)


def run(config):
    configure_workspace_imports(config)
    output = Path(config["output"])
    options = scientific_options(config)
    prior_status = output / "full_pipeline_status.json"
    if prior_status.is_file():
        previous = json.loads(prior_status.read_text(encoding="utf-8")).get("scientific_options")
        if previous is not None and previous != options:
            raise ValueError("Scientific options changed for an existing run. Choose a new output folder; report-section changes alone do not require new analysis.")
    output.mkdir(parents=True, exist_ok=True)
    lock = output / "full_pipeline.lock"
    try:
        descriptor = os.open(lock, os.O_CREAT | os.O_EXCL | os.O_WRONLY)
    except FileExistsError as error:
        raise RuntimeError(f"Another controller owns {lock}; inspect that recorded PID before removing a stale lock") from error
    with os.fdopen(descriptor, "w", encoding="utf-8") as stream:
        json.dump({"pid": os.getpid(), "started_utc": utc_now()}, stream)
    status_path = output / "full_pipeline_status.json"
    workbook = Path(config.get("expected_xlsx", final_report_folder(config) / "batch_summary_grouped_fingerprints.xlsx"))
    export_xlsx = bool(config.get("export_xlsx", True))
    status = {"status": "running", "stage": "starting", "started_utc": utc_now(),
              "pid": os.getpid(), "output": str(output), "expected_xlsx": str(workbook) if export_xlsx else None,
              "log": str(output / "full_pipeline.log"), "completed_stages": [],
              "parallel_jobs": config["jobs"], "fresh_start": bool(config.get("fresh_start")),
              "weights": config.get("weights", "ligand"),
              "scoring_options": config.get("weights_identity"),
              "hydrogen_scoring_options": config.get("hydrogen_weights_identity"),
              "analyses": config.get("analyses", {}),
              "scientific_options": options,
              "requested_count": config.get("requested_count"),
              "available_input_count": config.get("expected_count"),
              "missing_expected_inputs": config.get("missing_expected_inputs", []),
              "params": config["params"], "params_sha256": sha256(Path(config["params"]).read_bytes()).hexdigest()}

    def save_status():
        temporary = status_path.with_suffix(".json.tmp")
        temporary.write_text(json.dumps(status, indent=2) + "\n", encoding="utf-8")
        temporary.replace(status_path)

    environment = os.environ.copy()
    package_container = str(Path(__file__).resolve().parents[1])
    environment["PYTHONPATH"] = os.pathsep.join(dict.fromkeys([package_container, *filter(None, environment.get("PYTHONPATH", "").split(os.pathsep))]))
    environment["PYTHONUNBUFFERED"] = "1"
    for key in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
        environment[key] = "1"
    flags = subprocess.CREATE_NO_WINDOW if os.name == "nt" else 0
    save_status()
    try:
        with (output / "full_pipeline.log").open("a", encoding="utf-8", buffering=1) as log:
            def stage(name, command, allowed_exit_codes=(0,)):
                status.update(stage=name, stage_started_utc=utc_now())
                save_status()
                log.write(f"\n[{utc_now()}] START {name}\n{command!r}\n")
                log.flush()
                merge_output, merge_identity = None, None
                if name in {"combine_full_csv", "e3_source_inventory"}:
                    merge_output = Path(command[command.index("--output") + 1])
                    merge_identity = combination_identity(command)
                    if reuse_combination(merge_output, merge_identity):
                        status["completed_stages"].append(name + "_reused")
                        save_status()
                        log.write(f"[{utc_now()}] REUSED {name}; source/protocol and all saved table hashes match\n")
                        return
                process = subprocess.run(command, cwd=config.get("workspace", str(output.parent)), env=environment,
                                         stdin=subprocess.DEVNULL, stdout=log, stderr=subprocess.STDOUT,
                                         creationflags=flags, check=False)
                if process.returncode not in allowed_exit_codes:
                    raise subprocess.CalledProcessError(process.returncode, command)
                if merge_output is not None:
                    (merge_output / "full_pipeline_combination_identity.json").write_text(json.dumps({
                        "identity": merge_identity, "artifacts": combination_artifacts(merge_output)}, indent=2) + "\n", encoding="utf-8")
                status["completed_stages"].append(name)
                save_status()
                log.write(f"[{utc_now()}] COMPLETE {name}\n")

            has_valid_mutants = True
            for name, command in preflight_commands(config):
                stage(name, command)
            if not config.get("postprocessing_only"):
                hydrogen = None
                if config.get("hydrogen_output"):
                    stage("hydrogen_optimization", hydrogen_command(config), allowed_exit_codes=(0, 1))
                    hydrogen = json.loads((Path(config["hydrogen_output"]) / "hydrogen_batch_manifest.json").read_text())
                    if hydrogen["status"] not in {"complete", "complete_with_failures"}:
                        raise ValueError("OptH preparation did not complete; inspect its manifest and per-case diagnostics")
                    config = optimized_analysis_config(config)
                    status["optimized_structures"] = config["input_folder"]
                    save_status()
                if hydrogen is not None and not hydrogen["completed_count"]:
                    initialize_rejected_cohort(config, hydrogen)
                    has_valid_mutants = False
                    status["completed_stages"].append("preparation_rejections_recorded")
                    save_status()
                else:
                    stage("scientific_analysis", batch_command(config))
                    if hydrogen is not None:
                        merge_preparation_failures(config, hydrogen)
                    existing = json.loads((output / "batch_manifest.json").read_text(encoding="utf-8"))
                    if existing["status"] not in {"complete", "complete_with_failures"}:
                        raise ValueError("Scientific batch did not complete; inspect batch_manifest.json")
                    has_valid_mutants = int(existing["completed_count"]) > 0
            else:
                existing = json.loads((output / "batch_manifest.json").read_text(encoding="utf-8"))
                if existing["status"] not in {"complete", "complete_with_failures"}:
                    raise ValueError("Scientific stage is not complete; refusing postprocessing-only mode")
                if config.get("hydrogen_output"):
                    hydrogen = json.loads((Path(config["hydrogen_output"]) / "hydrogen_batch_manifest.json").read_text())
                    if hydrogen["status"] not in {"complete", "complete_with_failures"}:
                        raise ValueError("Saved OptH stage is not complete; refusing postprocessing-only mode")
                    merge_preparation_failures(config, hydrogen)
                    existing = json.loads((output / "batch_manifest.json").read_text(encoding="utf-8"))
                    status["completed_stages"].append("hydrogen_optimization_reused")
                status["completed_stages"].append("scientific_analysis")
                has_valid_mutants = int(existing["completed_count"]) > 0
            if has_valid_mutants:
                if config.get("skip_discovery_studio"):
                    status["unavailable_tools"] = ["Discovery Studio"]
                    manifest = json.loads((output / "batch_manifest.json").read_text(encoding="utf-8"))
                    manifest["discovery_studio"] = {"status": "unavailable", "reason": "Explicitly skipped; DS-based E2 features were not assessed"}
                    (output / "batch_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
                    save_status()
                else:
                    ds_identity = discovery_studio_identity(config)
                    identity_file = output / "full_pipeline_ds_identity.json"
                    reuse = identity_file.is_file() and json.loads(identity_file.read_text(encoding="utf-8")) == ds_identity
                    previous_ds = json.loads((output / "batch_manifest.json").read_text(encoding="utf-8")).get("discovery_studio", {})
                    stage("discovery_studio", discovery_studio_command(config, reuse=reuse))
                    if reuse or config.get("reuse_discovery_studio"):
                        saved_manifest = json.loads((output / "batch_manifest.json").read_text(encoding="utf-8"))
                        for key in ("build_info", "installation_root", "version"):
                            if key in previous_ds:
                                saved_manifest["discovery_studio"].setdefault(key, previous_ds[key])
                        (output / "batch_manifest.json").write_text(json.dumps(saved_manifest, indent=2) + "\n", encoding="utf-8")
                    identity_file.write_text(json.dumps(ds_identity, indent=2) + "\n", encoding="utf-8")
                stage("fingerprints_and_scores", fingerprint_command(config))
            for name, command in finalization_commands(config):
                stage(name, command)
            if config.get("select_report"):
                command = [sys.executable, "-u", "-m", "fv_energy_analysis.workbook_validation", str(final_report_folder(config)),
                           "--formats", ",".join(config["report_formats"])]
                for key, flag in (("report_sections", "--sections"), ("report_columns", "--columns"), ("report_output", "--report-output")):
                    if config.get(key) is not None:
                        value = config[key]
                        command.extend([flag, ",".join(value) if isinstance(value, list) else str(value)])
                stage("selected_report", command)
                selected_folder = Path(config.get("report_output", final_report_folder(config) / "selected_report"))
                status["selected_report_folder"] = str(selected_folder)
                workbook = selected_folder / "batch_summary_grouped_fingerprints.xlsx"
            elif export_xlsx:
                stage("excel_export", [sys.executable, "-u", "-m", "fv_energy_analysis.workbook_validation",
                                       str(final_report_folder(config)), "--output", str(workbook)])
                if not workbook.is_file():
                    raise FileNotFoundError(f"Exporter did not produce the expected workbook: {workbook}")
            if config.get("combined_output"):
                status["combined_output"] = config["combined_output"]
                combined_path = Path(config["combined_output"]) / "combined_manifest.json"
                combined = json.loads(combined_path.read_text(encoding="utf-8"))
                status["combined_status"] = combined["status"]
                status["combined_csv"] = str(Path(config["combined_output"]) / "batch_summary_grouped_fingerprints.csv")
            status["final_report_folder"] = str(final_report_folder(config))
            status["final_csv"] = str(final_report_folder(config) / "batch_summary_grouped_fingerprints.csv")
            manifest = json.loads((output / "batch_manifest.json").read_text(encoding="utf-8"))
            failures = int(manifest.get("failure_count", 0))
            final_status = "complete_with_failures" if failures else "complete_with_unavailable_tools" if status.get("unavailable_tools") else "complete"
            status.update(status=final_status, stage="finished",
                          completed_utc=utc_now(), scientific_failure_count=failures,
                          completed_count=manifest.get("completed_count"), workbook=str(workbook) if export_xlsx else None)
            save_status()
            log.write(f"[{utc_now()}] FINISHED {status['status']}\n")
    except BaseException as error:
        status.update(status="failed", completed_utc=utc_now(), error=str(error))
        save_status()
        with (output / "full_pipeline.log").open("a", encoding="utf-8") as log:
            traceback.print_exc(file=log)
        raise
    finally:
        lock.unlink(missing_ok=True)


def scientific_options(config):
    """Keep presentation changes separate from settings that change results."""
    defaults = {"weights": "ligand", "legacy_pre_talaris": None, "hydrogen_weights": None,
                "hydrogen_legacy_pre_talaris": None, "optimize_hydrogens": True,
                "scheme": "imgt", "cdr_definition": "imgt", "fit_region": "fv", "rmsd_atoms": "backbone",
                "rotamer_rmsd_cutoff": 3.0, "hcdr3_tip_sequence": "YYGSS",
                "network_energy_presence_cutoff": 0.1, "network_energy_change_cutoff": 0.25,
                "arpeggio_scope": "neighbourhood", "arpeggio_neighbourhood_radius": 6.0,
                "parent_chain_map": None, "mutant_chain_map": None,
                "parent_ligand_chain": None, "mutant_ligand_chain": None, "ligand_chain": "d"}
    result = {key: config.get(key, default) for key, default in defaults.items()}
    from .pyrosetta_score import resolve_legacy_pre_talaris
    result["hydrogen_weights"] = config.get("hydrogen_weights") or result["weights"]
    result["legacy_pre_talaris"] = resolve_legacy_pre_talaris(result["weights"], config.get("legacy_pre_talaris"))
    result["hydrogen_legacy_pre_talaris"] = resolve_legacy_pre_talaris(result["hydrogen_weights"], config.get("hydrogen_legacy_pre_talaris", config.get("legacy_pre_talaris")))
    for key in ("weights_identity", "hydrogen_weights_identity"):
        if key in config:
            result[key] = config[key]
    selected = config.get("analyses", {})
    result["analyses"] = {key: selected.get(key, True) for key in ("arpeggio", "plip", "network", "rotamers")}
    for key in ("parent", "params", "ligand_template_pdb"):
        if config.get(key) and Path(config[key]).is_file():
            result[key + "_sha256"] = sha256(Path(config[key]).read_bytes()).hexdigest()
    return result


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("config", type=Path, nargs="?")
    parser.add_argument("--describe-options", action="store_true", help="List accepted inputs, analysis switches, defaults and reporting controls")
    parser.add_argument("--preflight", action="store_true", help="Validate configuration paths without running analysis")
    arguments = parser.parse_args(argv)
    if arguments.describe_options:
        from .pipeline_options import describe_options
        print(json.dumps(describe_options(), indent=2))
        return 0
    if arguments.config is None:
        parser.error("Supply a configuration file or --describe-options")
    configuration = read_config(arguments.config)
    if arguments.preflight:
        environment = os.environ.copy()
        environment["PYTHONPATH"] = os.pathsep.join([str(Path(__file__).resolve().parents[1]), environment.get("PYTHONPATH", "")])
        for name, command in preflight_commands(configuration):
            print(f"Checking {name}...", flush=True)
            subprocess.run(command, cwd=configuration["workspace"], env=environment, check=True)
        print("Configuration and dependencies are valid; no analysis was started.")
    else:
        run(configuration)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
