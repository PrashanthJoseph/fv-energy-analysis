"""Selection and provenance for optional scientific analysis stages."""
from __future__ import annotations

import json
from pathlib import Path


SKIP_OPTIONS = {
    "arpeggio": "skip_arpeggio",
    "plip": "skip_plip",
    "network": "skip_network_analysis",
    "rotamers": "skip_rotamer_analysis",
    "e2_architecture": "skip_e2_architecture",
}


def omitted_e2_architecture():
    return {role: {"e2_architecture_status": "not_run"} for role in ("parent", "mutant")}


def option(options, name, default=None):
    return options.get(name, default) if isinstance(options, dict) else getattr(options, name, default)


def analysis_status(options, *, pending=False):
    return {name: "not_run" if option(options, flag, False) else "pending" if pending else "complete"
            for name, flag in SKIP_OPTIONS.items()}


def add_analysis_arguments(parser):
    for flag in SKIP_OPTIONS.values():
        parser.add_argument("--" + flag.replace("_", "-"), action="store_true",
                            help="Explicitly omit this analysis; unavailable values are not zero counts")


def analysis_flags(options):
    return ["--" + flag.replace("_", "-") for flag in SKIP_OPTIONS.values()
            if option(options, flag, False)]


def validate_analysis_reuse(directory, options):
    manifest = json.loads((Path(directory) / "workflow_manifest.json").read_text(encoding="utf-8"))
    saved = manifest.get("analysis_status", {})
    expected = analysis_status(options)
    if any(saved.get(name, "complete") != status for name, status in expected.items()):
        raise ValueError("Saved analysis selection differs from requested optional stages; choose a new output folder")


def validate_comparison_reuse(directory, options, candidate=None):
    """Do not silently reuse derived metrics computed with different settings."""
    from hashlib import sha256
    from .paths import native_path

    directory = Path(directory)
    path = directory / "final/summary.json"
    summary = json.loads(path.read_text(encoding="utf-8")) if path.is_file() else {}
    checks = {
        "numbering_scheme": ("scheme", "imgt"),
        "cdr_definition": ("cdr_definition", "imgt"),
        "fit_region": ("fit_region", "fv"),
        "rmsd_atoms": ("rmsd_atoms", "backbone"),
        "rotamer_rmsd_cutoff_A": ("rotamer_rmsd_cutoff", 3.0),
        "hcdr3_tip_reference_sequence": ("hcdr3_tip_sequence", "YYGSS"),
        "network_energy_presence_cutoff_REU": ("network_energy_presence_cutoff", 0.1),
        "network_energy_change_cutoff_REU": ("network_energy_change_cutoff", 0.25),
    }
    for field, (argument, historical_default) in checks.items():
        expected = option(options, argument, historical_default)
        if summary.get(field, historical_default) != expected:
            raise ValueError(f"Saved comparison setting differs: {argument}; choose a new output or reuse mutants in a new comparison")
    if candidate is not None:
        manifest = json.loads((directory / "workflow_manifest.json").read_text(encoding="utf-8"))
        candidate = Path(candidate)
        digest = sha256(candidate.read_bytes()).hexdigest()
        saved_digest = manifest.get("mutant_input_sha256")
        if saved_digest and saved_digest != digest:
            raise ValueError("Saved mutant input content differs; choose a new output folder")
        recorded = native_path(manifest.get("mutant", str(candidate)))
        if recorded.resolve() != candidate.resolve() and (not recorded.is_file() or sha256(recorded.read_bytes()).hexdigest() != digest):
            raise ValueError("Saved mutant input differs from requested candidate")


def record_comparison_status(final, options):
    path = Path(final) / "summary.json"
    # A real comparison always writes this file; lightweight runner doubles may not.
    if path.is_file():
        summary = json.loads(path.read_text(encoding="utf-8"))
        summary["analysis_status"] = analysis_status(options)
        path.write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")


def scoring_protocol(params, weights, legacy_pre_talaris=None):
    """Retain legacy built-in metadata while binding custom/overridden scoring."""
    from hashlib import sha256
    from .pyrosetta_score import SCORING_PROTOCOL_VERSION, scoring_weights_identity

    value = {"version": SCORING_PROTOCOL_VERSION, "weights": str(weights),
             "params_sha256": sha256(Path(params).read_bytes()).hexdigest()}
    if legacy_pre_talaris is not None or Path(weights).is_file() or str(weights) != "ligand":
        from importlib.util import find_spec
        package = find_spec("pyrosetta")
        database = Path(package.origin).parent / "database" if package and package.origin else None
        identity = scoring_weights_identity(weights, legacy_pre_talaris, database=database)
        value.update(legacy_pre_talaris=identity["legacy_pre_talaris"], weights_identity=identity)
        from importlib.metadata import PackageNotFoundError, version
        try:
            value["pyrosetta_version"] = version("pyrosetta")
        except PackageNotFoundError:
            pass  # Mocked unit tests may exercise protocol construction without the engine.
    return value
