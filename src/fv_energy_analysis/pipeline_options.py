"""Small, explicit configuration contract for the complete energy workflow."""
from __future__ import annotations

import math


ANALYSES = {
    "arpeggio": True,
    "plip": True,
    "discovery_studio": True,
    "e2_architecture": True,
    "e3_discrimination": False,
    "network": True,
    "rotamers": True,
    "hcdr3_plddt": True,
    "boltz2_scores": True,
    "gnina_scores": True,
}
CHOICES = {
    "scheme": ("imgt", "chothia", "kabat", "aho"),
    "cdr_definition": ("imgt", "chothia", "kabat", "north"),
    "fit_region": ("fv", "framework"),
    "rmsd_atoms": ("backbone", "ca"),
    "arpeggio_scope": ("neighbourhood", "full"),
}
NUMERICAL = {
    "rotamer_rmsd_cutoff": (3.0, True),
    "network_energy_presence_cutoff": (0.1, True),
    "network_energy_change_cutoff": (0.25, True),
    "arpeggio_neighbourhood_radius": (6.0, False),
    "e3_overlap_warning": (0.4, False),
}


def normalize_options(config: dict) -> dict:
    """Validate knobs without importing scientific applications or running work."""
    config = dict(config)
    for key in ("optimize_hydrogens", "export_xlsx", "recursive", "postprocessing_only",
                "stop_on_error", "e3_discrimination", "skip_discovery_studio",
                "fresh_start", "reuse_discovery_studio", "include_empty_metadata_sections"):
        if key in config and not isinstance(config[key], bool):
            raise ValueError(f"{key} must be true or false")
    for key in ("jobs", "hydrogen_seed", "expected_count"):
        if key in config and (isinstance(config[key], bool) or not isinstance(config[key], int)
                              or config[key] < (1 if key == "jobs" else 0)):
            raise ValueError(f"{key} must be a {'positive' if key == 'jobs' else 'non-negative'} integer")
    for key in ("weights", "hydrogen_weights"):
        if key == "hydrogen_weights" and config.get(key) is None:
            config.pop(key, None)
            continue
        if key in config and (not isinstance(config[key], str) or not config[key].strip()):
            raise ValueError(f"{key} must name a Rosetta weight set or an existing .wts file")
    requested = config.get("analyses", {})
    if not isinstance(requested, dict):
        raise ValueError("analyses must be an object of named true/false switches")
    unknown = set(requested) - set(ANALYSES)
    if unknown:
        raise ValueError("Unknown analyses: " + ", ".join(sorted(unknown)))
    if any(not isinstance(value, bool) for value in requested.values()):
        raise ValueError("Each analyses switch must be true or false")
    options = {**ANALYSES, **requested}
    # Preserve existing configuration files while letting an explicit named
    # analysis switch take precedence over the older top-level switches.
    if "discovery_studio" not in requested:
        options["discovery_studio"] = not config.get("skip_discovery_studio", False)
    if "e2_architecture" not in requested:
        options["e2_architecture"] = options["discovery_studio"]
    if "e3_discrimination" not in requested:
        options["e3_discrimination"] = config.get("e3_discrimination", False)
    if options["e2_architecture"] and not options["discovery_studio"]:
        raise ValueError("E2 architecture requires analyses.discovery_studio=true")
    if options["e3_discrimination"] and not options["e2_architecture"]:
        raise ValueError("The integrated E3 report requires DS-based E2 architecture; enable it or run the standalone E3 diagnostic")
    config["analyses"] = options
    config["skip_discovery_studio"] = not options["discovery_studio"]
    config["e3_discrimination"] = options["e3_discrimination"]

    for key, allowed in CHOICES.items():
        if key in config and config[key] not in allowed:
            raise ValueError(f"{key} must be one of: {', '.join(allowed)}")
    for key, (default, allow_zero) in NUMERICAL.items():
        value = config.get(key, default)
        if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value):
            raise ValueError(f"{key} must be a finite number")
        if value < 0 or (not allow_zero and value == 0):
            raise ValueError(f"{key} must be {'non-negative' if allow_zero else 'positive'}")
        if key in config:
            config[key] = float(value)
    for key in ("legacy_pre_talaris", "hydrogen_legacy_pre_talaris"):
        if config.get(key) is not None and not isinstance(config[key], bool):
            raise ValueError(f"{key} must be true, false or null (automatic)")
    for key in ("report_sections", "report_columns"):
        if config.get(key) is not None:
            value = config[key]
            if not isinstance(value, list) or any(not isinstance(v, str) or not v.strip() for v in value):
                raise ValueError(f"{key} must be a list of names, or null for all")
            if len(set(value)) != len(value):
                raise ValueError(f"{key} contains duplicate names")
    formats = config.get("report_formats", ["csv", "xlsx"] if config.get("export_xlsx", True) else ["csv"])
    if not isinstance(formats, list) or not formats or set(formats) - {"csv", "xlsx"}:
        raise ValueError("report_formats must select csv, xlsx, or both")
    config["report_formats"] = list(dict.fromkeys(formats))
    config["export_xlsx"] = "xlsx" in formats
    config["select_report"] = config.get("report_sections") is not None or config.get("report_columns") is not None or bool(config.get("report_output"))
    return config


def describe_options() -> dict:
    """Human-readable machine output for --describe-options."""
    return {
        "required_inputs": ["parent", "input_folder OR mutant", "output", "params", "ligand_template_pdb"],
        "complex_formats": {"with_OptH": ["pdb"], "without_OptH": ["pdb", "cif", "mmcif"]},
        "optional_sources": ["cif_folder", "parent_cif", "scores_csv", "ligand_source", "ligand_smiles"],
        "always_run": ["ligand validation", "Rosetta scoring and H-bonds", "antibody numbering and global-frame RMSDs"],
        "analyses": ANALYSES,
        "scoring": {"weights": "ligand (default), another installed Rosetta weight set, or a custom .wts file",
                    "legacy_pre_talaris": "null=automatic; true/false explicit; required for a custom .wts file",
                    "hydrogen_weights": "defaults to weights; explicit alternative allowed and recorded separately",
                    "hydrogen_legacy_pre_talaris": "defaults to legacy_pre_talaris"},
        "choices": CHOICES,
        "numeric_defaults": {key: value[0] for key, value in NUMERICAL.items()},
        "execution": {"jobs": 2, "recursive": False, "input_glob": "*.pdb", "optimize_hydrogens": True},
        "report": {"report_sections": "list of section names; null=all",
                   "report_columns": "list of column keys; combined with sections as a union; null=no extra keys",
                   "report_formats": ["csv", "xlsx"],
                   "report_output": "defaults to selected_report/ beneath the completed report folder"},
    }
