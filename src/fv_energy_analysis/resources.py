"""Locate runtime files in a source checkout or an installed distribution."""

from __future__ import annotations

import os
from pathlib import Path
import sysconfig


def resource_path(relative: str) -> Path:
    """Return an owned runtime resource, including files installed by pip."""
    package_parent = Path(__file__).resolve().parent.parent
    candidates = (
        package_parent.parent / relative,
        package_parent / "share" / "fv-energy-analysis" / relative,
        Path(sysconfig.get_path("data")) / "share" / "fv-energy-analysis" / relative,
    )
    for candidate in candidates:
        if candidate.is_file():
            return candidate
    searched = ", ".join(str(path) for path in dict.fromkeys(candidates))
    raise FileNotFoundError(f"Package resource {relative!r} was not found. Searched: {searched}")


def molfile_to_params_path(explicit: str | Path | None = None) -> Path:
    """Resolve the separately licensed Rosetta helper only when it is needed."""
    configured = explicit or os.environ.get("FV_MOLFILE_TO_PARAMS")
    candidate = (
        Path(configured).expanduser().resolve()
        if configured
        else Path(__file__).resolve().parents[2] / "vendor" / "rosetta_scripts" / "molfile_to_params.py"
    )
    if not candidate.is_file():
        raise FileNotFoundError(
            "Rosetta molfile_to_params.py was not found. Supply --molfile-to-params "
            "or set FV_MOLFILE_TO_PARAMS to your licensed Rosetta script, or reuse "
            "existing --params and --ligand-template-pdb files. "
            f"Checked: {candidate}"
        )
    return candidate
