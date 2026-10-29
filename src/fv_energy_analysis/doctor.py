"""Check whether the full analysis environment is ready."""

from __future__ import annotations

import argparse
from datetime import datetime, timezone
import importlib.metadata
import importlib.util
import json
from pathlib import Path
import shutil
import sys

from .resources import molfile_to_params_path, resource_path


PACKAGES = {
    "numpy": "numpy",
    "gemmi": "gemmi",
    "biopython": "Bio",
    "lxml": "lxml",
    "rdkit": "rdkit",
    "pdbe-arpeggio": "arpeggio",
    "plip": "plip",
    "openbabel-wheel": "openbabel",
    "abnumber": "abnumber",
    "anarcii": "anarcii",
    "pyrosetta": "pyrosetta",
}
PROFILES = {
    "core": {"numpy": "numpy"},
    "reporting": {"numpy": "numpy", "gemmi": "gemmi", "openpyxl": "openpyxl"},
    "analysis": PACKAGES,
}


def environment_report(profile: str = "analysis", discovery_studio_root: str | None = None,
                       *, skip_tools: tuple[str, ...] = ()) -> tuple[list[str], list[str]]:
    lines, problems = [], []
    lines.append(f"Python: {sys.version.split()[0]} ({sys.executable})")
    lines.append(f"Dependency profile: {profile}")
    packages = dict(PROFILES[profile])
    if "arpeggio" in skip_tools:
        packages.pop("pdbe-arpeggio", None)
    if "plip" in skip_tools:
        packages.pop("plip", None)
    if {"arpeggio", "plip"} <= set(skip_tools):
        packages.pop("openbabel-wheel", None)
        packages.pop("lxml", None)
    for distribution, module in packages.items():
        available = importlib.util.find_spec(module) is not None
        try:
            version = importlib.metadata.version(distribution) if available else "missing"
        except importlib.metadata.PackageNotFoundError:
            version = "installed (version unavailable)" if available else "missing"
        lines.append(f"{distribution}: {version}")
        if not available:
            problems.append(f"Missing Python package: {distribution}")
    if profile == "analysis":
        for name in ("pdbe-arpeggio", "plip"):
            if (name == "pdbe-arpeggio" and "arpeggio" in skip_tools) or name in skip_tools:
                lines.append(f"{name} executable: not required (analysis disabled)")
                continue
            executable = shutil.which(name)
            lines.append(f"{name} executable: {executable or 'missing'}")
            if not executable:
                problems.append(f"{name} executable is not on PATH")
        # Distribution metadata alone does not establish binary binding compatibility.
        if not {"arpeggio", "plip"} <= set(skip_tools):
            try:
                from openbabel import openbabel
                if "plip" not in skip_tools:
                    from plip.structure.preparation import PDBComplex
                    assert PDBComplex
                lines.append(f"Open Babel binding check: {openbabel.OBReleaseVersion()}")
            except (ImportError, OSError, RuntimeError, AssertionError) as error:
                problems.append(f"PLIP/Open Babel runtime check failed: {error}")

    required_files = (
        "scripts/prepare_ligand_complexes.py",
        "scripts/pdb_to_arpeggio_cif.py",
    )
    for relative in required_files:
        try:
            resource_path(relative)
            lines.append(f"package resource {relative}: ok")
        except FileNotFoundError as error:
            lines.append(f"package resource {relative}: missing")
            problems.append(str(error))
    try:
        params_script = molfile_to_params_path()
        lines.append(f"Rosetta molfile_to_params.py: {params_script}")
    except FileNotFoundError:
        lines.append(
            "Rosetta molfile_to_params.py: not configured "
            "(needed only when generating new ligand params; set FV_MOLFILE_TO_PARAMS)"
        )
    if discovery_studio_root:
        from .paths import native_path
        ds_root = native_path(discovery_studio_root)
        perl = ds_root / "bin" / "perl.bat"
        lines.append(f"Native Discovery Studio Perl: {perl}")
        if not perl.is_file():
            problems.append(f"Native Discovery Studio executable was not found: {perl}")
    else:
        lines.append("Native Discovery Studio: not checked (separate licensed Windows installation)")
    return lines, problems


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--profile", choices=tuple(PROFILES), default="analysis")
    parser.add_argument("--discovery-studio-root", help="Also check a separately installed native Windows DS root")
    parser.add_argument("--output-json", type=Path, help="Save checks and resolved distribution versions for reproducibility")
    parser.add_argument("--skip-tool", action="append", choices=("arpeggio", "plip"), default=[],
                        help="Do not require a detector explicitly disabled in this run")
    args = parser.parse_args(argv)
    lines, problems = environment_report(args.profile, args.discovery_studio_root, skip_tools=tuple(args.skip_tool))
    if args.output_json:
        inventory = sorted((d.metadata.get("Name", "unknown"), d.version)
                           for d in importlib.metadata.distributions())
        args.output_json.parent.mkdir(parents=True, exist_ok=True)
        args.output_json.write_text(json.dumps({
            "checked_utc": datetime.now(timezone.utc).isoformat(),
            "profile": args.profile, "ready": not problems,
            "python": sys.version, "executable": sys.executable,
            "checks": lines, "problems": problems,
            "distributions": dict(inventory),
        }, indent=2) + "\n", encoding="utf-8")
    print("\n".join(lines))
    if problems:
        print("\nEnvironment is NOT ready:")
        print("\n".join(f"- {problem}" for problem in problems))
        return 1
    print(f"\nEnvironment checks passed for the {args.profile} profile.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
