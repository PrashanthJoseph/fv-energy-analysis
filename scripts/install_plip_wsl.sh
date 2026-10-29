#!/usr/bin/env bash
set -euo pipefail

project_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
venv_dir="${1:-$project_dir/.venv}"
python_bin="$venv_dir/bin/python"

if [[ ! -x "$python_bin" ]]; then
    echo "Python environment not found: $venv_dir" >&2
    echo "Create it with scripts/bootstrap_wsl.sh or pass its path as the first argument." >&2
    exit 2
fi

plip_version="3.0.1"
echo "Installing/verifying PLIP in: $venv_dir"

# Read the bindings' requirements from installed package metadata so this
# standalone installer cannot drift from pyproject.toml's analysis extra.
"$python_bin" - <<'PY'
from importlib.metadata import PackageNotFoundError, requires
import subprocess
import sys

try:
    dependencies = requires("fv-energy-analysis") or []
except PackageNotFoundError:
    raise SystemExit("Install fv-energy-analysis[analysis] before adding PLIP.")
runtime = [item.split(";", 1)[0].strip() for item in dependencies
           if item.lower().startswith(("lxml", "openbabel-wheel"))]
if len(runtime) != 2:
    raise SystemExit("Expected lxml and openbabel-wheel in package dependency metadata.")
subprocess.check_call([sys.executable, "-m", "pip", "install", *runtime])
PY

# PLIP 3.0.1 declares openbabel>=3.2.0 by distribution name. PyPI currently
# supplies the importable binary bindings separately as openbabel-wheel. Avoid
# trying to compile a second Open Babel copy and use the verified bindings.
"$python_bin" -m pip install --no-deps "plip==$plip_version"

"$python_bin" - <<'PY'
from openbabel import openbabel
from plip.basic import config
from plip.structure.preparation import PDBComplex

assert openbabel.OBReleaseVersion()
assert PDBComplex
print(f"Open Babel Python bindings: {openbabel.OBReleaseVersion()}")
print(f"PLIP: {config.__version__}")
PY

"$venv_dir/bin/plip" --help >/dev/null
echo "PLIP installation verified. This script is idempotent and only needs to be run once per environment."
