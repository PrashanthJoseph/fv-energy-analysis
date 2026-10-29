#!/usr/bin/env bash
set -euo pipefail

project_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
venv_dir="${1:-$project_dir/.venv}"

if [[ ! -x "$venv_dir/bin/python" ]]; then
    python3 -m venv "$venv_dir"
fi
source "$venv_dir/bin/activate"
python -m pip install --upgrade pip setuptools wheel
python -m pip install -e "$project_dir[analysis,test,reporting]"
bash "$project_dir/scripts/install_plip_wsl.sh" "$venv_dir"
pyrosetta_requirement="pyrosetta${FV_PYROSETTA_VERSION:+==$FV_PYROSETTA_VERSION}"
python -m pip install "$pyrosetta_requirement" \
  --find-links https://west.rosettacommons.org/pyrosetta/quarterly/release
python -m fv_energy_analysis.doctor --profile analysis \
  --output-json "$venv_dir/fv-energy-environment.json"
python -m pip freeze --all > "$venv_dir/fv-energy-requirements.lock.txt"

echo ""
echo "Installation complete. Activate with: source '$venv_dir/bin/activate'"
echo "Resolved versions: $venv_dir/fv-energy-requirements.lock.txt"
echo "Environment checks: $venv_dir/fv-energy-environment.json"
