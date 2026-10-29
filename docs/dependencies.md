# Dependencies and reproducible installation

`pyproject.toml` is the source of truth for ordinary Python dependencies.
`requirements-analysis.txt` delegates to its `analysis` and `reporting` extras;
it does not maintain a second list of versions.

| Component | Installation | Required for |
|---|---|---|
| NumPy | Base package dependency | Geometry, tables and analysis modules |
| Gemmi, RDKit, Biopython, AbNumber, ANARCII, PDBe Arpeggio, lxml, Open Babel bindings | `.[analysis]` | Scientific preparation, numbering, fingerprints and optional SASA |
| openpyxl and Gemmi | `.[reporting]` | Portable XLSX export and original-CIF confidence |
| pytest and Biopython | `.[test]` | Automated tests, including C18/SASA numerical checks |
| PLIP 3.0.1 | `scripts/install_plip_wsl.sh` | Independent ligand interaction fingerprint |
| PyRosetta | Official quarterly wheel index in the bootstrap; or an authorized local build | OptH and Rosetta scoring |
| Native BIOVIA Discovery Studio | Separate licensed Windows installation | DS fingerprints and DS-derived E2 recognition features |
| Rosetta `molfile_to_params.py` | Existing licensed Rosetta installation, configured explicitly | Generating new ligand parameters only |

PyRosetta and Discovery Studio are not redistributed in this repository.
The complete four-tool workflow uses a Windows controller for native DS and
a Linux/WSL scientific Python environment for PyRosetta and the other detectors.
Linux-only runs can explicitly omit DS; DS-dependent E2 features then remain
unavailable. This is not the complete DS-based production screen.

## Scientific environment

From the source checkout under Linux or WSL:

Python 3.12 is the tested scientific environment. The lightweight package and
CI also support Python 3.10; installing PyRosetta additionally requires a wheel
available for the chosen Python version and operating system.

```bash
bash scripts/bootstrap_wsl.sh
source .venv/bin/activate
fv-energy-doctor --profile analysis
```

The bootstrap creates or reuses that virtual environment, installs the package
with its analysis, reporting and test extras, installs PLIP, and obtains
PyRosetta from the official RosettaCommons quarterly release index. It does
not install native Discovery Studio or bypass any license requirement.

PLIP's dependency metadata requests an `openbabel` distribution which differs
from the available binary distribution name, `openbabel-wheel`. The installer
uses the package's declared lxml/Open Babel requirements, installs pinned
PLIP with `--no-deps`, then checks the actual Open Babel bindings, PLIP
preparation API and executable. Do not replace this with an unrestricted
`pip install plip`, which may attempt to build another Open Babel copy.
Because of this distribution-name mismatch, `pip check` can report a PLIP
metadata requirement even when the explicitly checked bindings work.

The bootstrap records:

- `.venv/fv-energy-requirements.lock.txt`: the exact resolved environment from
  `pip freeze --all`;
- `.venv/fv-energy-environment.json`: Python and distribution versions plus
  the dependency checks.

These describe one platform/environment. They are not a universal lock file
for Windows, Linux, every Python version or a different PyRosetta release.
The editable-project path in `pip freeze` also belongs to that checkout.
Archive both records with a scientific run's configuration. The analysis
manifests separately record the scoring protocol and relevant code/input hashes.

To select an authorized exact PyRosetta build instead of the quarterly
index's available build, set `FV_PYROSETTA_VERSION` when running the bootstrap.
An existing installed PyRosetta is retained unless a different version is
requested. A recipient without download access must install an authorized
compatible wheel in the same environment before running the doctor.

## Windows controller and reporting environment

Use a separate Windows Python environment, because native DS is run on Windows:

```powershell
py -3.12 -m venv .venv-windows
.\.venv-windows\Scripts\python -m pip install ".[reporting]"
.\.venv-windows\Scripts\fv-energy-doctor --profile reporting
```

The full runner's JSON configuration identifies the WSL Python executable and
the native DS installation root. It does not install or assume a Codex runtime,
Node.js, or the original author's environment paths. See the portable example
and command in the main README.

For a Linux-only subset, set `science_backend` to `native` and explicitly set
`skip_discovery_studio` to `true`; `science_python` defaults to the interpreter
running the controller. Do not describe that subset as a DS-validated E2 screen.

For a Windows DS installation, optionally check its executable path too:

```powershell
.\.venv-windows\Scripts\fv-energy-doctor --profile reporting --discovery-studio-root "C:/Program Files/BIOVIA/Discovery Studio 2024"
```

## Lightweight development and saved-report export

```bash
python -m pip install -e '.[test,reporting]'
python -m pytest -q
fv-workbook /path/to/completed/batch
```

These tasks do not require PyRosetta or native DS. Scientific tests which need
additional dependencies or local reference coordinates are explicitly skipped
when those prerequisites are unavailable. CI tests an installed wheel outside
the original checkout so ignored local helpers cannot hide missing resources.
