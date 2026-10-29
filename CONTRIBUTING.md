# Development

Use Python 3.10 or newer. From the source checkout:

```bash
python -m pip install -e '.[test,reporting]'
python -m pytest -q
```

The core tests and Excel export do not require PyRosetta, Arpeggio, PLIP or
Discovery Studio. Those tools are needed for their corresponding analysis
stages and integration tests.

## Keep changes small

- Tool readers normalize their own outputs; they do not run other tools.
- Scientific calculations belong in the analysis modules.
- `fingerprint_counts.py` reads saved detector outputs.
- `fingerprint_reports.py` owns the summary columns and report generation.
- `spreadsheet.py` formats the optional Excel output.
- `workbook_validation.py` validates saved data, exports the current complete
  workbook and verifies typed cells, filters, panes and failure rows.
- `full_pipeline.py` owns the portable whole-folder stage orchestration.
- `fingerprints.py` is the command-line entry point.
- Cache the parent and reuse completed candidate artifacts. Report changes do
  not require rescoring or re-running detectors.

Use one canonical grouped summary. Keep full detector data in the individual
candidate outputs. Avoid adding new export variants for presentation changes.
Add tests for scientific mappings, file formats and package resources.

## Build and release

```bash
python -m pip install build
python -m build
```

CI builds the source distribution and wheel, then runs the tests using the
installed wheel from a directory outside the original checkout. Tests must
import packaged modules; never import ignored files under `tmp`.

Runtime scripts and report definitions are included in the wheel.
Data and run outputs, local environments, diagnostic
logs and licensed Rosetta utilities are excluded. Install Rosetta/PyRosetta
separately according to their terms. For new ligand parameters, provide
`--molfile-to-params PATH` or set `FV_MOLFILE_TO_PARAMS`; existing ligand
parameters can be supplied directly.

The owner has left the original project license undecided. Do not add an
open-source license without that decision. Preserve the
[third-party notices](THIRD_PARTY_NOTICES.md). PyRosetta and native
Discovery Studio are external installations and are not bundled. Dependency
extras, installer behavior and environment snapshots are described in
[docs/dependencies.md](docs/dependencies.md).
