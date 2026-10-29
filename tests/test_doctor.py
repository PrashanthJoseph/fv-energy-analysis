"""Dependency checks distinguish portable reporting from licensed analysis tools."""
from pathlib import Path
import json

from fv_energy_analysis import doctor


def available_environment(monkeypatch):
    monkeypatch.setattr(doctor.importlib.util, "find_spec", lambda name: object())
    monkeypatch.setattr(doctor.importlib.metadata, "version", lambda name: "1.0")
    monkeypatch.setattr(doctor, "resource_path", lambda name: Path(__file__))
    monkeypatch.setattr(doctor, "molfile_to_params_path", lambda: Path(__file__))


def test_reporting_does_not_require_proprietary_tools_or_detector_executables(monkeypatch):
    available_environment(monkeypatch)
    monkeypatch.setattr(doctor.shutil, "which", lambda name: None)
    lines, problems = doctor.environment_report("reporting")
    assert not problems
    assert "openpyxl: 1.0" in lines
    assert not any(line.startswith("pyrosetta:") for line in lines)
    assert any("Discovery Studio: not checked" in line for line in lines)


def test_missing_reporting_package_is_reported(monkeypatch):
    available_environment(monkeypatch)
    monkeypatch.setattr(doctor.importlib.util, "find_spec", lambda name: None if name == "openpyxl" else object())
    _, problems = doctor.environment_report("reporting")
    assert problems == ["Missing Python package: openpyxl"]


def test_explicit_missing_discovery_studio_is_an_error(monkeypatch, tmp_path):
    available_environment(monkeypatch)
    _, problems = doctor.environment_report("core", str(tmp_path / "missing_DS"))
    assert len(problems) == 1 and "Discovery Studio executable" in problems[0]


def test_json_snapshot_captures_resolved_versions_and_profile(monkeypatch, tmp_path):
    available_environment(monkeypatch)
    class Distribution:
        metadata = {"Name": "example-runtime"}
        version = "2.3.4"
    monkeypatch.setattr(doctor.importlib.metadata, "distributions", lambda: [Distribution()])
    output = tmp_path / "environment.json"
    assert doctor.main(["--profile", "core", "--output-json", str(output)]) == 0
    result = json.loads(output.read_text())
    assert result["ready"] is True and result["profile"] == "core"
    assert result["distributions"] == {"example-runtime": "2.3.4"}


def test_disabled_detectors_do_not_require_their_packages_or_executables(monkeypatch):
    available_environment(monkeypatch)
    excluded = {"plip", "arpeggio", "openbabel", "lxml"}
    monkeypatch.setattr(doctor.importlib.util, "find_spec", lambda name: None if name in excluded else object())
    monkeypatch.setattr(doctor.shutil, "which", lambda name: None)
    lines, problems = doctor.environment_report("analysis", skip_tools=("arpeggio", "plip"))
    assert not problems
    assert "plip executable: not required (analysis disabled)" in lines
