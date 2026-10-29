"""User-facing configuration switches and their computational dependencies."""
import json
from pathlib import Path

import pytest

from fv_energy_analysis import full_pipeline
from fv_energy_analysis.pipeline_options import describe_options, normalize_options


def test_default_analysis_and_report_contract():
    config = normalize_options({})
    assert config["analyses"]["arpeggio"] is True
    assert config["analyses"]["e2_architecture"] is True
    assert config["analyses"]["e3_discrimination"] is False
    assert config["report_formats"] == ["csv", "xlsx"]
    assert config["select_report"] is False


def test_disabled_native_ds_defaults_e2_off_but_never_fakes_e3():
    config = normalize_options({"analyses": {"discovery_studio": False}})
    assert config["skip_discovery_studio"] is True
    assert config["analyses"]["e2_architecture"] is False
    with pytest.raises(ValueError, match="E2 architecture requires"):
        normalize_options({"analyses": {"discovery_studio": False, "e2_architecture": True}})
    with pytest.raises(ValueError, match="integrated E3 report requires"):
        normalize_options({"analyses": {"discovery_studio": False, "e3_discrimination": True}})


@pytest.mark.parametrize("config", [
    {"analyses": {"plipp": True}}, {"analyses": {"plip": "false"}},
    {"jobs": 1.5}, {"jobs": True}, {"recursive": "false"},
    {"fit_region": "local"}, {"rmsd_atoms": "sidechain"},
    {"rotamer_rmsd_cutoff": -1}, {"arpeggio_neighbourhood_radius": 0},
    {"network_energy_change_cutoff": float("nan")},
    {"legacy_pre_talaris": "true"}, {"report_formats": ["pdf"]},
    {"report_sections": "Rosetta"}, {"report_columns": ["candidate", "candidate"]},
])
def test_invalid_knobs_are_rejected(config):
    with pytest.raises(ValueError):
        normalize_options(config)


def test_view_selection_does_not_disable_computation():
    config = normalize_options({"report_sections": ["rosetta", "e2"], "report_formats": ["csv"]})
    assert config["select_report"] is True
    assert config["export_xlsx"] is False
    assert config["analyses"]["arpeggio"] is True
    assert normalize_options({"report_sections": []})["select_report"] is True


def test_forwarded_scientific_choices_and_disabled_detectors():
    config = {"science_backend": "native", "science_python": "/venv/python", "parent": "/in/parent.pdb",
              "input_folder": "/in", "output": "/out", "params": "/in/L.params",
              "ligand_template_pdb": "/in/L.pdb", "jobs": 2, "weights": "ref2015",
              "legacy_pre_talaris": False, "scheme": "kabat", "cdr_definition": "north",
              "fit_region": "framework", "rmsd_atoms": "ca", "rotamer_rmsd_cutoff": 2.0,
              "network_energy_presence_cutoff": 0.2, "network_energy_change_cutoff": 0.5,
              "analyses": {"arpeggio": False, "plip": False, "network": False, "rotamers": False}}
    command = full_pipeline.batch_command(config)
    for flag, value in (("--weights", "ref2015"), ("--scheme", "kabat"), ("--fit-region", "framework"),
                        ("--cdr-definition", "north"), ("--rmsd-atoms", "ca"), ("--rotamer-rmsd-cutoff", "2.0")):
        assert command[command.index(flag) + 1] == value
    assert {"--no-legacy-pre-talaris", "--skip-arpeggio", "--skip-plip", "--skip-network-analysis", "--skip-rotamer-analysis"} <= set(command)


def test_single_mutant_selects_only_its_exact_file(tmp_path):
    for name in ("parent.pdb", "one.pdb", "other.pdb", "L.params", "L.pdb"):
        (tmp_path / name).write_text("fixture")
    path = tmp_path / "config.json"
    path.write_text(json.dumps({"parent": "parent.pdb", "mutant": "one.pdb", "output": "out",
        "params": "L.params", "ligand_template_pdb": "L.pdb", "science_backend": "native",
        "analyses": {"discovery_studio": False}}))
    config = full_pipeline.read_config(path)
    assert config["input_glob"] == "one.pdb" and config["recursive"] is False
    assert config["input_folder"] == str(tmp_path)


def test_options_are_available_without_paths_or_external_tools(capsys):
    assert full_pipeline.main(["--describe-options"]) == 0
    value = json.loads(capsys.readouterr().out)
    assert value["analyses"] == describe_options()["analyses"]
    assert "report_sections" in value["report"]


def test_report_selection_does_not_change_scientific_reuse_identity():
    base = {"weights": "ligand"}
    assert full_pipeline.scientific_options(base) == full_pipeline.scientific_options({
        **base, "report_sections": ["rosetta"], "report_formats": ["csv"]})
    assert full_pipeline.scientific_options(base) == full_pipeline.scientific_options({
        **base, "legacy_pre_talaris": True})
    assert full_pipeline.scientific_options(base) != full_pipeline.scientific_options({
        **base, "fit_region": "framework"})


def test_custom_weight_files_need_explicit_legacy_choice_and_are_hashed(tmp_path):
    for name in ("parent.pdb", "one.pdb", "L.params", "L.pdb"):
        (tmp_path / name).write_text("fixture")
    weight_file = tmp_path / "custom.wts"
    weight_file.write_text("fa_atr 1.0\n")
    config = {"parent": "parent.pdb", "mutant": "one.pdb", "output": "out",
              "params": "L.params", "ligand_template_pdb": "L.pdb", "science_backend": "native",
              "analyses": {"discovery_studio": False}, "weights": "custom.wts"}
    path = tmp_path / "config.json"
    path.write_text(json.dumps(config))
    with pytest.raises(ValueError, match="explicitly for custom"):
        full_pipeline.read_config(path)
    path.write_text(json.dumps({**config, "legacy_pre_talaris": False}))
    original = full_pipeline.read_config(path)
    assert original["weights"] == str(weight_file)
    weight_file.write_text("fa_atr 0.8\n")
    edited = full_pipeline.read_config(path)
    assert original["weights_identity"]["weights_sha256"] != edited["weights_identity"]["weights_sha256"]


def test_disabled_e2_is_forwarded_to_scientific_and_reporting_commands():
    config = {"science_backend": "native", "science_python": "/venv/python", "parent": "/in/p.pdb",
              "input_folder": "/in", "output": "/out", "params": "/in/L.params",
              "ligand_template_pdb": "/in/L.pdb", "jobs": 2,
              "analyses": {"discovery_studio": False, "e2_architecture": False}}
    assert "--skip-e2-architecture" in full_pipeline.batch_command(config)
    report = full_pipeline.fingerprint_command(config)
    assert "--skip-e2-architecture" in report and "--skip-discovery-studio" in report
