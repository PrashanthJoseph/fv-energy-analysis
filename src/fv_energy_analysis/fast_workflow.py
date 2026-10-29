"""Cached single-candidate workflow used by persistent parallel batch workers."""

from __future__ import annotations

import json
import os
import shutil
import subprocess
import sys
import time
from hashlib import sha256
from datetime import datetime, timezone
from pathlib import Path
from types import SimpleNamespace
from typing import Any

from .cli import compare
from .estrogen_architecture import analyze_estrogen_architecture
from .plip import summarize_plip_xml
from .paths import native_path
from .pyrosetta_score import score_and_extract_hbonds, SCORING_PROTOCOL_VERSION
from .resources import resource_path
from .arpeggio import (run_arpeggio, mutation_neighbourhood_selection,
                       filter_arpeggio_contacts, validate_arpeggio_coverage)
from .validation import validate_ligand_against_params
from .analysis_stages import analysis_status, record_comparison_status, validate_analysis_reuse, scoring_protocol, omitted_e2_architecture
from .workflow import Runner, _convert_structure_input, _write_automatic_summary


def _link_or_copy(source: Path, destination: Path) -> None:
    if source.resolve() == destination.resolve():
        if not source.is_file():
            raise FileNotFoundError(source)
        return
    destination.parent.mkdir(parents=True, exist_ok=True)
    if destination.exists():
        destination.unlink()
    try:
        os.link(source, destination)
    except OSError:
        shutil.copy2(source, destination)


def _recorded_pdb(result_dir: Path, role: str) -> Path:
    """Use the saved comparison's exact structure, never an arbitrary glob match."""
    summary = json.loads((result_dir / "final" / "summary.json").read_text(encoding="utf-8"))
    recorded = native_path(summary[role])
    path = recorded if recorded.is_file() else result_dir / "prepared" / recorded.name
    if not path.is_file():
        raise FileNotFoundError(f"Recorded prepared {role} PDB is missing: {summary[role]}")
    return path


def _cache_paths(cache: Path) -> dict[str, Path]:
    manifest = json.loads((cache / "workflow_manifest.json").read_text(encoding="utf-8"))
    analysis_input = native_path(str(manifest.get("parent_analysis_input", "")))
    paths = {
        "prepared": _recorded_pdb(cache, "parent"),
        "analysis_input": analysis_input,
        "pairs": cache / "parent_pair_energies.sc",
        "interface": cache / "parent_interface.sc",
        "hbonds": cache / "parent_hbonds.csv",
        "arpeggio": cache / "arpeggio_parent" / "parent_arpeggio.json",
        "plip": cache / "plip_ligand_interactions.json",
        "validation": cache / "ligand_atom_validation.json",
        "architecture": cache / "estrogen_binding_architecture.json",
    }
    status = manifest.get("analysis_status", {})
    return {name: path for name, path in paths.items()
            if name not in {"arpeggio", "plip"} or status.get(name, "complete") != "not_run"}


def validate_scoring_reuse(result_dir: Path, params: Path, weights: str, legacy_pre_talaris=None) -> None:
    """Never mix legacy/incompatible scientific scores with corrected scores."""
    manifest = json.loads((result_dir / 'workflow_manifest.json').read_text())
    expected = scoring_protocol(params, weights, legacy_pre_talaris)
    if manifest.get('scoring_protocol') != expected:
        raise ValueError(f'Cannot reuse legacy or incompatible scoring from {result_dir}; use a new output folder and rescore with the corrected ligand parameters')


def run_cached_candidate(config: dict[str, Any]) -> dict[str, Any]:
    """Run one mutant while reusing immutable parent artifacts.

    ProcessPoolExecutor reuses this process for multiple calls, so PyRosetta and
    ANARCII's sequence cache remain warm across mutants assigned to the worker.
    """
    started = time.time()
    candidate = Path(config["candidate"])
    output = Path(config["output"])
    output.mkdir(parents=True, exist_ok=True)
    prepared = output / "prepared"
    final = output / "final"
    prepared.mkdir(exist_ok=True)
    (output / "arpeggio_mutant").mkdir(exist_ok=True)
    (output / "arpeggio_parent").mkdir(exist_ok=True)
    (output / "plip_mutant").mkdir(exist_ok=True)
    runner = Runner(output / "workflow.log")
    with runner.log_path.open("a", encoding="utf-8") as log:
        log.write(f"\n[{datetime.now(timezone.utc).isoformat()}] CANDIDATE ATTEMPT\n")
    cache = Path(config["parent_cache"])
    cached = _cache_paths(cache)
    params = Path(config["params"])
    template = Path(config["ligand_template_pdb"])
    prepare_script = resource_path("scripts/prepare_ligand_complexes.py")
    cif_script = resource_path("scripts/pdb_to_arpeggio_cif.py")
    manifest_path = output / "workflow_manifest.json"
    manifest: dict[str, Any] = {
        "status": "running",
        "started_utc": datetime.now(timezone.utc).isoformat(),
        "parent": config["parent"],
        "mutant": str(candidate),
        "mutant_input_sha256": sha256(candidate.read_bytes()).hexdigest() if candidate.is_file() else None,
        "output": str(output),
        "parent_cache": str(cache),
        "persistent_pyrosetta_worker_pid": os.getpid(),
        "ligand_resname": config["ligand_resname"],
        "analysis_status": analysis_status(config, pending=True),
    }
    manifest_path.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    try:
        mutant_ligand_chain = config.get("mutant_ligand_chain") or config["ligand_chain"]
        mutant = _convert_structure_input(
            candidate,
            prepared / "mutant_input.pdb",
            config.get("mutant_chain_map", {}),
            mutant_ligand_chain,
        )
        if mutant != prepared / "mutant_input.pdb":
            shutil.copyfile(mutant, prepared / "mutant_input.pdb")
            mutant = prepared / "mutant_input.pdb"
        runner.run([
            sys.executable, str(prepare_script),
            "--reference", str(cached["analysis_input"]),
            "--complex", str(mutant),
            "--generated-ligand-pdb", str(template),
            "--params", str(params),
            "--ligand-chain", mutant_ligand_chain,
            "--output-ligand-chain", config["rosetta_ligand_chain"],
            "--output", str(prepared),
            *(["--ligand-smiles", config["ligand_smiles"]] if config.get("ligand_smiles") else []),
        ])
        mutant_prepared = prepared / f"{mutant.stem}_rosetta.pdb"

        parent_prepared = prepared / "parent_rosetta.pdb"
        _link_or_copy(cached["prepared"], parent_prepared)
        parent_artifacts = {
            "pairs": output / "parent_pair_energies.sc",
            "interface": output / "parent_interface.sc",
            "hbonds": output / "parent_hbonds.csv",
        }
        for key, destination in parent_artifacts.items():
            _link_or_copy(cached[key], destination)
        if not config.get('skip_arpeggio'):
            _link_or_copy(cached['arpeggio'], output/'arpeggio_parent'/'parent_arpeggio.json')

        parent_validation = json.loads(cached["validation"].read_text(encoding="utf-8"))["parent"]
        mutant_validation = validate_ligand_against_params(
            mutant_prepared, params, config["rosetta_ligand_chain"], template=template,
        )
        validation = {"parent": parent_validation, "mutant": mutant_validation}
        (output / "ligand_atom_validation.json").write_text(
            json.dumps(validation, indent=2) + "\n", encoding="utf-8"
        )
        if not mutant_validation["valid"]:
            raise ValueError("Mutant ligand atom validation failed")

        score_and_extract_hbonds(
            mutant_prepared,
            output / "mutant_pair_energies.sc",
            output / "mutant_interface.sc",
            output / "mutant_hbonds.csv",
            config["rosetta_ligand_chain"],
            params=[str(params)],
            weights=config["weights"],
            **({"legacy_pre_talaris": config["legacy_pre_talaris"]} if config.get("legacy_pre_talaris") is not None else {}),
        )

        if not config.get("skip_arpeggio"):
            mutant_cif = prepared / "mutant_arpeggio.cif"
            runner.run([sys.executable, str(cif_script), str(mutant_prepared), str(mutant_cif)])
            arpeggio = shutil.which("pdbe-arpeggio")
            if not arpeggio:
                raise RuntimeError("pdbe-arpeggio is not installed")
            selection = None
            if config.get("arpeggio_scope", "full") == "neighbourhood":
                selection = mutation_neighbourhood_selection(
                    parent_prepared, mutant_prepared, config["fv_chains"].split(","),
                    config["rosetta_ligand_chain"], config.get("arpeggio_neighbourhood_radius", 6.0),
                )
                filter_arpeggio_contacts(cached["arpeggio"], output / "arpeggio_parent" / "parent_arpeggio.json",
                                        selection["selected_residues"])
            manifest["arpeggio_selection"] = selection or {"scope": "full"}
            mutant_json = run_arpeggio(
                runner, arpeggio, output / "arpeggio_mutant", mutant_cif,
                **({"selectors": selection["selectors"]} if selection else {}),
            )
            if selection:
                filter_arpeggio_contacts(mutant_json, mutant_json, selection["selected_residues"])

        else:
            mutant_json = None
            manifest["arpeggio_selection"] = {"scope": "not_run"}

        if not config.get("skip_plip"):
            plip = shutil.which("plip")
            if not plip:
                raise RuntimeError("PLIP is not installed")
            runner.run([
                plip, "-f", str(mutant_prepared), "-o", str(output / "plip_mutant"),
                "-x", "-t", "--nohydro", "--name", "mutant_plip",
            ])
            parent_plip = json.loads(cached["plip"].read_text(encoding="utf-8"))["parent"]
            plip_summary = {
                "parent": parent_plip,
                "mutant": summarize_plip_xml(
                    output / "plip_mutant" / "mutant_plip.xml", config["ligand_resname"]
                ),
            }
            (output / "plip_ligand_interactions.json").write_text(
                json.dumps(plip_summary, indent=2) + "\n", encoding="utf-8"
            )

        else:
            plip_summary = None

        compare(SimpleNamespace(
            parent=str(parent_prepared), mutant=str(mutant_prepared),
            fv_chains=config["fv_chains"].split(","),
            ligand_chains=[config["rosetta_ligand_chain"]],
            scheme=config["scheme"], cdr_definition=config["cdr_definition"],
            fit_region=config["fit_region"], rmsd_atoms=config["rmsd_atoms"],
            rotamer_rmsd_cutoff=config["rotamer_rmsd_cutoff"],
            hcdr3_tip_sequence=config.get("hcdr3_tip_sequence", "YYGSS"),
            skip_arpeggio=config.get("skip_arpeggio", False), skip_plip=config.get("skip_plip", False),
            skip_network_analysis=config.get("skip_network_analysis", False),
            skip_rotamer_analysis=config.get("skip_rotamer_analysis", False),
            skip_e2_architecture=config.get("skip_e2_architecture", False),
            network_energy_presence_cutoff=config["network_energy_presence_cutoff"],
            network_energy_change_cutoff=config["network_energy_change_cutoff"],
            parent_rosetta=str(parent_artifacts["pairs"]),
            mutant_rosetta=str(output / "mutant_pair_energies.sc"),
            parent_arpeggio=str(output / "arpeggio_parent" / "parent_arpeggio.json") if not config.get("skip_arpeggio") else None,
            mutant_arpeggio=str(mutant_json) if mutant_json is not None else None,
            parent_hbonds=str(parent_artifacts["hbonds"]),
            mutant_hbonds=str(output / "mutant_hbonds.csv"),
            parent_interface=str(parent_artifacts["interface"]),
            mutant_interface=str(output / "mutant_interface.sc"), output=str(final),
        ))

        record_comparison_status(final, config)
        if not config.get("skip_e2_architecture"):
            parent_architecture = json.loads(cached["architecture"].read_text(encoding="utf-8"))["parent"]
            architecture = {
                "parent": parent_architecture,
                "mutant": analyze_estrogen_architecture(
                    mutant_prepared, mutant_json,
                    ligand_resname=config["ligand_resname"],
                    fv_chains=set(config["fv_chains"].split(",")),
                    params_path=params,
                    template_path=template,
                ),
            }
        else:
            architecture = omitted_e2_architecture()
        (output / "estrogen_binding_architecture.json").write_text(
            json.dumps(architecture, indent=2) + "\n", encoding="utf-8"
        )
        _write_automatic_summary(output, None if config.get("skip_e2_architecture") else architecture, plip_summary)
        manifest.update({
            "status": "complete",
            "completed_utc": datetime.now(timezone.utc).isoformat(),
            "elapsed_seconds": round(time.time() - started, 3),
            "params": str(params), "ligand_template_pdb": str(template), "commands": runner.commands,
            "scoring_protocol": scoring_protocol(params, config["weights"], config.get("legacy_pre_talaris")),
            "analysis_status": analysis_status(config),
            "summary": str(output / "automatic_summary.md"),
            "optimization": {
                "parent_reused": True,
                "pyrosetta_persistent_per_worker": True,
                "parent_numbering_cached_per_worker": True,
            },
        })
        status, error = "complete", ""
    except Exception as exc:
        manifest.update({
            "status": "failed", "completed_utc": datetime.now(timezone.utc).isoformat(),
            "elapsed_seconds": round(time.time() - started, 3), "error": str(exc),
            "commands": runner.commands,
        })
        status, error = "failed", str(exc)
    finally:
        manifest_path.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    return {
        "candidate": str(candidate), "result_folder": str(output),
        "status": status, "error": error, "worker_pid": os.getpid(),
    }


def run_rebased_candidate(config: dict[str, Any]) -> dict[str, Any]:
    """Recompare a saved mutant analysis against a newly cached parent.

    Mutant preparation, PyRosetta scoring, Arpeggio, PLIP, and ligand
    validation are immutable inputs here.  Only parent-relative structural,
    energetic, network, and reporting products are regenerated.
    """
    started = time.time()
    candidate = Path(config["candidate"])
    source = Path(config["source_result"])
    validate_scoring_reuse(source, Path(config['params']), config['weights'], config.get('legacy_pre_talaris'))
    validate_analysis_reuse(source, config)
    output = Path(config["output"])
    output.mkdir(parents=True, exist_ok=True)
    prepared = output / "prepared"
    final = output / "final"
    prepared.mkdir(exist_ok=True)
    (output / "arpeggio_parent").mkdir(exist_ok=True)
    (output / "arpeggio_mutant").mkdir(exist_ok=True)
    cache = Path(config["parent_cache"])
    cached = _cache_paths(cache)
    manifest_path = output / "workflow_manifest.json"
    manifest: dict[str, Any] = {
        "status": "running",
        "started_utc": datetime.now(timezone.utc).isoformat(),
        "parent": config["parent"],
        "mutant": str(candidate),
        "mutant_input_sha256": sha256(candidate.read_bytes()).hexdigest() if candidate.is_file() else None,
        "output": str(output),
        "parent_cache": str(cache),
        "mutant_artifacts_reused_from": str(source),
        "ligand_resname": config.get("ligand_resname"),
        "persistent_worker_pid": os.getpid(),
        "analysis_status": analysis_status(config, pending=True),
    }
    manifest_path.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    try:
        required_source = [
            source / "mutant_pair_energies.sc",
            source / "mutant_interface.sc",
            source / "mutant_hbonds.csv",
            source / "arpeggio_mutant" / "mutant_arpeggio.json",
            source / "plip_ligand_interactions.json",
            source / "ligand_atom_validation.json",
            source / "estrogen_binding_architecture.json",
        ]
        if config.get("skip_arpeggio"):
            required_source.remove(source / "arpeggio_mutant" / "mutant_arpeggio.json")
        if config.get("skip_plip"):
            required_source.remove(source / "plip_ligand_interactions.json")
        missing = [str(path) for path in required_source if not path.is_file()]
        if missing:
            raise FileNotFoundError("Saved mutant analysis is incomplete: " + "; ".join(missing))

        saved_mutant = _recorded_pdb(source, "mutant")
        mutant_prepared = prepared / saved_mutant.name
        _link_or_copy(saved_mutant, mutant_prepared)
        parent_prepared = prepared / cached["prepared"].name
        _link_or_copy(cached["prepared"], parent_prepared)

        artifact_links = {
            cached["pairs"]: output / "parent_pair_energies.sc",
            cached["interface"]: output / "parent_interface.sc",
            cached["hbonds"]: output / "parent_hbonds.csv",
            **({cached["arpeggio"]: output / "arpeggio_parent" / "parent_arpeggio.json"} if not config.get("skip_arpeggio") else {}),
            source / "mutant_pair_energies.sc": output / "mutant_pair_energies.sc",
            source / "mutant_interface.sc": output / "mutant_interface.sc",
            source / "mutant_hbonds.csv": output / "mutant_hbonds.csv",
            **({source / "arpeggio_mutant" / "mutant_arpeggio.json": output / "arpeggio_mutant" / "mutant_arpeggio.json"}
               if not config.get("skip_arpeggio") else {}),
        }
        for source_path, destination in artifact_links.items():
            _link_or_copy(source_path, destination)
        if not config.get("skip_arpeggio"):
            selection = None
            if config.get("arpeggio_scope", "full") == "neighbourhood":
                selection = mutation_neighbourhood_selection(
                    parent_prepared, mutant_prepared, config["fv_chains"].split(","),
                    config["rosetta_ligand_chain"], config.get("arpeggio_neighbourhood_radius", 6.0),
                )
            validate_arpeggio_coverage(source, selection)
            manifest["arpeggio_selection"] = selection or {"scope": "full"}
            if selection:
                for role in ("parent", "mutant"):
                    contacts = output / f"arpeggio_{role}" / f"{role}_arpeggio.json"
                    filter_arpeggio_contacts(contacts, contacts, selection["selected_residues"])

        else:
            manifest["arpeggio_selection"] = {"scope": "not_run"}

        source_validation = json.loads(
            (source / "ligand_atom_validation.json").read_text(encoding="utf-8")
        )
        parent_validation = json.loads(cached["validation"].read_text(encoding="utf-8"))["parent"]
        validation = {"parent": parent_validation, "mutant": source_validation["mutant"]}
        (output / "ligand_atom_validation.json").write_text(
            json.dumps(validation, indent=2) + "\n", encoding="utf-8"
        )

        compare(SimpleNamespace(
            parent=str(parent_prepared), mutant=str(mutant_prepared),
            fv_chains=config["fv_chains"].split(","),
            ligand_chains=[config["rosetta_ligand_chain"]],
            scheme=config["scheme"], cdr_definition=config["cdr_definition"],
            fit_region=config["fit_region"], rmsd_atoms=config["rmsd_atoms"],
            rotamer_rmsd_cutoff=config["rotamer_rmsd_cutoff"],
            hcdr3_tip_sequence=config.get("hcdr3_tip_sequence", "YYGSS"),
            skip_arpeggio=config.get("skip_arpeggio", False), skip_plip=config.get("skip_plip", False),
            skip_network_analysis=config.get("skip_network_analysis", False),
            skip_rotamer_analysis=config.get("skip_rotamer_analysis", False),
            skip_e2_architecture=config.get("skip_e2_architecture", False),
            network_energy_presence_cutoff=config["network_energy_presence_cutoff"],
            network_energy_change_cutoff=config["network_energy_change_cutoff"],
            parent_rosetta=str(output / "parent_pair_energies.sc"),
            mutant_rosetta=str(output / "mutant_pair_energies.sc"),
            parent_arpeggio=str(output / "arpeggio_parent" / "parent_arpeggio.json") if not config.get("skip_arpeggio") else None,
            mutant_arpeggio=str(output / "arpeggio_mutant" / "mutant_arpeggio.json") if not config.get("skip_arpeggio") else None,
            parent_hbonds=str(output / "parent_hbonds.csv"),
            mutant_hbonds=str(output / "mutant_hbonds.csv"),
            parent_interface=str(output / "parent_interface.sc"),
            mutant_interface=str(output / "mutant_interface.sc"), output=str(final),
        ))

        record_comparison_status(final, config)
        if not config.get("skip_e2_architecture"):
            source_architecture = json.loads(
                (source / "estrogen_binding_architecture.json").read_text(encoding="utf-8")
            )
            parent_architecture = json.loads(cached["architecture"].read_text(encoding="utf-8"))["parent"]
            architecture = {
                "parent": parent_architecture,
                "mutant": source_architecture["mutant"],
            }
        else:
            architecture = omitted_e2_architecture()
        (output / "estrogen_binding_architecture.json").write_text(
            json.dumps(architecture, indent=2) + "\n", encoding="utf-8"
        )

        if not config.get("skip_plip"):
            source_plip = json.loads(
                (source / "plip_ligand_interactions.json").read_text(encoding="utf-8")
            )
            parent_plip = json.loads(cached["plip"].read_text(encoding="utf-8"))["parent"]
            plip_summary = {"parent": parent_plip, "mutant": source_plip["mutant"]}
            (output / "plip_ligand_interactions.json").write_text(
                json.dumps(plip_summary, indent=2) + "\n", encoding="utf-8"
            )
        else:
            plip_summary = None

        _write_automatic_summary(output, None if config.get("skip_e2_architecture") else architecture, plip_summary)
        manifest.update({
            "status": "complete",
            "completed_utc": datetime.now(timezone.utc).isoformat(),
            "elapsed_seconds": round(time.time() - started, 3),
            "params": config['params'],
            "scoring_protocol": scoring_protocol(config["params"], config["weights"], config.get("legacy_pre_talaris")),
            "analysis_status": analysis_status(config),
            "optimization": {
                "parent_scored_once": True,
                "mutant_pyrosetta_reused": True,
                "mutant_arpeggio_reused": not config.get("skip_arpeggio", False),
                "mutant_plip_reused": not config.get("skip_plip", False),
                "comparison_numbering_cached_per_worker": True,
            },
        })
        status, error = "complete", ""
    except Exception as exc:
        manifest.update({
            "status": "failed",
            "completed_utc": datetime.now(timezone.utc).isoformat(),
            "elapsed_seconds": round(time.time() - started, 3),
            "error": str(exc),
        })
        status, error = "failed", str(exc)
    finally:
        manifest_path.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    return {
        "candidate": str(candidate), "result_folder": str(output),
        "status": status, "error": error, "worker_pid": os.getpid(),
    }
