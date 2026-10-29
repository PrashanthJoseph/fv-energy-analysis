"""Batch BIOVIA Discovery Studio ligand-interaction analysis on saved results."""

from __future__ import annotations

import argparse
import csv
import json
import os
import re
import shutil
import subprocess
from collections import Counter, defaultdict
from pathlib import Path

from .report import write_csv, write_json
from .resources import resource_path
from .paths import native_path
from .fingerprint_counts import _ligand_resname
from .fingerprint_reports import generate as generate_fingerprints


DEFAULT_DS_ROOT = Path(r"C:\Program Files\BIOVIA\Discovery Studio 2024")
TYPE_NAMES = {
    "conventionalHBondType": "conventional_hbond",
    "carbonHBondType": "carbon_hbond",
    "waterBridgeHBondType": "water_bridge_hbond",
    "waterHBondType": "water_hbond",
    "saltBridgeType": "salt_bridge",
    "attractiveChargeType": "attractive_charge",
    "piCationType": "pi_cation",
    "piAnionType": "pi_anion",
    "piPiStackedType": "pi_pi_stacked",
    "piPiTShapedType": "pi_pi_t_shaped",
    "piAlkylType": "pi_alkyl",
    "piSigmaType": "pi_sigma",
    "piDonorType": "donor_pi",
    "piLonePairType": "lone_pair_pi",
    "piSulfurType": "pi_sulfur",
    "piAmideType": "pi_amide",
    "alkylType": "alkyl",
    "sulfurType": "sulfur",
    "halogenClBrIType": "halogen_cl_br_i",
    "halogenFluorineType": "halogen_f",
    "metalAcceptorType": "metal_acceptor",
    "unfavorableVdwBumpType": "unfavorable_vdw_bump",
    "unfavorableDonorDonorType": "unfavorable_donor_donor",
    "unfavorableAcceptorAcceptorType": "unfavorable_acceptor_acceptor",
    "unfavorablePositivePositiveType": "unfavorable_positive_positive",
    "unfavorableNegativeNegativeType": "unfavorable_negative_negative",
    "unfavorableMetalDonorType": "unfavorable_metal_donor",
}


def _windows_path(value: str | Path) -> Path:
    return native_path(value)


def _read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8-sig") as handle:
        return list(csv.DictReader(handle))


def _prepared_pdb(result_dir: Path, parent: bool = False) -> Path:
    prepared = result_dir / "prepared"
    summary = result_dir / 'final' / 'summary.json'
    recorded = _windows_path(json.loads(summary.read_text(encoding='utf-8'))['parent' if parent else 'mutant'])
    path = recorded if recorded.is_file() else prepared / recorded.name
    if not path.is_file():
        raise FileNotFoundError(f"Recorded prepared PDB is missing: {recorded}")
    return path


def _copy_to_short_stage(source: Path, destination: Path) -> None:
    """Copy through the Windows extended-path namespace when needed."""
    source_text = str(source.resolve())
    source_for_copy = Path("\\\\?\\" + source_text) if len(source_text) >= 248 else source
    shutil.copyfile(source_for_copy, destination)


def _protein_and_ligand(record: dict[str, str]) -> tuple[str, str]:
    if record.get("from_role") == "receptor" and record.get("to_role") == "ligand":
        return record["from_atoms"], record["to_atoms"]
    if record.get("to_role") == "receptor" and record.get("from_role") == "ligand":
        return record["to_atoms"], record["from_atoms"]
    return record.get("from_atoms", ""), record.get("to_atoms", "")


def _protein_residue(atom_group: str) -> str:
    first = atom_group.split(";", 1)[0]
    match = re.match(r"^([^:]+):([A-Za-z0-9]{3})(-?\d+[A-Za-z]?):", first)
    if not match:
        return first or "unknown"
    return f"{match.group(1)}:{match.group(3)}:{match.group(2).upper()}"


def _atom_group_key(group: str) -> str:
    """Match corresponding named atoms, independent of residue name/group order."""
    return ";".join(sorted(re.sub(r"^([^:]+):[A-Za-z0-9]{3}(-?\d+[A-Za-z]?):", r"\1:\2:", atom.strip())
                           for atom in group.split(";") if atom.strip()))


def _normalize(records: list[dict[str, str]]) -> list[dict[str, object]]:
    normalized: list[dict[str, object]] = []
    for record in records:
        protein_atoms, ligand_atoms = _protein_and_ligand(record)
        interaction_type = TYPE_NAMES.get(record["interaction_type"], record["interaction_type"])
        normalized.append({
            **record,
            "interaction_type_raw": record["interaction_type"],
            "interaction_type": interaction_type,
            "distance_A": float(record["distance_A"]),
            "protein_residue": _protein_residue(protein_atoms),
            "protein_atoms": protein_atoms,
            "ligand_atoms": ligand_atoms,
            "atom_group_pair_id": f"{_atom_group_key(protein_atoms)}<->{_atom_group_key(ligand_atoms)}",
            "residue_pair_id": f"{_protein_residue(protein_atoms)}<->{_protein_residue(ligand_atoms)}",
        })
    return normalized


def _summary(records: list[dict[str, object]], prefix: str) -> dict[str, object]:
    counts = Counter(str(row["interaction_type"]) for row in records)
    residues: dict[str, set[str]] = defaultdict(set)
    pairs: dict[str, set[str]] = defaultdict(set)
    for row in records:
        interaction_type = str(row["interaction_type"])
        residues[interaction_type].add(str(row["protein_residue"]))
        pairs[interaction_type].add(str(row["atom_group_pair_id"]))
    output: dict[str, object] = {f"{prefix}_ds_total_interaction_annotations": len(records)}
    for interaction_type in sorted(set(TYPE_NAMES.values()) | set(counts)):
        output[f"{prefix}_ds_{interaction_type}_annotations"] = counts.get(interaction_type, 0)
        output[f"{prefix}_ds_{interaction_type}_unique_residues"] = len(residues.get(interaction_type, set()))
        output[f"{prefix}_ds_{interaction_type}_unique_atom_group_pairs"] = len(pairs.get(interaction_type, set()))
    return output


def _changes(
    candidate: str,
    result_folder: str,
    parent: list[dict[str, object]],
    mutant: list[dict[str, object]],
) -> list[dict[str, object]]:
    def keyed(rows: list[dict[str, object]]) -> dict[tuple[str, str, str], dict[str, object]]:
        return {
            (str(row["interaction_type"]), _atom_group_key(str(row["protein_atoms"])), _atom_group_key(str(row["ligand_atoms"]))): row
            for row in rows
        }
    parent_by, mutant_by = keyed(parent), keyed(mutant)
    output: list[dict[str, object]] = []
    for event, keys in (
        ("formed", sorted(set(mutant_by) - set(parent_by))),
        ("broken", sorted(set(parent_by) - set(mutant_by))),
    ):
        source = mutant_by if event == "formed" else parent_by
        for key in keys:
            row = source[key]
            output.append({
                "candidate": candidate,
                "result_folder": result_folder,
                "event": event,
                "interaction_type": row["interaction_type"],
                "interaction_name": row["interaction_name"],
                "protein_residue": row["protein_residue"],
                "protein_atoms": row["protein_atoms"],
                "ligand_atoms": row["ligand_atoms"],
                "distance_A": row["distance_A"],
                "atom_group_pair_id": row["atom_group_pair_id"],
                "note": "Discovery Studio geometry only; no interaction energy is assigned",
            })
    return output


def _arpeggio_cation_pi(result_dir: Path, ligand_resname: str) -> bool:
    path = result_dir / "arpeggio_mutant" / "mutant_arpeggio.json"
    if not path.is_file():
        return False
    for row in json.loads(path.read_text(encoding="utf-8")):
        contacts = {str(value).upper() for value in row.get("contact", [])}
        if "CATIONPI" not in contacts:
            continue
        endpoints = (row.get("bgn", {}), row.get("end", {}))
        if any(str(endpoint.get("label_comp_id", "")).upper() == ligand_resname.upper() for endpoint in endpoints):
            return True
    return False


def _positive_detectors(*pairs: tuple[str, bool]) -> tuple[str, int]:
    names = [name for name, present in pairs if present]
    return ";".join(names) or "none", len(names)


def run(args: argparse.Namespace) -> Path:
    if os.name != "nt":
        raise RuntimeError("Discovery Studio post-processing must run in native Windows Python")
    batch = Path(args.batch_results).resolve()
    manifest_path = batch / "batch_manifest.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    if getattr(args, "reports_only", False):
        statuses = _read_csv(batch / "discovery_studio_status.csv")
        completed = [item for item in manifest["candidates"] if item["status"] == "complete"]
        if len(statuses) != len(completed) + 1 or any(row["status"] != "complete" for row in statuses):
            raise ValueError("Saved Discovery Studio results are incomplete; refusing a reports-only resume")
        for item in completed:
            if not (_windows_path(item["result_folder"]) / "discovery_studio_interactions.csv").is_file():
                raise FileNotFoundError(f"Saved Discovery Studio contacts missing: {item['candidate']}")
        if not (batch / "discovery_studio_parent_interactions.csv").is_file():
            raise FileNotFoundError("Saved parent Discovery Studio contacts missing")
        ds_by_candidate = {row["candidate"]: row for row in _read_csv(batch / "discovery_studio_summary.csv")}
        merged_path = batch / "batch_summary_with_discovery_studio.csv"
        write_csv(merged_path, [{**row, **ds_by_candidate.get(row["candidate"], {})}
                              for row in _read_csv(batch / "batch_summary.csv")])
        fingerprint_path = None if getattr(args, "skip_fingerprint_reports", False) else generate_fingerprints(batch)
        manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
        manifest["discovery_studio"] = {"status": "complete", "saved_contacts_reused": True,
            "completed_structures": len(statuses), "summary_csv": str(batch / "discovery_studio_summary.csv"),
            "merged_batch_summary_csv": str(merged_path), "add_hydrogens": False,
            "interaction_fingerprint_summary_csv": str(fingerprint_path) if fingerprint_path else None}
        manifest_path.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
        return merged_path
    candidates = [item for item in manifest.get("candidates", []) if item.get("status") == "complete"]
    if not candidates:
        summary_path = _windows_path(manifest.get("summary_csv", batch / "batch_summary.csv"))
        if not summary_path.is_file():
            summary_path = batch / "batch_summary.csv"
        candidates = [
            {
                "candidate": row.get("candidate_path") or row.get("candidate", ""),
                "result_folder": row.get("result_folder", ""),
                "status": row.get("status", ""),
                "recovered_from_summary_csv": True,
            }
            for row in _read_csv(summary_path)
            if row.get("status") == "complete"
        ]
    if not candidates:
        raise RuntimeError("No completed candidate results were found")

    ds_root = Path(args.discovery_studio_root)
    perl_bat = ds_root / "bin" / "perl.bat"
    script = resource_path("scripts/discovery_studio_batch_interactions.pl")
    if not perl_bat.is_file():
        raise FileNotFoundError(f"Discovery Studio Perl wrapper not found: {perl_bat}")
    if not script.is_file():
        raise FileNotFoundError(f"Discovery Studio batch script not found: {script}")

    parent_cache = _windows_path(manifest.get("optimization", {}).get("parent_cache", ""))
    parent_pdb = _prepared_pdb(parent_cache, parent=True)
    # tempfile.mkdtemp can create app-container ACLs when invoked by the
    # Microsoft Store Python build.  An ordinary project-local directory is
    # both short enough for Discovery Studio's legacy path handling and
    # writable by the invoking Windows user.
    stage = Path.cwd() / "tmp" / "discovery_studio_batch" / f"run_{os.getpid()}"
    stage.mkdir(parents=True, exist_ok=False)
    raw_csv, status_csv = stage / "raw.csv", stage / "status.csv"
    ds_manifest = stage / "manifest.tsv"
    entries: list[tuple[str, str, str, str, Path]] = []
    staged_parent = stage / "parent.pdb"
    _copy_to_short_stage(parent_pdb, staged_parent)
    entries.append(("parent", "parent", "parent", str(parent_cache), staged_parent))
    for index, item in enumerate(candidates, 1):
        result_dir = _windows_path(item["result_folder"])
        candidate = Path(str(item["candidate"])).name
        source = _prepared_pdb(result_dir)
        staged = stage / f"m{index:04d}.pdb"
        _copy_to_short_stage(source, staged)
        entries.append((f"m{index:04d}", "mutant", candidate, str(result_dir), staged))
    with ds_manifest.open("w", encoding="utf-8", newline="") as handle:
        handle.write("structure_id\trole\tcandidate\tresult_folder\tpdb_path\n")
        for entry in entries:
            handle.write("\t".join(str(value) for value in entry) + "\n")

    command = [
        "cmd.exe", "/d", "/c", str(perl_bat), str(script),
        str(ds_manifest), str(raw_csv), str(status_csv),
    ]
    completed = subprocess.run(command, text=True, capture_output=True, check=False)
    (batch / "discovery_studio_batch.log").write_text(
        completed.stdout + completed.stderr, encoding="utf-8"
    )
    if completed.returncode:
        raise RuntimeError(f"Discovery Studio batch failed; see {batch / 'discovery_studio_batch.log'}")

    statuses = _read_csv(status_csv)
    failures = [row for row in statuses if row["status"] != "complete"]
    if failures:
        write_csv(batch / "discovery_studio_failures.csv", failures)
        raise RuntimeError(f"Discovery Studio failed for {len(failures)} structures")
    raw_records = _read_csv(raw_csv)
    normalized = _normalize(raw_records)
    by_id: dict[str, list[dict[str, object]]] = defaultdict(list)
    for record in normalized:
        by_id[str(record["structure_id"])].append(record)
    parent_records = by_id["parent"]
    write_csv(batch / "discovery_studio_parent_interactions.csv", parent_records)
    write_json(batch / "discovery_studio_parent_interactions.json", parent_records)

    master_rows: list[dict[str, object]] = []
    all_changes: list[dict[str, object]] = []
    for index, item in enumerate(candidates, 1):
        structure_id = f"m{index:04d}"
        candidate = Path(str(item["candidate"])).name
        result_dir = _windows_path(item["result_folder"])
        records = by_id[structure_id]
        changes = _changes(candidate, str(result_dir), parent_records, records)
        write_csv(result_dir / "discovery_studio_interactions.csv", records)
        write_json(result_dir / "discovery_studio_interactions.json", records)
        write_csv(result_dir / "discovery_studio_interaction_changes.csv", changes)
        all_changes.extend(changes)

        row: dict[str, object] = {
            "candidate": candidate,
            "result_folder": str(result_dir),
            **_summary(parent_records, "parent"),
            **_summary(records, "mutant"),
            "ds_formed_atom_group_interactions": sum(change["event"] == "formed" for change in changes),
            "ds_broken_atom_group_interactions": sum(change["event"] == "broken" for change in changes),
        }
        for interaction_type in sorted(set(TYPE_NAMES.values())):
            parent_count = int(row[f"parent_ds_{interaction_type}_annotations"])
            mutant_count = int(row[f"mutant_ds_{interaction_type}_annotations"])
            row[f"delta_ds_{interaction_type}_annotations"] = mutant_count - parent_count

        automatic = json.loads((result_dir / "automatic_summary.json").read_text(encoding="utf-8"))
        plip = automatic.get("plip_ligand_interactions", {}).get("mutant", {}).get("interaction_counts", {})
        mutant_ds_pi_cation = int(row["mutant_ds_pi_cation_annotations"]) > 0
        mutant_arpeggio_pi_cation = _arpeggio_cation_pi(result_dir, _ligand_resname(automatic, result_dir))
        mutant_plip_pi_cation = int(plip.get("pi_cation_interactions", 0)) > 0
        names, count = _positive_detectors(
            ("DiscoveryStudio", mutant_ds_pi_cation),
            ("Arpeggio", mutant_arpeggio_pi_cation),
            ("PLIP", mutant_plip_pi_cation),
        )
        row.update({
            "mutant_ds_pi_cation_detected": mutant_ds_pi_cation,
            "mutant_arpeggio_pi_cation_detected": mutant_arpeggio_pi_cation,
            "mutant_plip_pi_cation_detected": mutant_plip_pi_cation,
            "mutant_pi_cation_positive_detectors": names,
            "mutant_pi_cation_detector_count": count,
            "mutant_pi_cation_consensus_2plus": count >= 2,
            "mutant_ds_pi_anion_detected": int(row["mutant_ds_pi_anion_annotations"]) > 0,
        })
        master_rows.append(row)

    write_csv(batch / "discovery_studio_summary.csv", master_rows)
    write_csv(batch / "discovery_studio_interaction_changes.csv", all_changes)
    write_csv(batch / "discovery_studio_status.csv", statuses)
    definitions = resource_path("docs/discovery_studio_columns.md")
    shutil.copyfile(definitions, batch / "DISCOVERY_STUDIO_COLUMNS.md")

    summary_path = _windows_path(manifest.get("summary_csv", batch / "batch_summary.csv"))
    base_rows = _read_csv(summary_path)
    ds_by_candidate = {str(row["candidate"]): row for row in master_rows}
    merged = [{**row, **ds_by_candidate.get(row.get("candidate", ""), {})} for row in base_rows]
    merged_path = batch / "batch_summary_with_discovery_studio.csv"
    write_csv(merged_path, merged)

    fingerprint_path = None if getattr(args, "skip_fingerprint_reports", False) else generate_fingerprints(batch)
    # Fingerprint generation may add configuration/provenance. Do not overwrite
    # those changes using the pre-generation manifest snapshot.
    manifest = json.loads(manifest_path.read_text(encoding='utf-8'))

    manifest["discovery_studio"] = {
        "status": "complete",
        "version": None,  # Release version is not detected; do not infer it from the folder name.
        "installation_root": str(ds_root.resolve()),
        "build_info": (ds_root / "buildinfo.txt").read_text(encoding="utf-8").strip() if (ds_root / "buildinfo.txt").is_file() else None,
        "add_hydrogens": False,
        "completed_structures": len(statuses),
        "summary_csv": str(batch / "discovery_studio_summary.csv"),
        "interaction_changes_csv": str(batch / "discovery_studio_interaction_changes.csv"),
        "merged_batch_summary_csv": str(merged_path),
        "interaction_fingerprint_summary_csv": str(fingerprint_path) if fingerprint_path else None,
    }
    manifest_path.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    if not args.keep_stage:
        shutil.rmtree(stage)
    return merged_path


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--batch-results", required=True)
    parser.add_argument("--reports-only", action="store_true", help="Reuse and verify saved detector outputs; do not rerun Discovery Studio")
    parser.add_argument("--skip-fingerprint-reports", action="store_true", help="Leave fingerprint generation to a following reporting stage")
    parser.add_argument("--discovery-studio-root", default=str(DEFAULT_DS_ROOT))
    parser.add_argument("--keep-stage", action="store_true")
    return parser


def main(argv: list[str] | None = None) -> int:
    output = run(build_parser().parse_args(argv))
    print(f"Discovery Studio batch analysis complete: {output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
