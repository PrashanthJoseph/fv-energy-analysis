"""Transfer original ModelCIF pLDDT to canonical heavy-chain residues 97–113."""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
import re
import shutil
from statistics import fmean

from .fingerprint_counts import _read_csv, _windows_path
from .report import write_csv, write_json
from .resources import resource_path
from .structure import THREE_TO_ONE, read_pdb, residue_sequence


RESIDUES = tuple(range(97, 114))
SECTION = "HCDR3 pLDDT"


def _matching_cifs(candidate: str, sources: dict[str, list[Path]]) -> list[Path]:
    """Resolve exact model names, including a strictly duplicated job wrapper.

    Some relaxed exports are named ``job__job_model_N_relaxed.pdb`` while
    their original CIF is ``job_model_N.cif``. Only that repeated-job form
    is unwrapped; arbitrary prefixes and different model numbers never match.
    All equivalent matches are retained so the caller can reject ambiguity.
    """
    stem = Path(candidate).stem.removesuffix("_relaxed").casefold()
    stems = [stem]
    repeated = re.fullmatch(r"(?P<job>.+)__(?P=job)_model_(?P<model>\d+)", stem)
    if repeated:
        stems.append(f"{repeated['job']}_model_{repeated['model']}")
    return list(dict.fromkeys(path for key in stems for path in sources.get(key, [])))


def confidence_columns() -> list[dict[str, object]]:
    fields = [("hcdr3_plddt_mean", "Mean pLDDT"), ("hcdr3_plddt_min", "Minimum pLDDT")]
    fields += [(f"hcdr3_plddt_{number}", f"H{number} pLDDT") for number in RESIDUES]
    return [{"key": key, "label": label, "section": SECTION, "kind": "number"} for key, label in fields]


def extract_hcdr3_plddt(cif_path: Path, normalized_pdb: Path) -> dict[str, object]:
    """Match complete Fv sequences before transferring heavy-chain local pLDDT.

    Numbers 97–113 refer to the normalized analysis PDB, not an independently
    assigned antibody-numbering scheme. Arbitrary crystallographic B factors
    are never interpreted as pLDDT. If canonical L is present, its complete
    sequence must also match uniquely, preventing cross-variant confidence
    transfer when several light-chain mutants share an identical H sequence.
    """
    try:
        import gemmi
    except ImportError as error:
        raise RuntimeError("CIF confidence extraction requires gemmi; install fv-energy-analysis[reporting].") from error

    block = gemmi.cif.read_file(str(cif_path)).sole_block()
    metrics = block.find(["_ma_qa_metric.id", "_ma_qa_metric.type", "_ma_qa_metric.mode"])
    metric_ids = [row[0] for row in metrics if row[1].lower() == "plddt" and row[2].lower() == "local"]
    if len(metric_ids) != 1:
        raise ValueError(f"Expected one explicitly declared local pLDDT metric in {cif_path}")

    atoms = block.find([
        "_atom_site.label_atom_id", "_atom_site.label_comp_id", "_atom_site.label_asym_id",
        "_atom_site.label_seq_id", "_atom_site.pdbx_PDB_model_num",
    ])
    if not len(atoms):
        raise ValueError(f"No labeled atom-site records in {cif_path}")
    model_id = atoms[0][4]
    chains: dict[str, dict[int, str]] = {}
    for atom, name, chain, seq_id, model in atoms:
        if model == model_id and atom == "CA" and name in THREE_TO_ONE and seq_id not in (".", "?"):
            chains.setdefault(chain, {})[int(seq_id)] = name
    target_structure = read_pdb(normalized_pdb)
    target = target_structure.chain_residues("H")
    target_sequence = residue_sequence(target)
    if not target_sequence:
        raise ValueError(f"Canonical chain H is absent from {normalized_pdb}")
    matches = [
        (chain, sorted(residues.items())) for chain, residues in chains.items()
        if "".join(THREE_TO_ONE[name] for _, name in sorted(residues.items())) == target_sequence
    ]
    if len(matches) != 1:
        raise ValueError(f"Expected one exact full-chain sequence match to canonical H in {cif_path}; found {len(matches)}")
    source_chain, source_residues = matches[0]
    light = target_structure.chain_residues("L")
    if light:
        light_sequence = residue_sequence(light)
        light_matches = [
            chain for chain, residues in chains.items()
            if chain != source_chain
            and "".join(THREE_TO_ONE[name] for _, name in sorted(residues.items())) == light_sequence
        ]
        if len(light_matches) != 1:
            raise ValueError(
                f"Expected one exact full-chain sequence match to canonical L in {cif_path}; "
                f"found {len(light_matches)}. Both Fv sequences must match before confidence transfer."
            )
    local = block.find([
        "_ma_qa_metric_local.model_id", "_ma_qa_metric_local.label_asym_id",
        "_ma_qa_metric_local.label_seq_id", "_ma_qa_metric_local.metric_id",
        "_ma_qa_metric_local.metric_value",
    ])
    scores: dict[int, float] = {}
    for model, chain, seq_id, metric, value in local:
        if model == model_id and chain == source_chain and metric == metric_ids[0]:
            number, score = int(seq_id), float(value)
            if not math.isfinite(score) or not 0 <= score <= 100:
                raise ValueError(f"Invalid pLDDT {value} for {chain}:{seq_id} in {cif_path}")
            if number in scores:
                raise ValueError(f"Duplicate pLDDT for {chain}:{seq_id} in {cif_path}")
            scores[number] = score
    transferred = []
    for residue, (seq_id, _) in zip(target, source_residues):
        if residue.key.number not in RESIDUES:
            continue
        if residue.key.insertion or seq_id not in scores:
            raise ValueError(f"Missing or ambiguous confidence for {residue.key.label()} in {cif_path}")
        transferred.append({
            "residue": residue.key.label(), "residue_number": residue.key.number,
            "residue_name": residue.name, "source_chain": source_chain,
            "source_seq_id": seq_id, "plddt": scores[seq_id],
        })
    if tuple(row["residue_number"] for row in transferred) != RESIDUES:
        raise ValueError(f"Incomplete canonical H97–H113 interval in {normalized_pdb}")
    values = {"hcdr3_plddt_mean": fmean(row["plddt"] for row in transferred),
              "hcdr3_plddt_min": min(row["plddt"] for row in transferred)}
    values.update({f"hcdr3_plddt_{row['residue_number']}": row["plddt"] for row in transferred})
    return {"source_cif": str(cif_path.resolve()), "source_cif_chain": source_chain,
            "canonical_chain": "H", "score_source": "_ma_qa_metric_local.metric_value",
            "residues": transferred, "values": values}


def _write_detail(result_dir: Path, confidence: dict[str, object]) -> None:
    write_json(result_dir / "hcdr3_plddt.json", confidence)
    write_csv(result_dir / "hcdr3_plddt.csv", confidence["residues"])
    summary_path = result_dir / "automatic_summary.json"
    summary = json.loads(summary_path.read_text(encoding="utf-8"))
    summary["hcdr3_confidence"] = confidence
    write_json(summary_path, summary)
    lines = ["## HCDR3 prediction confidence (H97–H113)", "",
             f"Original CIF: `{confidence['source_cif']}`. Source chain `{confidence['source_cif_chain']}` "
             "matches the complete canonical H sequence. Scores are original prediction pLDDT (0–100).", "",
             "| Canonical residue | Amino acid | CIF residue | pLDDT |", "|---|---|---|---:|"]
    lines += [f"| {r['residue']} | {r['residue_name']} | {r['source_chain']}:{r['source_seq_id']} | {r['plddt']:.3f} |"
              for r in confidence["residues"]]
    values = confidence["values"]
    lines += ["", f"Mean: {values['hcdr3_plddt_mean']:.3f}; minimum: {values['hcdr3_plddt_min']:.3f}."]
    path = result_dir / "automatic_summary.md"
    text = path.read_text(encoding="utf-8")
    begin, end = "<!-- hcdr3-confidence:begin -->", "<!-- hcdr3-confidence:end -->"
    block = begin + "\n" + "\n".join(lines) + "\n" + end
    if begin in text and end in text:
        prefix, remainder = text.split(begin, 1)
        _, suffix = remainder.split(end, 1)
        text = prefix.rstrip() + "\n\n" + block + suffix
    else:
        text = text.rstrip() + "\n\n" + block + "\n"
    path.write_text(text, encoding="utf-8")


def add_confidence(batch: Path, cif_folder: Path, parent_cif: Path | None = None) -> Path:
    """Enrich saved reports only; existing energies and fingerprints are reused."""
    batch, cif_folder = Path(batch).resolve(), Path(cif_folder).resolve()
    parent_cif = Path(parent_cif).resolve() if parent_cif else None
    payload_path = batch / "batch_summary_grouped_fingerprints.json"
    payload = json.loads(payload_path.read_text(encoding="utf-8"))
    if not cif_folder.is_dir():
        raise ValueError(f"CIF folder does not exist: {cif_folder}")
    sources: dict[str, list[Path]] = {}
    for path in cif_folder.rglob("*"):
        if path.suffix.lower() in (".cif", ".mmcif"):
            sources.setdefault(path.stem.casefold(), []).append(path)
    bases = {r["candidate"]: r for r in _read_csv(batch / "batch_summary.csv")}
    extracted = []
    for row in payload["rows"]:
        candidate = row["candidate"]
        matches = _matching_cifs(candidate, sources)
        if len(matches) != 1:
            raise ValueError(f"Expected one exact matching CIF for {candidate}; found {len(matches)}")
        result_dir = _windows_path(bases[candidate]["result_folder"])
        normalized = result_dir / "final" / "mutant_aligned_to_parent.pdb"
        confidence = extract_hcdr3_plddt(matches[0], normalized)
        extracted.append((result_dir, confidence))
        row.update(confidence["values"])
    parent_values = {column["key"]: None for column in confidence_columns()}
    if parent_cif:
        manifest = json.loads((batch / "batch_manifest.json").read_text(encoding="utf-8"))
        cache = _windows_path(manifest["optimization"]["parent_cache"])
        recorded = json.loads((cache / "final" / "summary.json").read_text(encoding="utf-8"))["parent"]
        prepared = _windows_path(recorded)
        if not prepared.is_file():
            prepared = cache / "prepared" / prepared.name
        parent = extract_hcdr3_plddt(Path(parent_cif), prepared)
        parent_values = parent["values"]
        write_json(batch / "parent_hcdr3_plddt.json", parent)
    payload["parent"].update(parent_values)
    payload["columns"] = [c for c in payload["columns"] if c["section"] != SECTION] + confidence_columns()
    payload["hcdr3_confidence"] = {"cif_folder": str(cif_folder), "parent_cif": str(parent_cif) if parent_cif else None,
                                  "matched_candidates": len(extracted), "residues_per_candidate": len(RESIDUES)}
    for result_dir, confidence in extracted:
        _write_detail(result_dir, confidence)
    write_json(payload_path, payload)
    csv_path = batch / "batch_summary_grouped_fingerprints.csv"
    try:
        write_csv(csv_path, payload["rows"])
    except PermissionError:
        csv_path = batch / "batch_summary_grouped_fingerprints_with_plddt.csv"
        write_csv(csv_path, payload["rows"])
    numeric_path = batch / "batch_interaction_fingerprint_numeric.csv"
    if numeric_path.is_file():
        values_by_name = {r["candidate"]: r for r in payload["rows"]}
        numeric = _read_csv(numeric_path)
        for row in numeric:
            row.update({c["key"]: values_by_name[row["candidate"]][c["key"]] for c in confidence_columns()})
        try:
            write_csv(numeric_path, numeric)
        except PermissionError:
            numeric_path = batch / "batch_interaction_fingerprint_numeric_with_plddt.csv"
            write_csv(numeric_path, numeric)
    manifest_path = batch / "batch_manifest.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    manifest["hcdr3_confidence"] = {**payload["hcdr3_confidence"], "summary_csv": str(csv_path),
                                   "numeric_csv": str(numeric_path) if numeric_path.is_file() else None}
    manifest.setdefault("interaction_fingerprints", {})["grouped_summary_csv"] = str(csv_path)
    if numeric_path.is_file():
        manifest["interaction_fingerprints"]["numeric_flat_csv"] = str(numeric_path)
    write_json(manifest_path, manifest)
    shutil.copyfile(resource_path("docs/interaction_fingerprint_columns.md"), batch / "INTERACTION_FINGERPRINT_COLUMNS.md")
    return csv_path


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--batch-results", required=True)
    parser.add_argument("--cif-folder", required=True)
    parser.add_argument("--parent-cif")
    args = parser.parse_args(argv)
    result = add_confidence(Path(args.batch_results), Path(args.cif_folder), Path(args.parent_cif) if args.parent_cif else None)
    print(f"HCDR3 pLDDT added: {result}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
