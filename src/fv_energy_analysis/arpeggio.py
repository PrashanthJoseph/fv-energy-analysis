from __future__ import annotations

import json
import math
from collections import defaultdict
from pathlib import Path

import numpy as np

from .models import ResidueKey, canonical_pair, pair_id
from .structure import read_pdb


def mutation_neighbourhood_selection(
    parent_pdb: Path, mutant_pdb: Path, fv_chains: list[str],
    ligand_chain: str, radius: float = 6.0,
) -> dict:
    """One residue query shared by both states; keep the full input structures.

    Include mutations, whole residues within radius of any mutation in either
    state, and the complete ligand. Arpeggio expands this query to its own
    contact neighbourhood; no atoms or chemical context are cropped.
    """
    if not math.isfinite(radius) or radius <= 0:
        raise ValueError("Arpeggio neighbourhood radius must be finite and positive")
    states = [read_pdb(parent_pdb), read_pdb(mutant_pdb)]
    mutations = set()
    for chain in fv_chains:
        old, new = [state.chain_residues(chain) for state in states]
        if not old or [r.key for r in old] != [r.key for r in new]:
            raise ValueError(f"Arpeggio requires matching residue IDs and lengths in chain {chain}")
        mutations.update(a.key for a, b in zip(old, new) if a.name != b.name)
    selected = set(mutations)
    for state in states:
        ligand = state.chain_residues(ligand_chain, protein_only=False)
        if not ligand:
            raise ValueError(f"Arpeggio ligand chain {ligand_chain} is empty")
        selected.update(r.key for r in ligand)
        seed_atoms = [atom.coord for key in mutations for atom in state.by_key[key].atoms.values()
                      if atom.element.upper() not in {"H", "D"}]
        if not seed_atoms:
            continue
        seeds = np.asarray(seed_atoms)
        for residue in state.residues:
            if residue.key.chain not in fv_chains or residue.hetero:
                continue
            atoms = [a.coord for a in residue.atoms.values() if a.element.upper() not in {"H", "D"}]
            if atoms and np.any(np.sum((np.asarray(atoms)[:, None, :] - seeds[None, :, :]) ** 2, axis=2) <= radius ** 2):
                selected.add(residue.key)
    return {
        "scope": "neighbourhood", "radius_A": radius,
        "mutated_residues": [key.label() for key in sorted(mutations)],
        "selected_residues": [key.label() for key in sorted(selected)],
        "selectors": [f"/{key.chain}/{key.number}{key.insertion}/" for key in sorted(selected)],
        "selection_basis": "union of parent and mutant mutation neighbours; ligand always included",
        "full_structure_retained": True,
    }


def filter_arpeggio_contacts(source: Path, destination: Path, selected_residues: list[str]) -> Path:
    """Keep contacts touching the shared query; replace safely even for hardlinks."""
    selected = set(selected_residues)
    records = json.loads(source.read_text(encoding="utf-8"))
    records = [record for record in records
               if any(_key(record[end]).label() in selected for end in ("bgn", "end") if end in record)]
    destination.parent.mkdir(parents=True, exist_ok=True)
    temporary = destination.with_suffix(destination.suffix + ".tmp")
    temporary.write_text(json.dumps(records, indent=2) + "\n", encoding="utf-8")
    temporary.replace(destination)
    return destination


def validate_arpeggio_coverage(result_dir: Path, selection: dict | None) -> None:
    """Full outputs are reusable; a local output must cover the new query."""
    manifest = json.loads((result_dir / "workflow_manifest.json").read_text(encoding="utf-8"))
    previous = manifest.get("arpeggio_selection", {})
    if previous.get("scope", "full") == "full":
        return
    if selection is None or not set(selection["selected_residues"]).issubset(previous.get("selected_residues", [])):
        raise ValueError("Saved Arpeggio neighbourhood does not cover this comparison; rerun Arpeggio, not Rosetta")


def run_arpeggio(runner, executable: str, output_dir: Path, cif_path: Path, *, selectors=None) -> Path:
    """Query residues against an intact complex; None explicitly requests full scope."""
    output_dir.mkdir(parents=True, exist_ok=True)
    command = [executable, "-o", str(output_dir), "-m"]
    if selectors is not None:
        if not selectors:
            raise ValueError("An empty Arpeggio query must not silently select the full structure")
        selection_file = output_dir / "selection.txt"
        selection_file.write_text("\n".join(selectors) + "\n", encoding="utf-8")
        command.extend(["-sf", str(selection_file)])
    command.append(str(cif_path))
    runner.run(
        command,
        env_overrides={"PYTHONHASHSEED": "0"},
    )
    return output_dir / f"{cif_path.stem}.json"


def _key(endpoint: dict) -> ResidueKey:
    chain = str(endpoint.get("auth_asym_id") or endpoint.get("label_asym_id") or "_")
    number = int(endpoint.get("auth_seq_id") or endpoint.get("label_seq_id"))
    insertion = str(endpoint.get("pdbx_PDB_ins_code") or "").strip()
    if insertion in {"?", "."}:
        insertion = ""
    return ResidueKey(chain, number, insertion)


def _interaction_names(entry: dict) -> list[str]:
    raw = [str(x) for x in entry.get("contact", [])]
    output: list[str] = []
    for value in raw:
        upper = value.upper()
        lower = value.lower()
        if upper in {"FT", "ET", "FE", "EF", "FF", "EE"} or entry.get("type") == "plane-plane":
            name = "pi_stack"
        elif upper == "CATIONPI":
            name = "cation_pi"
        elif upper.endswith("PI"):
            name = upper.lower().replace("pi", "_pi").strip("_")
        elif lower == "aromatic":
            name = "aromatic_contact"
        else:
            name = lower
        if name not in output:
            output.append(name)
    return output


def read_arpeggio(path: str | Path, fv_chains: set[str] | None = None) -> list[dict[str, object]]:
    data = json.loads(Path(path).read_text(encoding="utf-8"))
    if isinstance(data, dict):
        for candidate in ("interactions", "contacts", "data"):
            if isinstance(data.get(candidate), list):
                data = data[candidate]
                break
    if not isinstance(data, list):
        raise ValueError(f"Expected an Arpeggio JSON list in {path}")
    aggregate: dict[tuple[ResidueKey, ResidueKey], dict[str, object]] = {}
    counts: dict[tuple[ResidueKey, ResidueKey], dict[str, int]] = defaultdict(lambda: defaultdict(int))
    distances: dict[tuple[ResidueKey, ResidueKey], dict[str, list[float]]] = defaultdict(lambda: defaultdict(list))
    for entry in data:
        if not isinstance(entry, dict) or "bgn" not in entry or "end" not in entry:
            continue
        try:
            a, b = _key(entry["bgn"]), _key(entry["end"])
        except (TypeError, ValueError):
            continue
        if a == b or (fv_chains and a.chain not in fv_chains and b.chain not in fv_chains):
            continue
        pair = canonical_pair(a, b)
        if pair not in aggregate:
            endpoints = {a: entry["bgn"], b: entry["end"]}
            aggregate[pair] = {
                "pair_id": pair_id(*pair),
                "residue1": pair[0].label(),
                "residue1_name": endpoints[pair[0]].get("label_comp_id", ""),
                "residue2": pair[1].label(),
                "residue2_name": endpoints[pair[1]].get("label_comp_id", ""),
            }
        for name in _interaction_names(entry):
            counts[pair][name] += 1
            if isinstance(entry.get("distance"), (int, float)):
                distances[pair][name].append(float(entry["distance"]))
    rows: list[dict[str, object]] = []
    for pair, base in aggregate.items():
        row = dict(base)
        row["interaction_types"] = ";".join(sorted(counts[pair]))
        row["total_contact_count"] = sum(counts[pair].values())
        for name, count in counts[pair].items():
            row[f"count_{name}"] = count
            if distances[pair][name]:
                row[f"min_distance_{name}_A"] = min(distances[pair][name])
        rows.append(row)
    return rows
