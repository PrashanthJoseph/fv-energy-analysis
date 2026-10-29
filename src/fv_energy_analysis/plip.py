"""Read the ligand-interaction summary from a PLIP XML report."""

from __future__ import annotations

from collections import defaultdict
from pathlib import Path
import xml.etree.ElementTree as ET


INTERACTION_TAGS = {
    "hydrophobic_interactions": "hydrophobic_interaction",
    "hydrogen_bonds": "hydrogen_bond",
    "water_bridges": "water_bridge",
    "salt_bridges": "salt_bridge",
    "pi_stacks": "pi_stack",
    "pi_cation_interactions": "pi_cation_interaction",
    "halogen_bonds": "halogen_bond",
    "metal_complexes": "metal_complex",
}


def summarize_plip_xml(path: str | Path, ligand_resname: str = "E2L") -> dict[str, object]:
    """Return counts and participating protein residues for one ligand type."""
    root = ET.parse(path).getroot()
    counts: dict[str, int] = {name: 0 for name in INTERACTION_TAGS}
    residues: dict[str, set[str]] = defaultdict(set)
    binding_sites = 0

    for site in root.findall("bindingsite"):
        if (site.findtext("identifiers/hetid") or "").strip() != ligand_resname:
            continue
        binding_sites += 1
        interactions = site.find("interactions")
        if interactions is None:
            continue
        for group_name, item_name in INTERACTION_TAGS.items():
            group = interactions.find(group_name)
            items = [] if group is None else group.findall(item_name)
            counts[group_name] += len(items)
            for item in items:
                chain = (item.findtext("reschain") or "?").strip()
                number = (item.findtext("resnr") or "?").strip()
                residue_type = (item.findtext("restype") or "UNK").strip()
                residues[group_name].add(f"{chain}:{number}:{residue_type}")

    return {
        "plip_version": (root.findtext("plipversion") or "unknown").strip(),
        "ligand_resname": ligand_resname,
        "binding_site_count": binding_sites,
        "total_interaction_count": sum(counts.values()),
        "interaction_counts": counts,
        "residues_by_type": {
            name: sorted(residues.get(name, set())) for name in INTERACTION_TAGS
        },
    }
