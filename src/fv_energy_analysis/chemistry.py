from __future__ import annotations


PRIMARY_CLASS = {
    "A": "hydrophobic_aliphatic", "V": "hydrophobic_aliphatic", "L": "hydrophobic_aliphatic",
    "I": "hydrophobic_aliphatic", "M": "hydrophobic_aliphatic",
    "F": "aromatic", "W": "aromatic", "Y": "aromatic",
    "S": "polar_uncharged", "T": "polar_uncharged", "N": "polar_uncharged", "Q": "polar_uncharged",
    "D": "acidic_negative", "E": "acidic_negative",
    "K": "basic_positive", "R": "basic_positive", "H": "basic_or_titratable",
    "C": "special_thiol", "G": "special_glycine", "P": "special_proline",
}

PROPERTIES = {
    "A": {"hydrophobic", "small"}, "R": {"basic", "positive", "polar", "long"},
    "N": {"polar", "amide"}, "D": {"acidic", "negative", "polar"},
    "C": {"polarizable", "thiol", "small"}, "Q": {"polar", "amide", "long"},
    "E": {"acidic", "negative", "polar", "long"}, "G": {"small", "flexible"},
    "H": {"aromatic", "basic", "titratable", "polar"}, "I": {"hydrophobic", "branched"},
    "L": {"hydrophobic", "branched"}, "K": {"basic", "positive", "polar", "long"},
    "M": {"hydrophobic", "sulfur", "long"}, "F": {"aromatic", "hydrophobic"},
    "P": {"hydrophobic", "rigid"}, "S": {"polar", "hydroxyl", "small"},
    "T": {"polar", "hydroxyl", "branched"}, "W": {"aromatic", "hydrophobic", "large"},
    "Y": {"aromatic", "hydroxyl", "polar"}, "V": {"hydrophobic", "branched", "small"},
}


def substitution_chemistry(parent: str, mutant: str) -> dict[str, object]:
    p_class = PRIMARY_CLASS[parent]
    m_class = PRIMARY_CLASS[mutant]
    shared = sorted(PROPERTIES[parent] & PROPERTIES[mutant])
    lost = sorted(PROPERTIES[parent] - PROPERTIES[mutant])
    gained = sorted(PROPERTIES[mutant] - PROPERTIES[parent])
    return {
        "parent_chemical_class": p_class,
        "mutant_chemical_class": m_class,
        "same_primary_chemical_class": p_class == m_class,
        "parent_properties": ";".join(sorted(PROPERTIES[parent])),
        "mutant_properties": ";".join(sorted(PROPERTIES[mutant])),
        "shared_properties": ";".join(shared),
        "lost_properties": ";".join(lost),
        "gained_properties": ";".join(gained),
    }
