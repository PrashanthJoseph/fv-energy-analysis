"""Mutation-centred interaction-network change classification."""

from __future__ import annotations

from collections import Counter


ENERGY_TERMS = (
    "total", "fa_atr", "fa_rep", "fa_sol", "fa_elec", "lk_ball_wtd",
    "hbond_total", "hbond_sc", "hbond_bb_sc", "hbond_sr_bb", "hbond_lr_bb",
    "rama_prepro", "dslf_fa13", "nonbonded_total",
)

# Nominal side-chain formal charges near neutral pH. Histidine is reported as
# neutral/titratable (0) rather than assuming a protonation state not encoded by
# the residue name. E2L was parameterized with net charge zero.
AA_CHARGE = {"D": -1.0, "E": -1.0, "K": 1.0, "R": 1.0}
RESNAME_CHARGE = {
    "ASP": -1.0, "GLU": -1.0, "LYS": 1.0, "ARG": 1.0,
    "ALA": 0.0, "ASN": 0.0, "CYS": 0.0, "GLN": 0.0, "GLY": 0.0,
    "HIS": 0.0, "ILE": 0.0, "LEU": 0.0, "MET": 0.0, "PHE": 0.0,
    "PRO": 0.0, "SER": 0.0, "THR": 0.0, "TRP": 0.0, "TYR": 0.0,
    "VAL": 0.0, "E2L": 0.0, "LIG": 0.0,
}


def _float(row: dict[str, object], key: str) -> float:
    value = row.get(key, "")
    if value in {"", None}:
        return 0.0
    try:
        return float(value)
    except (TypeError, ValueError):
        return 0.0


def _types(row: dict[str, object], prefix: str) -> set[str]:
    values = {
        value for value in str(row.get(f"{prefix}_interaction_types", "")).split(";")
        if value and value != "0"
    }
    if _float(row, f"{prefix}_count_rosetta_hbond") > 0:
        values.add("rosetta_hbond")
    if abs(_float(row, f"{prefix}_rama_prepro")) > 1e-12:
        values.add('polymer_conformational')
    if abs(_float(row, f"{prefix}_dslf_fa13")) > 1e-12:
        values.add('disulfide')
    return values


def _type_counts(row: dict[str, object], prefix: str, names: set[str]) -> dict[str, int]:
    counts = {}
    for name in names:
        key = f"{prefix}_count_rosetta_hbond" if name == "rosetta_hbond" else f"{prefix}_count_{name}"
        counts[name] = int(round(_float(row, key)))
    return counts


def _other_residue(row: dict[str, object], prefix: str, mutation_label: str) -> tuple[str, str]:
    residue1 = str(row.get(f"{prefix}_residue1", ""))
    residue2 = str(row.get(f"{prefix}_residue2", ""))
    if residue1 == mutation_label:
        return residue2, str(row.get(f"{prefix}_residue2_name", ""))
    if residue2 == mutation_label:
        return residue1, str(row.get(f"{prefix}_residue1_name", ""))
    return "", ""


def mutation_network_changes(
    comparison: list[dict[str, object]],
    substitutions: list[dict[str, object]],
    energy_presence_cutoff: float = 0.1,
    energy_change_cutoff: float = 0.25,
) -> tuple[list[dict[str, object]], list[dict[str, object]], list[dict[str, object]]]:
    """Return detailed and per-mutation summarized interaction-network changes.

    A pair is considered present if it has a typed structural contact, a detected
    Rosetta H-bond, or an absolute pair energy at least ``energy_presence_cutoff``.
    This retains long-range energetic interactions that contact classifiers can
    miss while excluding numerical noise.
    """
    details: list[dict[str, object]] = []
    summaries: list[dict[str, object]] = []
    type_energy_rows: list[dict[str, object]] = []
    for mutation in substitutions:
        parent_label = str(mutation["parent_residue"])
        mutant_label = str(mutation["mutant_residue"])
        mutation_rows = [
            row for row in comparison
            if parent_label in {str(row.get("parent_residue1", "")), str(row.get("parent_residue2", ""))}
            or mutant_label in {str(row.get("mutant_residue1", "")), str(row.get("mutant_residue2", ""))}
        ]
        mutation_details: list[dict[str, object]] = []
        parent_mutation_charge = AA_CHARGE.get(str(mutation.get("parent_aa", "")), 0.0)
        mutant_mutation_charge = AA_CHARGE.get(str(mutation.get("mutant_aa", "")), 0.0)
        for row in mutation_rows:
            p_types, m_types = _types(row, "parent"), _types(row, "mutant")
            p_total, m_total = _float(row, "parent_total"), _float(row, "mutant_total")
            p_present = bool(p_types) or abs(p_total) >= energy_presence_cutoff
            m_present = bool(m_types) or abs(m_total) >= energy_presence_cutoff
            if not p_present and not m_present:
                continue
            if p_present and not p_types:
                p_types = {"energy_only"}
            if m_present and not m_types:
                m_types = {"energy_only"}
            all_types = p_types | m_types
            p_counts, m_counts = _type_counts(row, "parent", all_types), _type_counts(row, "mutant", all_types)
            formed_types = sorted(m_types - p_types)
            broken_types = sorted(p_types - m_types)
            retained_types = sorted(p_types & m_types)
            count_changes = [
                f"{name}:{p_counts[name]}->{m_counts[name]}"
                for name in sorted(all_types) if p_counts[name] != m_counts[name]
            ]
            energy_delta = m_total - p_total
            if not p_present and m_present:
                status = "formed"
            elif p_present and not m_present:
                status = "broken"
            elif formed_types or broken_types or count_changes or abs(energy_delta) >= energy_change_cutoff:
                status = "retained_changed"
            else:
                status = "retained_stable"
            if status == "formed":
                event_basis = (
                    "energy_threshold_only" if m_types == {"energy_only"}
                    else "typed_contact_or_hbond"
                )
            elif status == "broken":
                event_basis = (
                    "energy_threshold_only" if p_types == {"energy_only"}
                    else "typed_contact_or_hbond"
                )
            else:
                event_basis = "retained_pair"
            p_partner, p_partner_name = _other_residue(row, "parent", parent_label)
            m_partner, m_partner_name = _other_residue(row, "mutant", mutant_label)
            parent_partner_charge = RESNAME_CHARGE.get(p_partner_name, 0.0)
            mutant_partner_charge = RESNAME_CHARGE.get(m_partner_name, 0.0)
            detail: dict[str, object] = {
                "mutation": mutation["substitution"],
                "chain": mutation["chain"],
                "antibody_number": mutation["parent_antibody_number"],
                "region": mutation["parent_region"],
                "parent_mutated_residue": parent_label,
                "mutant_mutated_residue": mutant_label,
                "parent_partner": p_partner,
                "mutant_partner": m_partner,
                "parent_partner_name": p_partner_name,
                "mutant_partner_name": m_partner_name,
                "parent_mutation_nominal_charge": parent_mutation_charge,
                "mutant_mutation_nominal_charge": mutant_mutation_charge,
                "parent_partner_nominal_charge": parent_partner_charge,
                "mutant_partner_nominal_charge": mutant_partner_charge,
                "parent_pair_nominal_charge": parent_mutation_charge + parent_partner_charge,
                "mutant_pair_nominal_charge": mutant_mutation_charge + mutant_partner_charge,
                "delta_pair_nominal_charge": (
                    mutant_mutation_charge + mutant_partner_charge
                    - parent_mutation_charge - parent_partner_charge
                ),
                "pair_scope": row.get("mutant_pair_scope") or row.get("parent_pair_scope", ""),
                "network_status": status,
                "event_basis": event_basis,
                "parent_interaction_types": ";".join(sorted(p_types)) or "none",
                "mutant_interaction_types": ";".join(sorted(m_types)) or "none",
                "broken_interaction_types": ";".join(broken_types) or "none",
                "formed_interaction_types": ";".join(formed_types) or "none",
                "retained_interaction_types": ";".join(retained_types) or "none",
                "type_count_changes": ";".join(count_changes) or "none",
                "parent_total_contact_count": _float(row, "parent_total_contact_count"),
                "mutant_total_contact_count": _float(row, "mutant_total_contact_count"),
                "delta_total_contact_count": _float(row, "mutant_total_contact_count") - _float(row, "parent_total_contact_count"),
                "parent_pair_energy_strength": row.get("parent_pair_energy_strength", "absent") or "absent",
                "mutant_pair_energy_strength": row.get("mutant_pair_energy_strength", "absent") or "absent",
                "broken_type_pair_energies_REU": ";".join(
                    f"{name}:{p_total:.6f}" for name in broken_types
                ) or "none",
                "formed_type_pair_energies_REU": ";".join(
                    f"{name}:{m_total:.6f}" for name in formed_types
                ) or "none",
            }
            for term in ENERGY_TERMS:
                parent_value = _float(row, f"parent_{term}")
                mutant_value = _float(row, f"mutant_{term}")
                detail[f"parent_{term}"] = parent_value
                detail[f"mutant_{term}"] = mutant_value
                detail[f"delta_{term}"] = mutant_value - parent_value
            mutation_details.append(detail)
            details.append(detail)

        status_counts = Counter(str(row["network_status"]) for row in mutation_details)
        broken_type_counts: Counter[str] = Counter()
        formed_type_counts: Counter[str] = Counter()
        for row in mutation_details:
            broken_type_counts.update(
                value for value in str(row["broken_interaction_types"]).split(";") if value != "none"
            )
            formed_type_counts.update(
                value for value in str(row["formed_interaction_types"]).split(";") if value != "none"
            )
        summaries.append({
            "mutation": mutation["substitution"],
            "antibody_number": mutation["parent_antibody_number"],
            "region": mutation["parent_region"],
            "partner_count": len(mutation_details),
            "formed_pairs": status_counts["formed"],
            "broken_pairs": status_counts["broken"],
            "retained_changed_pairs": status_counts["retained_changed"],
            "retained_stable_pairs": status_counts["retained_stable"],
            "formed_interaction_types": ";".join(
                f"{name}:{count}" for name, count in sorted(formed_type_counts.items())
            ) or "none",
            "broken_interaction_types": ";".join(
                f"{name}:{count}" for name, count in sorted(broken_type_counts.items())
            ) or "none",
            "sum_pair_energy_delta_REU": sum(_float(row, "delta_total") for row in mutation_details),
            "sum_nonbonded_pair_energy_delta_REU": sum(_float(row, "delta_nonbonded_total") for row in mutation_details),
            "parent_mutation_nominal_charge": parent_mutation_charge,
            "mutant_mutation_nominal_charge": mutant_mutation_charge,
            "mutation_charge_change": mutant_mutation_charge - parent_mutation_charge,
        })

        for event, field, energy_field in (
            ("formed", "formed_interaction_types", "mutant_total"),
            ("broken", "broken_interaction_types", "parent_total"),
        ):
            interaction_types = sorted({
                interaction_type
                for row in mutation_details
                for interaction_type in str(row[field]).split(";")
                if interaction_type != "none"
            })
            for interaction_type in interaction_types:
                matching = [
                    row for row in mutation_details
                    if interaction_type in str(row[field]).split(";")
                ]
                aggregate: dict[str, object] = {
                    "mutation": mutation["substitution"],
                    "antibody_number": mutation["parent_antibody_number"],
                    "region": mutation["parent_region"],
                    "event": event,
                    "interaction_type": interaction_type,
                    "pair_count": len(matching),
                    "contributing_pairs": ";".join(sorted({
                        f"{row.get('mutant_mutated_residue') or row.get('parent_mutated_residue')}"
                        f"<->{row.get('mutant_partner') or row.get('parent_partner')}"
                        for row in matching
                    })),
                    "sum_corresponding_pair_energy_REU": sum(
                        _float(row, energy_field) for row in matching
                    ),
                    "sum_parent_total_REU": sum(_float(row, "parent_total") for row in matching),
                    "sum_mutant_total_REU": sum(_float(row, "mutant_total") for row in matching),
                    "sum_delta_total_REU": sum(_float(row, "delta_total") for row in matching),
                    "parent_mutation_nominal_charge": parent_mutation_charge,
                    "mutant_mutation_nominal_charge": mutant_mutation_charge,
                    "mutation_charge_change": mutant_mutation_charge - parent_mutation_charge,
                    "overlap_warning": "Pair energies repeat across types when one pair has multiple interaction annotations",
                }
                for term in ENERGY_TERMS[1:]:
                    aggregate[f"sum_parent_{term}"] = sum(_float(row, f"parent_{term}") for row in matching)
                    aggregate[f"sum_mutant_{term}"] = sum(_float(row, f"mutant_{term}") for row in matching)
                    aggregate[f"sum_delta_{term}"] = sum(_float(row, f"delta_{term}") for row in matching)
                type_energy_rows.append(aggregate)
    return details, summaries, type_energy_rows
