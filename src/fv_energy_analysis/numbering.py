from __future__ import annotations

from dataclasses import dataclass
from functools import lru_cache

from .models import Residue, ResidueKey
from .structure import residue_sequence


@dataclass(frozen=True)
class AntibodyAnnotation:
    chain_id: str
    chain_type: str
    sequence: str
    region_by_key: dict[ResidueKey, str]
    number_by_key: dict[ResidueKey, str]
    backend: str


def _load_chain_class():
    errors: list[str] = []
    try:
        from abnumber.future import Chain  # type: ignore
        return Chain, "abnumber.future/ANARCII"
    except Exception as exc:  # pragma: no cover - depends on optional package
        errors.append(f"abnumber.future: {exc}")
    try:
        from abnumber import Chain  # type: ignore
        return Chain, "abnumber/ANARCII"
    except Exception as exc:  # pragma: no cover - depends on optional package
        errors.append(f"abnumber: {exc}")
    raise RuntimeError(
        "AbNumber with ANARCII is required for CDR annotation. Install the antibody extra "
        "(`pip install -e .[antibody]`) in Linux/WSL. Import errors: " + "; ".join(errors)
    )


@lru_cache(maxsize=4096)
def _number_sequence(
    sequence: str, scheme: str, cdr_definition: str
) -> tuple[tuple[tuple[str, str, str], ...], str, str, int]:
    """Number a unique sequence once per persistent batch worker."""
    Chain, backend = _load_chain_class()
    attempts = [
        {"scheme": scheme, "cdr_definition": cdr_definition, "use_anarcii": True},
        {"scheme": scheme, "cdr_definition": cdr_definition},
    ]
    last_error: Exception | None = None
    numbered = None
    for kwargs in attempts:
        try:
            numbered = Chain(sequence, **kwargs)
            break
        except TypeError as exc:
            last_error = exc
    if numbered is None:
        raise RuntimeError(f"Unable to initialize AbNumber/ANARCII: {last_error}")

    entries: list[tuple[str, str, str]] = []
    for position, aa in numbered:
        region = position.get_region()
        entries.append((str(position), str(aa), str(region).upper()))
    variable_sequence = "".join(aa for _, aa, _ in entries)
    start = sequence.find(variable_sequence)
    if start < 0:
        raise ValueError(
            "ANARCII-numbered variable-domain sequence could not be "
            "mapped exactly back to the coordinate sequence. Check missing residues/alternate models."
        )
    return tuple(entries), str(numbered.chain_type), backend, start


def annotate_chain(
    residues: list[Residue],
    chain_id: str,
    scheme: str = "imgt",
    cdr_definition: str = "imgt",
) -> AntibodyAnnotation:
    sequence = residue_sequence(residues)
    entries_tuple, chain_type, backend, start = _number_sequence(
        sequence, scheme, cdr_definition
    )
    entries = list(entries_tuple)
    selected = residues[start:start + len(entries)]
    region_by_key = {res.key: region for res, (_, _, region) in zip(selected, entries)}
    number_by_key = {res.key: number for res, (number, _, _) in zip(selected, entries)}
    return AntibodyAnnotation(
        chain_id=chain_id,
        chain_type=chain_type,
        sequence=sequence,
        region_by_key=region_by_key,
        number_by_key=number_by_key,
        backend=backend,
    )


def annotate_fv(structure, chain_ids: list[str], scheme: str, cdr_definition: str) -> dict[str, AntibodyAnnotation]:
    output = {}
    for chain_id in chain_ids:
        residues = structure.chain_residues(chain_id)
        if not residues:
            raise ValueError(f"Fv chain {chain_id!r} is absent from {structure.path}")
        output[chain_id] = annotate_chain(residues, chain_id, scheme, cdr_definition)
    return output
