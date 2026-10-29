from __future__ import annotations

import shlex
from hashlib import sha256
from pathlib import Path


_PYROSETTA = None
_INITIALIZED_PARAMS: tuple[str, ...] | None = None
_INITIALIZED_PARAM_DIGESTS: tuple[str, ...] | None = None
_INITIALIZED_LEGACY: bool | None = None
_SCORE_WEIGHTS_ENGINE = None
_SCORE_WEIGHTS_SIGNATURE: dict[str, object] | None = None
SCORING_PROTOCOL_VERSION = 2


def _legacy_weights(weights: str) -> bool:
    """Historical ligand weights require historical energy-method defaults."""
    return Path(weights).stem in {"ligand", "ligand_soft_rep"}


def resolve_legacy_pre_talaris(weights: str, override: bool | None = None) -> bool:
    """Resolve automatic historical defaults or an explicit custom-file choice."""
    if override is not None and not isinstance(override, bool):
        raise ValueError("legacy_pre_talaris must be true, false or null (automatic)")
    return _legacy_weights(weights) if override is None else override


def scoring_weights_identity(
    weights: str | Path,
    legacy_pre_talaris: bool | None = None,
    *,
    database: str | Path | None = None,
) -> dict[str, object]:
    """Describe a Rosetta weight set without initializing the licensed engine.

    Existing custom files are resolved and hashed. A database may be supplied
    by a scientific worker to resolve and hash a built-in weight file too.
    Config readers can use this helper without PyRosetta installed.
    """
    value = str(weights).strip()
    if not value:
        raise ValueError("weights must name a Rosetta weight set or an existing .wts file")
    candidate = Path(value).expanduser()
    explicit_path = candidate.is_file()
    if not explicit_path and (candidate.is_absolute() or "/" in value or "\\" in value):
        raise FileNotFoundError(f"Custom Rosetta weights file was not found: {candidate}")
    weight_file = candidate.resolve() if explicit_path else None
    if weight_file is None and database is not None:
        folder = Path(database) / "scoring" / "weights"
        choices = (folder / value, folder / (value + ".wts"))
        weight_file = next((path for path in choices if path.is_file()), None)
        if weight_file is None:
            raise FileNotFoundError(f"Rosetta weights {value!r} were not found in {folder}")
    return {
        "weights": str(candidate.resolve()) if explicit_path else value,
        "legacy_pre_talaris": resolve_legacy_pre_talaris(value, legacy_pre_talaris),
        "weights_sha256": sha256(weight_file.read_bytes()).hexdigest() if weight_file else None,
    }


def _initialize(params: list[str] | None = None, legacy: bool = False) -> object:
    """Initialize PyRosetta once per process and reuse it for later structures."""
    global _PYROSETTA, _INITIALIZED_PARAMS, _INITIALIZED_PARAM_DIGESTS, _INITIALIZED_LEGACY
    resolved = tuple(str(Path(param).expanduser().resolve()) for param in (params or []))
    digests = tuple(sha256(Path(param).read_bytes()).hexdigest() for param in resolved)
    if _PYROSETTA is not None:
        if legacy != _INITIALIZED_LEGACY:
            raise RuntimeError("Legacy and modern scoring must use separate PyRosetta processes")
        if resolved != (_INITIALIZED_PARAMS or ()) or digests != (_INITIALIZED_PARAM_DIGESTS or ()):
            raise RuntimeError(
                "A persistent PyRosetta worker cannot change ligand params after initialization: "
                f"initialized={_INITIALIZED_PARAMS or ()}, requested={resolved}"
            )
        return _PYROSETTA
    try:
        import pyrosetta  # type: ignore
    except ImportError as exc:  # pragma: no cover - licensed optional dependency
        raise RuntimeError("PyRosetta is required for score-pyrosetta") from exc
    options = ["-mute", "all"]
    if legacy:
        options += ["-restore_pre_talaris_2013_behavior", "true"]
    for param in resolved:
        options += ["-extra_res_fa", param]
    pyrosetta.init(" ".join(shlex.quote(value) for value in options), silent=True)
    _PYROSETTA = pyrosetta
    _INITIALIZED_PARAMS = resolved
    _INITIALIZED_PARAM_DIGESTS = digests
    _INITIALIZED_LEGACY = legacy
    return _PYROSETTA


def score_and_extract_hbonds(
    pdb: str | Path,
    pair_output: str | Path,
    interface_output: str | Path,
    hbond_output: str | Path,
    ligand_chain: str,
    params: list[str] | None = None,
    weights: str = "ref2015",
    legacy_pre_talaris: bool | None = None,
) -> None:
    """Load one pose and produce pair, interface, and H-bond outputs in one session."""
    from .hbonds import write_pose_hbonds

    pyrosetta = _initialize(params, legacy=resolve_legacy_pre_talaris(weights, legacy_pre_talaris))
    pose = pyrosetta.pose_from_file(str(pdb))
    scorefxn = _scorefunction(pyrosetta, weights)
    description = Path(pdb).stem
    write_pair_energies(pose, scorefxn, pair_output, description)
    write_interface_metrics(pose, scorefxn, ligand_chain, interface_output, description)
    write_pose_hbonds(pose, hbond_output, scorefxn=scorefxn)


def _scorefunction(pyrosetta, weights: str):
    global _SCORE_WEIGHTS_ENGINE, _SCORE_WEIGHTS_SIGNATURE
    identity = scoring_weights_identity(weights)
    if _SCORE_WEIGHTS_ENGINE is pyrosetta and _SCORE_WEIGHTS_SIGNATURE != identity:
        raise RuntimeError(
            "A persistent PyRosetta worker cannot change scoring weights or custom-file contents; "
            "start a new worker process for the new score function"
        )
    scorefxn = pyrosetta.create_score_function(weights)
    options = scorefxn.energy_method_options()
    hbond_options = options.hbond_options()
    hbond_options.decompose_bb_hb_into_pair_energies(True)
    options.hbond_options(hbond_options)
    scorefxn.set_energy_method_options(options)
    _SCORE_WEIGHTS_ENGINE, _SCORE_WEIGHTS_SIGNATURE = pyrosetta, identity
    return scorefxn


def _active_terms(scorefxn) -> list[object]:
    return list(scorefxn.get_nonzero_weighted_scoretypes())


def _term_name(score_type: object) -> str:
    value = str(score_type)
    return value.rsplit(".", 1)[-1]


def _pdb_id(pdb_info, index: int) -> str:
    chain = str(pdb_info.chain(index)).strip() or "_"
    number = int(pdb_info.number(index))
    insertion = str(pdb_info.icode(index)).strip()
    return f"{number}{insertion}{chain}"


def complete_pair_energies_array(pose, scorefxn):
    """Include Rosetta's separately stored long-range two-body containers.

    The built-in residue_pair_energies_array reads only EnergyGraph; in ref2015
    it omits both rama_prepro and dslf_fa13. Match residue_energy_breakdown's
    upper-neighbor traversal and apply weights exactly once.
    """
    import numpy as np
    from pyrosetta.rosetta.core.scoring import EMapVector, methods
    scorefxn(pose)
    terms = _active_terms(scorefxn)
    names = [_term_name(term) for term in terms]
    short = pose.energies().residue_pair_energies_array()
    result = np.zeros(short.shape, dtype=[("total_score", float)] + [(name, float) for name in names])
    for name in names:
        if name in short.dtype.names:
            result[name] = short[name]
    for lr in range(1, int(methods.n_long_range_types) + 1):
        container = pose.energies().long_range_container(methods.LongRangeEnergyType(lr))
        if container is None or container.empty():
            continue
        for i in range(1, pose.size() + 1):
            iterator = container.const_upper_neighbor_iterator_begin(i)
            end = container.const_upper_neighbor_iterator_end(i)
            while iterator != end:
                j = int(iterator.upper_neighbor_id())
                energy = EMapVector()
                iterator.accumulate_energy(energy)
                for term, name in zip(terms, names):
                    result[name][i - 1, j - 1] += float(energy[term])
                iterator.pre_increment()
    for term, name in zip(terms, names):
        result['total_score'] += result[name] * float(scorefxn.weights()[term])
    return result


def write_pair_energies(pose, scorefxn, output: str | Path, description: str) -> None:
    """Write a residue_energy_breakdown-compatible weighted pair table."""
    scorefxn(pose)
    terms = _active_terms(scorefxn)
    weights = scorefxn.weights()
    names = [_term_name(term) for term in terms]
    pair_array = complete_pair_energies_array(pose, scorefxn)
    available = set(pair_array.dtype.names or ())
    pdb_info = pose.pdb_info()
    header = [
        "SCORE:", "pose_id", "resi1", "pdbid1", "restype1",
        "resi2", "pdbid2", "restype2", *names, "total", "description",
    ]
    lines = [" ".join(header)]
    # Self energies are separate records, not inter-residue interactions.
    for i in range(1, pose.total_residue() + 1):
        energy = pose.energies().onebody_energies(i)
        values = [float(energy[term]) * float(weights[term]) for term in terms]
        lines.append(" ".join([
            "SCORE:", description, str(i), _pdb_id(pdb_info, i), pose.residue(i).name3(),
            "--", "--", "onebody", *(f"{value:.6f}" for value in values),
            f"{sum(values):.6f}", f"{description}_{i}_onebody",
        ]))
    for i in range(1, pose.total_residue() + 1):
        for j in range(i + 1, pose.total_residue() + 1):
            pair_record = pair_array[i - 1, j - 1]
            values = [
                float(pair_record[name]) * float(weights[term]) if name in available else 0.0
                for term, name in zip(terms, names)
            ]
            total = float(pair_record["total_score"])
            if abs(total) < 1e-12 and not any(abs(value) >= 1e-12 for value in values):
                continue
            row = [
                "SCORE:", description, str(i), _pdb_id(pdb_info, i), pose.residue(i).name3(),
                str(j), _pdb_id(pdb_info, j), pose.residue(j).name3(),
                *(f"{value:.6f}" for value in values), f"{total:.6f}", f"{description}_{i}_{j}",
            ]
            lines.append(" ".join(row))
    Path(output).write_text("\n".join(lines) + "\n", encoding="utf-8")


def _weighted_totals(pose, scorefxn, terms: list[object]) -> dict[str, float]:
    scorefxn(pose)
    energy_map = pose.energies().total_energies()
    weights = scorefxn.weights()
    return {_term_name(term): float(energy_map[term]) * float(weights[term]) for term in terms}


def write_interface_metrics(
    pose,
    scorefxn,
    ligand_chain: str,
    output: str | Path,
    description: str,
) -> None:
    """Score bound minus rigid-body-separated states using the ligand's jump."""
    from pyrosetta.rosetta.protocols.rigid import RigidBodyTransMover  # type: ignore

    pdb_info = pose.pdb_info()
    ligand_residues = [
        index for index in range(1, pose.total_residue() + 1)
        if str(pdb_info.chain(index)).strip() == ligand_chain
    ]
    if not ligand_residues:
        raise ValueError(f"Ligand chain {ligand_chain!r} is absent from the PyRosetta pose")
    jump = int(pose.fold_tree().get_jump_that_builds_residue(ligand_residues[0]))
    if jump <= 0:
        raise ValueError(
            f"Could not identify a rigid-body jump that builds ligand chain {ligand_chain!r}; "
            "place the ligand after the protein chains in the input PDB"
        )
    terms = _active_terms(scorefxn)
    bound = _weighted_totals(pose, scorefxn, terms)
    bound_total = float(scorefxn(pose))
    separated = pose.clone()
    mover = RigidBodyTransMover(separated, jump)
    mover.step_size(1000.0)
    mover.apply(separated)
    apart = _weighted_totals(separated, scorefxn, terms)
    apart_total = float(scorefxn(separated))
    prefix = f"if_{ligand_chain}_"
    deltas = {prefix + name: bound[name] - apart[name] for name in bound}
    values = {
        "total_score": bound_total,
        "bound_total": bound_total,
        "separated_total": apart_total,
        f"interface_delta_{ligand_chain}": bound_total - apart_total,
        **deltas,
    }
    header = ["SCORE:", *values.keys(), "description"]
    row = ["SCORE:", *(f"{value:.6f}" for value in values.values()), description]
    Path(output).write_text(" ".join(header) + "\n" + " ".join(row) + "\n", encoding="utf-8")


def score_with_pyrosetta(
    pdb: str | Path,
    pair_output: str | Path,
    interface_output: str | Path | None = None,
    ligand_chain: str | None = None,
    params: list[str] | None = None,
    weights: str = "ref2015",
    legacy_pre_talaris: bool | None = None,
) -> None:
    pyrosetta = _initialize(params, legacy=resolve_legacy_pre_talaris(weights, legacy_pre_talaris))
    pose = pyrosetta.pose_from_file(str(pdb))
    scorefxn = _scorefunction(pyrosetta, weights)
    description = Path(pdb).stem
    write_pair_energies(pose, scorefxn, pair_output, description)
    if interface_output:
        if not ligand_chain:
            raise ValueError("--ligand-chain is required when --interface-output is supplied")
        write_interface_metrics(pose, scorefxn, ligand_chain, interface_output, description)
