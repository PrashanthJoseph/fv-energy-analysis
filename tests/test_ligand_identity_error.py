"""A rejected ligand must not break the batch's worker pool."""
from concurrent.futures import ProcessPoolExecutor
import multiprocessing
import pickle

import pytest

from fv_energy_analysis.ligand_chemistry import LigandIdentityError


DIAGNOSTIC = {"differences": ["stereochemistry"],
              "expected": {"inchikey": "expected"},
              "observed": {"inchikey": "observed"},
              "stereocentre_differences": [{"atom_name": "C18", "expected_CIP": "S", "observed_CIP": "R"}]}


def reject_ligand():
    raise LigandIdentityError(DIAGNOSTIC)


def subsequent_work():
    return "worker remains usable"


def test_ligand_error_preserves_structured_diagnostics_through_pickle():
    error = LigandIdentityError(DIAGNOSTIC)
    restored = pickle.loads(pickle.dumps(error))
    assert restored.diagnostics == DIAGNOSTIC
    assert str(restored) == str(error)


def test_rejected_ligand_does_not_destroy_process_pool():
    with ProcessPoolExecutor(max_workers=1, mp_context=multiprocessing.get_context("spawn")) as pool:
        with pytest.raises(LigandIdentityError) as rejected:
            pool.submit(reject_ligand).result(timeout=30)
        assert rejected.value.diagnostics == DIAGNOSTIC
        assert pool.submit(subsequent_work).result(timeout=30) == "worker remains usable"
