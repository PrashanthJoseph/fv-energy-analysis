import json
from pathlib import Path
import sqlite3

import pytest

from fv_energy_analysis.flex_ddg import parse_mutations, validate_input, read_scores, gam_score


def test_mutation_dictionary_combined_and_insertion_codes():
    mutations = parse_mutations({'H':['E99F','R100D'], 'L':'Y91AF'})
    assert len(mutations) == 3
    assert mutations[-1]['number'] == 91
    assert mutations[-1]['insertion'] == 'A'
    assert mutations[-1]['new'] == 'F'


@pytest.mark.parametrize('value', [{}, {'H':'E99E'}, {'H':['E99F','E99W']},
                                  {'Heavy':'E99F'}, {'H':'X99F'}, {'H':99}])
def test_invalid_mutations(value):
    with pytest.raises(ValueError):
        parse_mutations(value)


def test_original_amino_acid_checked_against_pdb(tmp_path):
    pytest.importorskip("rdkit", reason="Flex ddG input chemical validation requires the analysis extra")
    pdb, params = tmp_path/'in.pdb', tmp_path/'lig.params'
    pdb.write_text('ATOM      1  CA  GLU H  99       0.000   0.000   0.000  1.00  0.00           C  \n'
                   'HETATM    2  C1  E2L X   1       1.000   0.000   0.000  1.00  0.00           C  \n')
    params.write_text('ATOM C1 CH1 X 0.00\n')
    (tmp_path/'lig_0001.pdb').write_text(pdb.read_text().splitlines()[1] + '\n')
    assert validate_input(pdb, params, 'X', parse_mutations({'H':'E99F'}))['valid']
    with pytest.raises(ValueError, match='expected R, found E'):
        validate_input(pdb, params, 'X', parse_mutations({'H':'R99F'}))


def test_four_state_score_sign_and_missing_state(tmp_path):
    database = tmp_path/'ddG.db3'
    with sqlite3.connect(database) as db:
        db.executescript('''CREATE TABLE batches(batch_id INTEGER,name TEXT);
            CREATE TABLE score_types(batch_id INTEGER,score_type_id INTEGER,score_type_name TEXT);
            CREATE TABLE structure_scores(batch_id INTEGER,score_type_id INTEGER,score_value REAL);''')
        # WT binding=-10, mutant binding=-13 -> ddG=-3 (improved).
        for i, (state, value) in enumerate(zip(
                ['bound_wt','unbound_wt','bound_mut','unbound_mut'], [-110,-100,-115,-102]),1):
            db.execute('INSERT INTO batches VALUES (?,?)', (i,state+'_dbreport'))
            db.execute('INSERT INTO score_types VALUES (?,1,?)', (i,'total_score'))
            db.execute('INSERT INTO structure_scores VALUES (?,1,?)', (i,value))
    result = read_scores(database)
    assert result['ddg'] == -3
    assert result['wt_dg'] == -10
    assert result['mut_dg'] == -13
    with sqlite3.connect(database) as db:
        db.execute('DELETE FROM structure_scores WHERE batch_id=4')
    with pytest.raises(ValueError, match='Missing four-state scores'):
        read_scores(database)


def test_gam_sign_and_zero():
    assert gam_score({}) == 0
    assert gam_score({'fa_atr':1}) > 0
    assert gam_score({'fa_atr':-1}) == -gam_score({'fa_atr':1})
