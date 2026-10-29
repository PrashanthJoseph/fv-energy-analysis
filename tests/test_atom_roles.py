import pytest
from fv_energy_analysis.atom_roles import load_atom_roles, save_atom_roles


def test_role_cache_requires_current_params_and_template(tmp_path):
    params, template = tmp_path / "E2L.params", tmp_path / "E2L_0001.pdb"
    params.write_text("validated parameters")
    template.write_text("validated template")
    roles = {"o3_atom": "O1", "o17_atom": "O2", "c18_methyl_atom": "C10", "a_ring_atoms": ["C1", "C2", "C3", "C4", "C5", "C6"]}
    save_atom_roles(params, template, {"valid": True, "estrogen_atom_roles": roles})
    assert load_atom_roles(params, template) == roles
    template.write_text("changed template")
    with pytest.raises(ValueError, match="does not match"):
        load_atom_roles(params, template)


def test_invalid_chemistry_cannot_create_a_role_cache(tmp_path):
    params = tmp_path / "E2L.params"
    save_atom_roles(params, None, {"valid": False, "estrogen_atom_roles": {"o3_atom": "O1"}})
    assert not params.with_suffix('.atom_roles.json').exists()
