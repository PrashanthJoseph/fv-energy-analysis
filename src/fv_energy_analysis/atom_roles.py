"""Validated ligand atom-role cache readable without chemistry dependencies."""
from hashlib import sha256
import json
import os
from pathlib import Path


def _identity(params, template):
    params = Path(params)
    template = Path(template) if template else params.with_name(params.stem + "_0001.pdb")
    return params, {"params_sha256": sha256(params.read_bytes()).hexdigest(),
                    "template_sha256": sha256(template.read_bytes()).hexdigest()}


def save_atom_roles(params, template, chemistry):
    if not chemistry.get("valid") or not chemistry.get("estrogen_atom_roles"):
        return
    params, identity = _identity(params, template)
    path = params.with_suffix(".atom_roles.json")
    payload = {"schema": 1, **identity, "validated": True, "roles": chemistry["estrogen_atom_roles"]}
    if path.is_file() and json.loads(path.read_text()) == payload:
        return
    temporary = path.with_suffix(f".json.{os.getpid()}.tmp")
    temporary.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    temporary.replace(path)


def load_atom_roles(params, template):
    params, identity = _identity(params, template)
    path = params.with_suffix(".atom_roles.json")
    if not path.is_file():
        raise ValueError(f"Validated ligand atom-role cache missing: {path}; generate it in the preparation environment")
    payload = json.loads(path.read_text(encoding="utf-8"))
    if payload.get("schema") != 1 or payload.get("validated") is not True or any(payload.get(k) != v for k, v in identity.items()):
        raise ValueError("Ligand atom-role cache does not match current params/template; revalidate preparation")
    return payload["roles"]
