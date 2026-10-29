"""Persist a chosen new hydroxyl H without changing any other PDB field."""
import numpy as np


def replace_hydroxyl_hydrogen(data: bytes, xyz, name="H25") -> bytes:
    values = np.asarray(xyz, dtype=float)
    if values.shape != (3,) or not np.isfinite(values).all():
        raise ValueError("Expected three finite hydroxyl-H coordinates")
    lines = data.decode("utf-8").splitlines(keepends=True)
    selected = [i for i, line in enumerate(lines) if line.startswith("HETATM")
                and line[17:20] == "E3L" and line[21:22] == "X" and line[12:16].strip() == name]
    if len(selected) != 1:
        raise ValueError(f"Expected exactly one E3 hydroxyl proton {name}")
    index = selected[0]
    lines[index] = lines[index][:30] + "".join(f"{value:8.3f}" for value in values) + lines[index][54:]
    return "".join(lines).encode("utf-8")
