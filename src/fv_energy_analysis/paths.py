"""Resolve saved Windows/WSL paths in the environment reading the report."""
import os
from pathlib import Path
import re


def native_path(value: str | Path) -> Path:
    text = str(value)
    mounted = re.match(r'^/mnt/([A-Za-z])/(.*)$', text)
    drive = re.match(r'^([A-Za-z]):[\\/](.*)$', text)
    if os.name == 'nt' and mounted:
        return Path(f'{mounted.group(1).upper()}:/{mounted.group(2)}')
    if os.name != 'nt' and drive:
        return Path(f'/mnt/{drive.group(1).lower()}/' + drive.group(2).replace('\\','/'))
    return Path(text)
