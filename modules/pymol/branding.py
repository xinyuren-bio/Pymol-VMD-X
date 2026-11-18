"""Display name for the desktop application and its command prompt."""

import json
from pathlib import Path

APP_NAME = 'Pymol-VMD-X'
COMMAND_PROMPT = APP_NAME + '>'


def local_build_info():
    """Metadata written by the desktop update command."""
    try:
        return json.loads(Path(__file__).with_name("_local_build.json").read_text(encoding="utf-8"))
    except (OSError, ValueError):
        return {}
