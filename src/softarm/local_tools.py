from __future__ import annotations

import tomllib
from pathlib import Path
from typing import Any


def local_tool_config_path() -> Path:
    return Path(__file__).resolve().parents[2] / ".softarm.local.toml"


def load_local_tools(path: Path) -> dict[str, Any]:
    with path.open("rb") as stream:
        tools = tomllib.load(stream).get("tools")
    if not isinstance(tools, dict):
        raise ValueError(f"{path} has no [tools] table")
    return tools
