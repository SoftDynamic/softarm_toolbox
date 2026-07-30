from __future__ import annotations

import json
import os
from pathlib import Path
import shutil
import subprocess
import tempfile
import tomllib

import sympy as sp

from ..ast import decode, encode


def _kernel_path(explicit: str | None) -> Path:
    candidate = explicit or os.environ.get("SOFTARM_WOLFRAM_KERNEL")
    local = Path.cwd() / ".softarm.local.toml"
    if not candidate and local.is_file():
        with local.open("rb") as stream:
            candidate = tomllib.load(stream).get("tools", {}).get("wolfram_math")
    if not candidate:
        candidate = shutil.which("math") or shutil.which("wolframscript")
    path = Path(candidate) if candidate else Path("__missing_wolfram_kernel__")
    if not path.is_file():
        raise RuntimeError(
            "Wolfram backend requested but no kernel was found; pass --wolfram-kernel or set SOFTARM_WOLFRAM_KERNEL"
        )
    return path


def optimize(expressions: list[sp.Expr], kernel: str | None = None) -> list[sp.Expr]:
    executable = _kernel_path(kernel)
    bridge = Path(__file__).with_name("wolfram_bridge.wls")
    payload = {"expressions": [encode(expr) for expr in expressions]}
    with tempfile.TemporaryDirectory(prefix="softarm-wolfram-") as directory:
        root = Path(directory)
        input_path = root / "input.json"
        output_path = root / "output.json"
        input_path.write_text(json.dumps(payload), encoding="utf-8")
        completed = subprocess.run(
            [str(executable), "-script", str(bridge), str(input_path), str(output_path)],
            capture_output=True,
            text=True,
            timeout=900,
            check=False,
        )
        if completed.returncode != 0 or not output_path.is_file():
            message = (completed.stderr or completed.stdout).strip()
            raise RuntimeError(f"Wolfram backend failed: {message or 'no output produced'}")
        raw = output_path.read_text(encoding="utf-8")
        if not raw.strip():
            message = (completed.stderr or completed.stdout).strip()
            raise RuntimeError(f"Wolfram backend produced an empty AST: {message or 'no diagnostic'}")
        result = json.loads(raw)
    return [decode(node) for node in result["expressions"]]
