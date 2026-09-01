from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path
from typing import Any


def _import_casadi():
    try:
        import casadi as ca
    except ImportError as error:
        raise RuntimeError(
            "CasADi support requires the optional dependency; "
            "install softarm-toolbox[casadi]"
        ) from error
    return ca


@dataclass(frozen=True)
class CasadiBundle:
    """Lazy loader for a generated collection of CasADi Functions."""

    path: Path
    manifest: dict[str, Any]

    def function(self, name: str):
        try:
            relative = self.manifest["functions"][name]["file"]
        except KeyError as error:
            raise KeyError(f"unknown CasADi bundle function: {name}") from error
        return _import_casadi().Function.load(str(self.path / relative))

    @property
    def parameters(self) -> tuple[float, ...]:
        return tuple(float(item["default"]) for item in self.manifest["parameters"])


def load_casadi_bundle(path: str | Path) -> CasadiBundle:
    target = Path(path).resolve()
    manifest_path = target / "manifest.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    target_info = manifest.get("target")
    if not isinstance(target_info, dict) or target_info.get("name") != "casadi":
        raise ValueError(f"not a CasADi bundle: {target}")
    if target_info.get("schema_version") != 1:
        raise ValueError(
            f"unsupported CasADi bundle schema: {target_info.get('schema_version')!r}"
        )
    ca = _import_casadi()
    generated = str(target_info.get("casadi_version", ""))
    if generated.split(".")[:2] != ca.__version__.split(".")[:2]:
        raise RuntimeError(
            f"bundle uses CasADi {generated}, but runtime is {ca.__version__}"
        )
    return CasadiBundle(target, manifest)
