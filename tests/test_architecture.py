from __future__ import annotations

import ast
from pathlib import Path

ROOT = Path(__file__).parents[1]
DOMAIN_MODULES = (
    "actuation.py",
    "constraints.py",
    "derive.py",
    "dynamics.py",
    "geometry.py",
    "integration.py",
    "special.py",
)


def test_model_domain_does_not_import_infrastructure():
    forbidden = {"softarm.backends", "softarm.codegen", "softarm.cli"}
    for filename in DOMAIN_MODULES:
        path = ROOT / "src" / "softarm" / filename
        tree = ast.parse(path.read_text(encoding="utf-8"), filename)
        imports = []
        for node in ast.walk(tree):
            if isinstance(node, ast.Import):
                imports.extend(alias.name for alias in node.names)
            elif isinstance(node, ast.ImportFrom):
                module = node.module or ""
                if node.level:
                    module = "softarm." + module
                imports.append(module)
        assert not any(
            name == blocked or name.startswith(blocked + ".")
            for name in imports
            for blocked in forbidden
        ), f"{filename} imports infrastructure: {imports}"


def test_wolfram_bridge_contains_no_model_specific_vocabulary():
    bridge = (
        ROOT / "src" / "softarm" / "backends" / "wolfram_bridge.wls"
    ).read_text(encoding="utf-8").lower()
    for forbidden in (
        "euler", "cosserat", "tendon", "plane_point_contact"
    ):
        assert forbidden not in bridge


def test_legacy_backend_modules_are_removed():
    backend = ROOT / "src" / "softarm" / "backends"
    assert not (backend / "session.py").exists()
    assert not (backend / "sympy_backend.py").exists()
    assert not (backend / "wolfram_backend.py").exists()
