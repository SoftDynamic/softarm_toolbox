from __future__ import annotations

import json
import re
import shutil
import subprocess
import tomllib
from pathlib import Path

import pytest

from softarm.actuation import derive_actuation
from softarm.cli import main as cli_main
from softarm.codegen import generate_latex_document, generate_matlab_bundle
from softarm.codegen import latex as latex_module
from softarm.codegen.optimization import FunctionOptimizer
from softarm.config import (
    ActuationConfig,
    BaseConfig,
    ConstraintConfig,
    IntegrationConfig,
    ModelConfig,
    TendonChannelConfig,
    TendonSpanConfig,
)
from softarm.constraints import derive_constraint
from softarm.derive import derive

ROOT = Path(__file__).parents[1]


def _find_pdflatex() -> str | None:
    executable = shutil.which("pdflatex")
    local_config = ROOT / ".softarm.local.toml"
    if executable is not None or not local_config.is_file():
        return executable
    with local_config.open("rb") as stream:
        configured = tomllib.load(stream).get("tools", {}).get("pdflatex")
    if isinstance(configured, str) and Path(configured).is_file():
        return configured
    return None


PDFLATEX = _find_pdflatex()
REFERENCE_TEX = tuple(sorted(
    (ROOT / "examples" / "generated").glob("*/softarm_model.tex")
))


def _euler_config(**changes):
    values = dict(
        rod="euler_bernoulli",
        parameterization="ritz",
        segments=1,
        integration=IntegrationConfig("analytic"),
        ritz_x=(0.0, 0.0, 1.5, -0.5),
        ritz_y=(0.0, 0.0, 1.5, -0.5),
    )
    values.update(changes)
    return ModelConfig(**values)


def test_document_sections_are_model_specific_and_paper_friendly(tmp_path):
    fixed = derive(ModelConfig(
        rod="extensible_euler_bernoulli", parameterization="pcs", segments=1,
        inertia="lumped", integration=IntegrationConfig()
    ))
    first = tmp_path / "first"
    second = tmp_path / "second"
    generate_latex_document(fixed, first)
    generate_latex_document(fixed, second)
    text = (first / "softarm_model.tex").read_text(encoding="utf-8")
    assert text == (second / "softarm_model.tex").read_text(encoding="utf-8")
    assert "Piecewise-Constant-Strain Parameterization" in text
    assert "Ritz Parameterization" not in text
    assert "The arm base is fixed" in text
    assert r"\Sfun" in text and r"\Cfun" in text
    assert "SincSqrt" not in text and "CoscSqrt" not in text
    assert r"b_{x,1}" in text and r"L_{0,1}" in text
    assert "Exact Symbolic Appendix" not in text

    floating = derive(_euler_config(base=BaseConfig("floating_rpy")))
    generate_latex_document(floating, tmp_path / "floating")
    floating_text = (tmp_path / "floating/softarm_model.tex").read_text(encoding="utf-8")
    assert "Ritz Parameterization" in floating_text
    assert "Piecewise-Constant-Strain Parameterization" not in floating_text
    assert r"H_{WB}" in floating_text
    assert "ZYX convention" in floating_text
    assert r"M_B=" in floating_text


def test_actuation_constraint_and_parameter_order_match_bundle(tmp_path):
    channels = (
        TendonChannelConfig("left_1", "signed", (TendonSpanConfig(1, 0.02, 0.0),)),
        TendonChannelConfig("right_1", "signed", (TendonSpanConfig(1, 0.02, 3.14),)),
    )
    config = _euler_config(
        actuation=ActuationConfig("tendon", "strict", channels),
        constraint=ConstraintConfig("plane_point_contact", {
            "family": "plane_point_contact", "plane_normal": [0.0, 0.0, -1.0],
        }),
    )
    plant = derive(config)
    actuation = derive_actuation(plant)
    constraint = derive_constraint(plant)
    generate_matlab_bundle(
        plant, tmp_path, actuation=actuation, constraint=constraint
    )
    text = (tmp_path / "softarm_model.tex").read_text(encoding="utf-8")
    assert r"\section{Actuation}" in text
    assert "Strict actuator-coordinate acceleration" in text
    assert r"\section{Active Acceleration Constraint}" in text
    assert "plane-point implementation" in text
    manifest = json.loads((tmp_path / "manifest.json").read_text(encoding="utf-8"))
    names = [item["name"] for item in manifest["parameters"]]
    positions = [text.index(rf"\texttt{{{name.replace('_', r'\_')}}}") for name in names]
    assert positions == sorted(positions)


def test_default_document_does_not_invoke_appendix_cse(tmp_path, monkeypatch):
    plant = derive(_euler_config())

    def fail(*args, **kwargs):
        raise AssertionError("appendix CSE must remain lazy")

    monkeypatch.setattr(latex_module, "_cse_data", fail)
    generate_latex_document(plant, tmp_path)


def test_exact_appendix_is_present_and_cse_reconstructs_outputs(tmp_path):
    plant = derive(_euler_config())
    generate_latex_document(plant, tmp_path, include_appendix=True)
    text = (tmp_path / "softarm_model.tex").read_text(encoding="utf-8")
    assert "Exact Symbolic Appendix" in text
    assert "Common subexpressions" in text
    expressions = [
        plant.mass[0, 0], plant.mass[0, 1], plant.mass[1, 1],
        *list(plant.bias), *list(plant.end_jacobian),
    ]
    replacements, reduced = latex_module._cse_data(
        expressions, (len(expressions),), FunctionOptimizer()
    )
    assert latex_module._reconstruct_cse(replacements, reduced) == expressions


def test_cli_tex_appendix_flag_builds_same_document(tmp_path):
    config = tmp_path / "small_euler.toml"
    config.write_text(
        """
[model]
rod = "euler_bernoulli"
parameterization = "ritz"
segments = 1

[integration]
method = "analytic"

[ritz]
x = [0.0, 0.0, 1.5, -0.5]
y = [0.0, 0.0, 1.5, -0.5]
""".strip() + "\n",
        encoding="utf-8",
    )
    bundle = tmp_path / "bundle"
    assert cli_main([
        "build", str(config), "--target", "matlab", "--out", str(bundle),
        "--tex-appendix",
    ]) == 0
    assert "Exact Symbolic Appendix" in (
        bundle / "softarm_model.tex"
    ).read_text(encoding="utf-8")


@pytest.mark.skipif(PDFLATEX is None, reason="pdflatex is not configured")
def test_generated_document_compiles_with_pdflatex(tmp_path):
    plant = derive(_euler_config())
    path = generate_latex_document(plant, tmp_path)
    result = subprocess.run(
        [PDFLATEX, "-interaction=nonstopmode", "-halt-on-error", path.name],
        cwd=tmp_path,
        capture_output=True,
        text=True,
        timeout=60,
        check=False,
    )
    assert result.returncode == 0, result.stdout + result.stderr


@pytest.mark.skipif(PDFLATEX is None, reason="pdflatex is not configured")
@pytest.mark.parametrize("source", REFERENCE_TEX, ids=lambda path: path.parent.name)
def test_reference_documents_compile_with_pdflatex(source, tmp_path):
    result = subprocess.run(
        [
            PDFLATEX,
            "-interaction=nonstopmode",
            "-halt-on-error",
            "-file-line-error",
            str(source),
        ],
        cwd=tmp_path,
        capture_output=True,
        text=True,
        timeout=60,
        check=False,
    )
    assert result.returncode == 0, result.stdout + result.stderr


def test_cosserat_document_describes_pcs_and_lumped_inertia(tmp_path):
    plant = derive(ModelConfig(
        rod="cosserat", parameterization="pcs", segments=1, inertia="lumped",
        integration=IntegrationConfig(),
    ))
    generate_latex_document(plant, tmp_path)
    text = (tmp_path / "softarm_model.tex").read_text(encoding="utf-8")
    assert "Piecewise-Constant-Strain Parameterization" in text
    assert r"\kappa_{0,x,1}" in text
    assert r"\nu_{0,z,1}" in text
    assert r"\kappa_0_{x,1}" not in text
    assert re.search(r"\\qquad[A-Za-z]", text) is None
    assert r"\Tfun" in text
    assert "lumped inertia option" in text
    assert r"GA_{x,i}" in text and r"GJ_i" in text
