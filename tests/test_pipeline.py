from __future__ import annotations

from pathlib import Path

import pytest
import sympy as sp

from softarm import pipeline as pipeline_module
from softarm.cli import main as cli_main
from softarm.codegen.optimization import SympyCse
from softarm.config import IntegrationConfig, ModelConfig
from softarm.dynamics import SympyBatchDifferentiator
from softarm.pipeline import BuildError, BuildOptions, build_bundle


def _config() -> ModelConfig:
    return ModelConfig(
        rod="euler_bernoulli",
        parameterization="ritz",
        segments=1,
        integration=IntegrationConfig("analytic"),
        ritz_x=(0.0, 0.0, 1.5, -0.5),
        ritz_y=(0.0, 0.0, 1.5, -0.5),
    )


class RecordingKernel:
    events: list[str] = []
    fail_cse = False

    def __init__(self, kernel=None, timeout=600.0):
        self.events.append(f"start:{kernel}:{timeout}")

    def __enter__(self):
        self.events.append("enter")
        return self

    def __exit__(self, exc_type, exc, traceback):
        self.events.append("close")

    def differentiate(self, expressions, variables):
        self.events.append("differentiate")
        return SympyBatchDifferentiator().differentiate(expressions, variables)

    def factor_terms(self, expressions):
        self.events.append("factor_terms")
        return list(expressions)

    def cse(self, expressions, prefix="t", order="none"):
        self.events.append("wolfram_cse")
        if self.fail_cse:
            raise RuntimeError("experimental CSE failed")
        return SympyCse().cse(expressions, prefix, order)


def _fake_generate(
    plant,
    output,
    actuation=None,
    constraint=None,
    tex_appendix=False,
    optimizer=None,
):
    del plant, actuation, constraint, tex_appendix
    x = sp.Symbol("x")
    optimizer.optimize([x * (x + 1)], (1,))
    return Path(output)


def test_explicit_wolfram_steps_share_one_kernel(monkeypatch, tmp_path):
    RecordingKernel.events = []
    RecordingKernel.fail_cse = False
    monkeypatch.setattr(pipeline_module, "WolframKernel", RecordingKernel)
    monkeypatch.setattr(
        pipeline_module, "generate_matlab_bundle", _fake_generate
    )
    result = build_bundle(
        _config(),
        tmp_path,
        BuildOptions(
            backend="wolfram",
            wolfram_timeout=12.0,
            wolfram_factor_terms=True,
            wolfram_cse=True,
        ),
    )
    assert result == tmp_path
    assert RecordingKernel.events == [
        "start:None:12.0",
        "enter",
        "differentiate",
        "factor_terms",
        "wolfram_cse",
        "close",
    ]


def test_selected_wolfram_cse_failure_does_not_fallback(monkeypatch, tmp_path):
    RecordingKernel.events = []
    RecordingKernel.fail_cse = True
    monkeypatch.setattr(pipeline_module, "WolframKernel", RecordingKernel)
    monkeypatch.setattr(
        pipeline_module, "generate_matlab_bundle", _fake_generate
    )
    with pytest.raises(BuildError, match="No alternative strategy was attempted"):
        build_bundle(
            _config(),
            tmp_path,
            BuildOptions(backend="wolfram", wolfram_cse=True),
        )
    assert RecordingKernel.events.count("wolfram_cse") == 1
    assert RecordingKernel.events[-1] == "close"


def test_cli_rejects_wolfram_flags_with_sympy(capsys):
    assert cli_main([
        "build",
        "unused.toml",
        "--out",
        "unused",
        "--wolfram-cse",
    ]) == 2
    assert "require --backend wolfram" in capsys.readouterr().err


def test_factor_terms_error_guides_without_retry():
    error = BuildError(
        "code generation",
        "Wolfram FactorTerms",
        RuntimeError("Wolfram factor_terms failed: timed out"),
    )
    assert "No alternative strategy was attempted" in str(error)
    assert "--wolfram-factor-terms" in str(error)
