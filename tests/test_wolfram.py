from contextlib import nullcontext
from pathlib import Path

import numpy as np
import pytest
import sympy as sp

from softarm import pipeline as pipeline_module
from softarm.backends.wolfram import WolframKernel
from softarm.config import IntegrationConfig, ModelConfig
from softarm.derive import derive
from softarm.dynamics import assemble_bias
from softarm.pipeline import BuildOptions, build_bundle
from softarm.special import LAMBDA_MODULES, SincSqrt, SincSqrtD


def test_wolfram_end_to_end(monkeypatch, tmp_path):
    local = Path(__file__).parents[1] / ".softarm.local.toml"
    if not local.is_file():
        pytest.skip("no local Wolfram tool configuration")
    config = ModelConfig(
        family="euler", segments=2, integration=IntegrationConfig("analytic"),
        ritz_x=(0.0, 0.0, 1.5, -0.5), ritz_y=(0.0, 0.0, 1.5, -0.5),
    )
    x = sp.Symbol("x", real=True, nonnegative=True)
    plant = derive(ModelConfig(
        family="pcc", segments=1, inertia="lumped",
        integration=IntegrationConfig(),
    ))
    with WolframKernel() as kernel:
        monkeypatch.setattr(
            pipeline_module,
            "WolframKernel",
            lambda *args, **kwargs: nullcontext(kernel),
        )
        build_bundle(config, tmp_path, BuildOptions(backend="wolfram"))
        assert kernel.differentiate([SincSqrt(x**2)], [x]) == [
            2 * x * SincSqrtD(x**2)
        ]
        factored = kernel.factor_terms([x * (x + 1)])
        assert sp.simplify(factored[0] - x * (x + 1)) == 0
        replacements, reduced = kernel.cse(
            [x**2 + sp.sin(x**2), x**2 + sp.cos(x**2)]
        )
        candidate = assemble_bias(plant, kernel)
    assert (tmp_path / "softarm_mass.m").is_file()
    assert (tmp_path / "softarm_bias.m").is_file()
    assert (tmp_path / "softarm_forward_dynamics.m").is_file()
    assert replacements
    rebuilt = list(reduced)
    for symbol, expression in reversed(replacements):
        rebuilt = [item.xreplace({symbol: expression}) for item in rebuilt]
    assert all(sp.simplify(left - right) == 0 for left, right in zip(
        rebuilt,
        [x**2 + sp.sin(x**2), x**2 + sp.cos(x**2)],
        strict=True,
    ))
    symbols = list(plant.q) + list(plant.dq) + [
        item.symbol for item in plant.parameters
    ]
    before = sp.lambdify(symbols, plant.bias, [LAMBDA_MODULES, "numpy"])
    after = sp.lambdify(symbols, candidate, [LAMBDA_MODULES, "numpy"])
    values = [0.03, -0.02, 0.51, 0.1, -0.05, 0.02]
    values += [item.default for item in plant.parameters]
    np.testing.assert_allclose(
        before(*values), after(*values), rtol=1e-10, atol=1e-12
    )
