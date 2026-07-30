from pathlib import Path

import numpy as np
import pytest
import sympy as sp

from softarm.backends.wolfram_backend import optimize
from softarm.codegen import generate_matlab_bundle
from softarm.config import IntegrationConfig, ModelConfig
from softarm.derive import derive
from softarm.special import LAMBDA_MODULES


def test_wolfram_n3_kinematics_round_trip():
    local = Path(__file__).parents[1] / ".softarm.local.toml"
    if not local.is_file():
        pytest.skip("no local Wolfram tool configuration")
    plant = derive(ModelConfig(family="pcc", segments=3, inertia="lumped", integration=IntegrationConfig()))
    original = [plant.kinematics[row, column] for column in range(plant.kinematics.cols) for row in range(4)]
    transformed = optimize(original)
    symbols = list(plant.q) + [item.symbol for item in plant.parameters]
    before = sp.lambdify(symbols, original, [LAMBDA_MODULES, "numpy"])
    after = sp.lambdify(symbols, transformed, [LAMBDA_MODULES, "numpy"])
    q = [0.03, -0.02, 0.50, -0.01, 0.04, 0.48, 0.02, 0.01, 0.52]
    values = q + [item.default for item in plant.parameters]
    np.testing.assert_allclose(before(*values), after(*values), rtol=1e-10, atol=1e-12)


def test_wolfram_n3_full_euler_generation(tmp_path):
    local = Path(__file__).parents[1] / ".softarm.local.toml"
    if not local.is_file():
        pytest.skip("no local Wolfram tool configuration")
    plant = derive(ModelConfig(
        family="euler", segments=3, integration=IntegrationConfig("analytic"),
        ritz_x=(0.0, 0.0, 1.5, -0.5), ritz_y=(0.0, 0.0, 1.5, -0.5),
    ))
    generate_matlab_bundle(plant, tmp_path, backend="wolfram")
    assert (tmp_path / "softarm_mass.m").is_file()
    assert (tmp_path / "softarm_bias.m").is_file()
    assert (tmp_path / "softarm_forward_dynamics.m").is_file()
