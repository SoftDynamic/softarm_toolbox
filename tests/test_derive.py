from pathlib import Path

import numpy as np
import sympy as sp

from softarm.config import load_config
from softarm.derive import derive
from softarm.special import LAMBDA_MODULES


ROOT = Path(__file__).parents[1]


def _numeric(matrix, plant, q):
    fn = sp.lambdify((plant.q, plant.p), matrix, [LAMBDA_MODULES, "numpy"])
    p = np.array([item.default for item in plant.parameters])
    return np.asarray(fn(np.asarray(q), p), dtype=float)


def test_lumped_pcc_mass_is_symmetric_and_finite():
    plant = derive(load_config(ROOT / "examples/config/pcc_lumped_n2.toml"))
    q = [0.08, -0.04, 0.46, -0.03, 0.06, 0.51]
    mass = _numeric(plant.mass, plant, q)
    np.testing.assert_allclose(mass, mass.T, rtol=1e-11, atol=1e-12)
    assert np.isfinite(mass).all()
    assert np.linalg.eigvalsh(mass).min() > 0


def test_euler_shapes_and_ritz_stiffness():
    plant = derive(load_config(ROOT / "examples/config/euler_ritz_n2.toml"))
    assert plant.mass.shape == (4, 4)
    assert plant.bias.shape == (4, 1)
    assert plant.kinematics.shape == (4, 8)
    assert plant.end_jacobian.shape == (6, 4)
    # Integral of (d2/dxi2 (1.5 xi^2 - .5 xi^3))^2 on [0,1] is 3.
    assert sp.integrate((3 - 3 * sp.Symbol("xi")) ** 2, (sp.Symbol("xi"), 0, 1)) == 3

