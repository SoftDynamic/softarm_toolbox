import math

import numpy as np
import sympy as sp

from softarm.config import IntegrationConfig, ModelConfig
from softarm.derive import derive
from softarm.special import LAMBDA_MODULES


def _pcc_position(q, xi):
    bx, by, length = q
    x = xi * bx
    y = xi * by
    z = x*x + y*y
    if abs(z) < 1e-12:
        a = 1-z/6+z*z/120
        b = 0.5-z/24+z*z/720
    else:
        root = math.sqrt(z)
        a = math.sin(root)/root
        b = (1-math.cos(root))/z
    return np.array([length*xi*x*b, length*xi*y*b, length*xi*a])


def test_distributed_mass_matches_independent_numerical_quadrature():
    config = ModelConfig(
        family="pcc", segments=1, inertia="distributed", integration=IntegrationConfig("gauss", 4),
        parameters={
            "mass": 1.0, "Ixx": 0.0, "Iyy": 0.0, "Izz": 0.0,
            "tip_mass": 0.0, "tip_Ixx": 0.0, "tip_Iyy": 0.0, "tip_Izz": 0.0,
        },
    )
    plant = derive(config)
    q = np.array([0.12, -0.07, 0.53])
    p = np.array([item.default for item in plant.parameters])
    symbolic = np.asarray(sp.lambdify((plant.q, plant.p), plant.mass, [LAMBDA_MODULES, "numpy"])(q, p), dtype=float)

    nodes, weights = np.polynomial.legendre.leggauss(30)
    independent = np.zeros((3, 3))
    step = 2e-6
    for node, weight in zip((nodes+1)/2, weights/2, strict=True):
        jacobian = np.zeros((3, 3))
        for column in range(3):
            delta = np.zeros(3); delta[column] = step
            jacobian[:, column] = (_pcc_position(q+delta, node)-_pcc_position(q-delta, node))/(2*step)
        independent += weight * jacobian.T @ jacobian
    np.testing.assert_allclose(symbolic, independent, rtol=2e-5, atol=2e-7)


def test_coriolis_term_satisfies_energy_identity():
    config = ModelConfig(family="pcc", segments=1, inertia="lumped", integration=IntegrationConfig())
    plant = derive(config)
    q = np.array([0.09, -0.04, 0.51])
    dq = np.array([0.03, -0.02, 0.01])
    p = np.array([item.default for item in plant.parameters])
    conservative = sp.Matrix([plant.potential]).jacobian(plant.q).T
    coriolis = plant.bias - conservative - plant.damping*plant.dq
    mdot = sp.zeros(len(plant.q))
    for index, coordinate in enumerate(plant.q):
        mdot += plant.mass.diff(coordinate) * plant.dq[index]
    residual = (plant.dq.T * (coriolis-sp.Rational(1,2)*mdot*plant.dq))[0]
    evaluate = sp.lambdify((plant.q, plant.dq, plant.p), residual, [LAMBDA_MODULES, "numpy"])
    assert abs(float(evaluate(q, dq, p))) < 1e-10
