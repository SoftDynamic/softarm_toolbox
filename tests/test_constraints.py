import numpy as np
import sympy as sp

from softarm.config import BaseConfig, ConstraintConfig, IntegrationConfig, ModelConfig
from softarm.constraints import derive_constraint
from softarm.derive import derive
from softarm.special import LAMBDA_MODULES


def _contact(friction: float = 0.3):
    config = ModelConfig(
        rod="extensible_euler_bernoulli",
        parameterization="pcs",
        segments=1,
        inertia="lumped",
        integration=IntegrationConfig(),
        base=BaseConfig("floating_rpy"),
        constraint=ConstraintConfig("plane_point_contact", {
            "family": "plane_point_contact",
            "plane_point": [0.0, 0.0, 0.5],
            "plane_normal": [0.0, 0.0, -1.0],
            "friction": friction,
            "friction_velocity": 0.01,
        }),
    )
    plant = derive(config)
    constraint = derive_constraint(plant)
    assert constraint is not None
    return plant, constraint


def test_plane_contact_shapes_and_constraint_bias():
    plant, constraint = _contact()
    assert constraint.coordinates.shape == (1, 1)
    assert constraint.jacobian.shape == (1, len(plant.q))
    assert constraint.reaction_map.shape == (len(plant.q), 1)
    parameters = plant.parameters + constraint.parameters
    p = np.array([item.default for item in parameters])
    q = np.zeros(len(plant.q))
    q[8] = 0.5
    dq = np.zeros(len(plant.q))
    value = sp.lambdify(
        (plant.q, sp.Matrix([item.symbol for item in parameters])),
        constraint.coordinates, [LAMBDA_MODULES, "numpy"],
    )
    bias = sp.lambdify(
        (plant.q, plant.dq, sp.Matrix([item.symbol for item in parameters])),
        constraint.velocity_bias, [LAMBDA_MODULES, "numpy"],
    )
    assert abs(float(np.asarray(value(q, p)).item())) < 1e-12
    assert abs(float(np.asarray(bias(q, dq, p)).item())) < 1e-12
    step = 1e-7
    q_plus = q.copy()
    q_plus[2] += step
    q_minus = q.copy()
    q_minus[2] -= step
    numerical = (np.asarray(value(q_plus, p))-np.asarray(value(q_minus, p)))/(2*step)
    jacobian = sp.lambdify(
        (plant.q, sp.Matrix([item.symbol for item in parameters])),
        constraint.jacobian, [LAMBDA_MODULES, "numpy"],
    )
    np.testing.assert_allclose(numerical.item(), jacobian(q, p)[0, 2], rtol=1e-8)
    assert constraint.is_feasible(sp.Matrix([1.0]))
    assert not constraint.is_feasible(sp.Matrix([-1.0]))


def test_plane_contact_friction_is_dissipative_and_mu_zero_is_normal():
    for friction in (0.0, 0.4):
        plant, constraint = _contact(friction)
        parameters = plant.parameters + constraint.parameters
        parameter_symbols = sp.Matrix([item.symbol for item in parameters])
        p = np.array([item.default for item in parameters])
        q = np.zeros(len(plant.q))
        q[8] = 0.5
        dq = np.zeros(len(plant.q))
        dq[0] = 0.2
        power = sp.lambdify(
            (plant.q, plant.dq, parameter_symbols),
            (plant.dq.T * constraint.reaction_map)[0],
            [LAMBDA_MODULES, "numpy"],
        )
        actual = float(power(q, dq, p))
        assert actual <= 1e-12
        if friction == 0.0:
            expected = constraint.jacobian.T
            compare = sp.lambdify(
                (plant.q, plant.dq, parameter_symbols),
                constraint.reaction_map - expected,
                [LAMBDA_MODULES, "numpy"],
            )
            np.testing.assert_allclose(compare(q, dq, p), 0.0, atol=1e-12)
