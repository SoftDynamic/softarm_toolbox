from dataclasses import replace
from pathlib import Path

import numpy as np
import pytest
import sympy as sp

from softarm.config import BaseConfig, IntegrationConfig, ModelConfig, load_config
from softarm.derive import derive, register_model
from softarm.special import LAMBDA_MODULES

ROOT = Path(__file__).parents[1]


def test_external_combination_registration_runs_its_validator():
    events: list[str] = []

    def validator(config):
        events.append(f"validate:{config.rod}:{config.parameterization}")

    def builder(config):
        return derive(replace(
            config,
            rod="euler_bernoulli",
            parameterization="pcs",
            inertia="lumped",
        ))

    register_model("custom_rod", "custom_parameterization", builder, validator=validator)
    config = ModelConfig(
        rod="custom_rod", parameterization="custom_parameterization", segments=1,
        inertia="lumped",
    )
    assert derive(config).config == config
    assert events == ["validate:custom_rod:custom_parameterization"]
    with pytest.raises(ValueError, match="already registered"):
        register_model("custom_rod", "custom_parameterization", builder)


def _numeric(matrix, plant, q):
    fn = sp.lambdify((plant.q, plant.p), matrix, [LAMBDA_MODULES, "numpy"])
    p = np.array([item.default for item in plant.parameters])
    return np.asarray(fn(np.asarray(q), p), dtype=float)


def test_lumped_extensible_kirchhoff_pcs_mass_is_symmetric_and_finite():
    plant = derive(load_config(ROOT / "examples/config/extensible_kirchhoff_pcs_lumped_n2.toml"))
    q = [0.08, -0.04, 0.46, -0.03, 0.06, 0.51]
    mass = _numeric(plant.mass, plant, q)
    np.testing.assert_allclose(mass, mass.T, rtol=1e-11, atol=1e-12)
    assert np.isfinite(mass).all()
    assert np.linalg.eigvalsh(mass).min() > 0


def test_euler_shapes_and_ritz_stiffness():
    plant = derive(load_config(ROOT / "examples/config/euler_bernoulli_ritz_n2.toml"))
    assert plant.mass.shape == (4, 4)
    assert plant.bias.shape == (4, 1)
    assert plant.kinematics.shape == (4, 8)
    assert plant.end_jacobian.shape == (6, 4)
    # Integral of (d2/dxi2 (1.5 xi^2 - .5 xi^3))^2 on [0,1] is 3.
    assert sp.integrate((3 - 3 * sp.Symbol("xi")) ** 2, (sp.Symbol("xi"), 0, 1)) == 3


def test_euler_bernoulli_pcs_has_only_fixed_length_bending_coordinates():
    plant = derive(load_config(ROOT / "examples/config/euler_bernoulli_pcs_n2.toml"))
    assert plant.arm_coordinate_names == ["bx1", "by1", "bx2", "by2"]
    assert plant.mass.shape == (4, 4)
    defaults = {item.symbol: item.default for item in plant.parameters}
    straight = plant.end_transform.subs(defaults).subs({item: 0 for item in plant.q})
    assert np.isclose(float(straight[2, 3]), 0.95)


def test_extensible_kirchhoff_ritz_has_three_displacement_coordinates():
    plant = derive(load_config(ROOT / "examples/config/extensible_kirchhoff_ritz_n2.toml"))
    assert plant.arm_coordinate_names == ["ax1", "ay1", "az1", "ax2", "ay2", "az2"]
    assert plant.mass.shape == (6, 6)
    axial_stiffness = sp.diff(plant.potential, plant.arm_q[2], 2)
    defaults = {item.symbol: item.default for item in plant.parameters}
    assert np.isclose(float(axial_stiffness.subs(defaults)), 100.0)


def test_floating_base_has_coupled_coordinates_and_wrench_map():
    fixed = derive(ModelConfig(
        rod="extensible_kirchhoff", parameterization="pcs", segments=1,
        inertia="lumped", integration=IntegrationConfig()
    ))
    floating = derive(ModelConfig(
        rod="extensible_kirchhoff", parameterization="pcs", segments=1,
        inertia="lumped", integration=IntegrationConfig(),
        base=BaseConfig("floating_rpy"),
    ))
    assert floating.base_coordinate_names == [
        "base_x", "base_y", "base_z", "base_roll", "base_pitch", "base_yaw"
    ]
    assert floating.arm_coordinate_names == fixed.arm_coordinate_names
    assert floating.mass.shape == (9, 9)
    assert floating.end_jacobian.shape == (6, 9)
    assert floating.vehicle_wrench_map.shape == (9, 6)
    q = np.zeros(len(floating.q))
    q[8] = 0.5
    p = np.array([item.default for item in floating.parameters])
    evaluate = sp.lambdify(
        (floating.q, floating.p), floating.mass, [LAMBDA_MODULES, "numpy"]
    )
    mass = np.asarray(evaluate(q, p), dtype=float)
    assert np.linalg.norm(mass[:6, 6:]) > 0
    np.testing.assert_allclose(mass, mass.T, atol=1e-12)
    assert floating.mass.diff(floating.base_q[0]) == sp.zeros(9)
    assert floating.mass.diff(floating.base_q[1]) == sp.zeros(9)
    assert floating.mass.diff(floating.base_q[2]) == sp.zeros(9)

    q_arm = np.array([0.02, -0.01, 0.5])
    fixed_p = np.array([item.default for item in fixed.parameters])
    fixed_mass = sp.lambdify(
        (fixed.q, fixed.p), fixed.mass, [LAMBDA_MODULES, "numpy"]
    )(q_arm, fixed_p)
    q[6:] = q_arm
    floating_mass = np.asarray(evaluate(q, p), dtype=float)
    np.testing.assert_allclose(floating_mass[6:, 6:], fixed_mass, rtol=1e-10, atol=1e-12)


def test_cosserat_pcs_lumped_shapes_energy_and_nominal_mass():
    plant = derive(ModelConfig(
        rod="cosserat", parameterization="pcs", segments=1, inertia="lumped",
        integration=IntegrationConfig(),
    ))
    assert plant.arm_coordinate_names == ["kx1", "ky1", "kz1", "vx1", "vy1", "vz1"]
    assert plant.mass.shape == (6, 6)
    assert plant.kinematics.shape == (4, 4)
    assert plant.end_jacobian.shape == (6, 6)
    defaults = {item.symbol: item.default for item in plant.parameters}
    reference = {coordinate: 0.0 for coordinate in plant.q}
    mass_function = sp.lambdify(
        (plant.q, plant.p), plant.mass, [LAMBDA_MODULES, "numpy"]
    )
    mass = np.asarray(
        mass_function(
            np.zeros(len(plant.q)),
            np.array([item.default for item in plant.parameters]),
        ),
        dtype=float,
    )
    np.testing.assert_allclose(mass, mass.T, atol=1e-13)
    assert np.linalg.eigvalsh(mass).min() > 0
    elastic_hessian = sp.hessian(
        plant.potential.subs({sp.Symbol("gravity", real=True): 0}), plant.arm_q
    )
    expected = np.diag([0.6, 0.6, 0.1, 10.0, 10.0, 25.0])
    np.testing.assert_allclose(
        np.asarray(elastic_hessian.subs(defaults).subs(reference), dtype=float),
        expected,
        atol=1e-13,
    )
