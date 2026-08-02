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


@pytest.mark.parametrize("config_name", [
    "euler_bernoulli_ritz_n2.toml",
    "euler_bernoulli_pcs_n2.toml",
    "extensible_euler_bernoulli_ritz_n2.toml",
    "extensible_euler_bernoulli_pcs_lumped_n2.toml",
    "euler_bernoulli_pac_distributed_n2.toml",
    "extensible_euler_bernoulli_pac_distributed_n2.toml",
    "cosserat_pcs_lumped_n1.toml",
])
def test_material_kinematics_matches_section_ends(config_name):
    plant = derive(load_config(ROOT / "examples/config" / config_name))
    assert plant._material_coordinate is not None
    assert plant._material_kinematics is not None
    q = np.linspace(0.01, 0.02, len(plant.q))
    p = np.array([item.default for item in plant.parameters])
    evaluate_material = sp.lambdify(
        (plant.q, plant.p, plant._material_coordinate),
        plant._material_kinematics,
        [LAMBDA_MODULES, "numpy"],
    )
    material = np.asarray(evaluate_material(q, p, 0.37), dtype=float)
    material_end = np.asarray(evaluate_material(q, p, 1.0), dtype=float)
    endpoint = _numeric(plant.kinematics, plant, q)
    assert np.isfinite(material).all()
    np.testing.assert_allclose(material_end, endpoint, rtol=1e-11, atol=1e-12)


def test_material_kinematics_respects_floating_base_and_mount():
    mount_xyz = (0.12, -0.23, 0.34)
    mount_rpy = (0.17, -0.11, 0.08)
    plant = derive(ModelConfig(
        rod="extensible_euler_bernoulli", parameterization="pcs", segments=1,
        inertia="lumped", integration=IntegrationConfig(),
        base=BaseConfig("floating_rpy", mount_xyz, mount_rpy),
    ))
    q = np.array([1.1, -2.2, 3.3, 0.21, -0.31, 0.42, 0.0, 0.0, 0.5])
    p = np.array([item.default for item in plant.parameters])
    evaluate = sp.lambdify(
        (plant.q, plant.p, plant._material_coordinate),
        plant._material_kinematics,
        [LAMBDA_MODULES, "numpy"],
    )
    actual = np.asarray(evaluate(q, p, 0.4), dtype=float)[:3, 3]

    def rotation(rpy):
        roll, pitch, yaw = rpy
        cx, sx = np.cos(roll), np.sin(roll)
        cy, sy = np.cos(pitch), np.sin(pitch)
        cz, sz = np.cos(yaw), np.sin(yaw)
        return np.array([
            [cz * cy, cz * sy * sx - sz * cx, cz * sy * cx + sz * sx],
            [sz * cy, sz * sy * sx + cz * cx, sz * sy * cx - cz * sx],
            [-sy, cy * sx, cy * cx],
        ])

    base_rotation = rotation(q[3:6])
    expected = q[:3] + base_rotation @ (
        np.asarray(mount_xyz) + rotation(mount_rpy) @ np.array([0.0, 0.0, 0.2])
    )
    np.testing.assert_allclose(actual, expected, rtol=1e-11, atol=1e-12)


def test_lumped_extensible_euler_bernoulli_pcs_mass_is_symmetric_and_finite():
    plant = derive(
        load_config(ROOT / "examples/config/extensible_euler_bernoulli_pcs_lumped_n2.toml")
    )
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


def test_extensible_euler_bernoulli_ritz_has_three_displacement_coordinates():
    plant = derive(
        load_config(ROOT / "examples/config/extensible_euler_bernoulli_ritz_n2.toml")
    )
    assert plant.arm_coordinate_names == ["ax1", "ay1", "az1", "ax2", "ay2", "az2"]
    assert plant.mass.shape == (6, 6)
    axial_stiffness = sp.diff(plant.potential, plant.arm_q[2], 2)
    defaults = {item.symbol: item.default for item in plant.parameters}
    assert np.isclose(float(axial_stiffness.subs(defaults)), 100.0)


def test_pac_coordinates_reference_lengths_and_anisotropic_hankel_energy():
    euler = derive(ModelConfig(
        rod="euler_bernoulli", parameterization="pac", segments=1,
        inertia="lumped", integration=IntegrationConfig(),
        parameters={"length": 0.5, "EI_x": 2.0, "EI_y": 3.0, "GJ": 0.4},
    ))
    assert euler.arm_coordinate_names == ["c0_1", "c1_1", "phi1"]
    defaults = {item.symbol: item.default for item in euler.parameters}
    elastic = euler.potential.subs(sp.Symbol("gravity", real=True), 0)
    hessian = sp.hessian(elastic, euler.arm_q)
    expected = np.array([[6.0, 3.0, 0.0], [3.0, 2.0, 0.0], [0.0, 0.0, 0.8]])
    np.testing.assert_allclose(
        np.asarray(hessian.subs(defaults).subs({item: 0 for item in euler.arm_q}), dtype=float),
        expected,
        atol=1e-13,
    )

    extensible = derive(ModelConfig(
        rod="extensible_euler_bernoulli", parameterization="pac", segments=1,
        inertia="lumped", integration=IntegrationConfig(),
    ))
    assert extensible.arm_coordinate_names == ["c0_1", "c1_1", "phi1", "l1"]
    reference = [0.0, 0.0, 0.0, 0.5]
    mass = _numeric(extensible.mass, extensible, reference)
    np.testing.assert_allclose(mass, mass.T, atol=1e-13)
    assert np.linalg.eigvalsh(mass).min() > 0
    transform = _numeric(extensible.end_transform, extensible, reference)
    np.testing.assert_allclose(transform[:3, 3], [0.0, 0.0, 0.5], atol=1e-14)


def test_floating_base_has_coupled_coordinates_and_wrench_map():
    fixed = derive(ModelConfig(
        rod="extensible_euler_bernoulli", parameterization="pcs", segments=1,
        inertia="lumped", integration=IntegrationConfig()
    ))
    floating = derive(ModelConfig(
        rod="extensible_euler_bernoulli", parameterization="pcs", segments=1,
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
