from __future__ import annotations

import json
import time
from dataclasses import replace
from pathlib import Path

import numpy as np
import pytest
import sympy as sp

from softarm.actuation import derive_actuation
from softarm.codegen import generate_matlab_bundle
from softarm.config import DynamicsConfig, IntegrationConfig, ModelConfig, load_config
from softarm.derive import derive
from softarm.models import RecursivePlant
from softarm.pipeline import BuildError, BuildOptions, build_bundle
from softarm.recursive import local_kernel_for
from softarm.special import LAMBDA_MODULES

ROOT = Path(__file__).parents[1]


def _recursive_config(segments: int = 1) -> ModelConfig:
    return ModelConfig(
        rod="euler_bernoulli",
        parameterization="pcs",
        segments=segments,
        dynamics=DynamicsConfig("recursive"),
        inertia="lumped",
        integration=IntegrationConfig(),
    )


def test_recursive_derivation_keeps_dynamics_algorithmic():
    plant = derive(_recursive_config(segments=20))
    assert isinstance(plant, RecursivePlant)
    assert plant.config.dynamics.formulation == "recursive"
    assert len(plant.arm_q) == 40
    assert not hasattr(plant, "mass")


def test_recursive_bundle_preserves_public_matlab_api(tmp_path):
    plant = derive(_recursive_config())
    generate_matlab_bundle(plant, tmp_path)
    for filename in (
        "softarm_mass.m",
        "softarm_bias.m",
        "softarm_inverse_dynamics.m",
        "softarm_forward_dynamics.m",
        "softarm_state_rhs.m",
        "softarm_kinematics.m",
        "softarm_centerline.m",
        "softarm_end_jacobian.m",
        "softarm_vehicle_wrench_map.m",
        "softarm_applied_force.m",
        "manifest.json",
        "softarm_model.tex",
    ):
        assert (tmp_path / filename).is_file()
    manifest = json.loads((tmp_path / "manifest.json").read_text(encoding="utf-8"))
    assert manifest["model"]["dynamics_formulation"] == "recursive"
    assert "softarm_inverse_dynamics" in (
        tmp_path / "softarm_mass.m"
    ).read_text(encoding="utf-8")
    wrapper = (tmp_path / "softarm_recursive_section.m").read_text(encoding="utf-8")
    assert "softarm_recursive_section_template" in wrapper


def test_twenty_section_bundle_reuses_one_local_template(tmp_path):
    start = time.perf_counter()
    plant = derive(_recursive_config(segments=20))
    generate_matlab_bundle(plant, tmp_path)
    elapsed = time.perf_counter() - start
    assert (tmp_path / "softarm_recursive_section_template.m").is_file()
    assert not (tmp_path / "softarm_recursive_section_20.m").exists()
    assert elapsed < 30.0


@pytest.mark.parametrize(
    ("name", "segments", "base_mode", "actuator_count"),
    [
        ("euler_bernoulli_pcs_recursive_n20", 20, "fixed", 0),
        ("extensible_euler_bernoulli_pcs_recursive_flying_n5", 5, "floating_rpy", 2),
    ],
)
def test_recursive_reference_bundle_manifest(name, segments, base_mode, actuator_count):
    config_path = ROOT / "examples" / "config" / f"{name}.toml"
    bundle_path = ROOT / "examples" / "generated" / name
    config = load_config(config_path)
    manifest = json.loads((bundle_path / "manifest.json").read_text(encoding="utf-8"))
    assert config.dynamics.formulation == "recursive"
    assert manifest["model"]["dynamics_formulation"] == "recursive"
    assert manifest["model"]["segments"] == segments
    assert manifest["model"]["base_mode"] == base_mode
    channels = [] if manifest["actuation"] is None else manifest["actuation"]["channels"]
    assert len(channels) == actuator_count
    assert (bundle_path / "softarm_inverse_dynamics.m").is_file()
    assert (bundle_path / "softarm_recursive_section_template.m").is_file()


def test_recursive_legacy_ritz_uses_affine_template(tmp_path):
    plant = derive(ModelConfig(
        rod="euler_bernoulli",
        parameterization="ritz",
        segments=2,
        dynamics=DynamicsConfig("recursive"),
        integration=IntegrationConfig(),
        ritz_x=(0.0, 0.0, 1.5, -0.5),
        ritz_y=(0.0, 0.0, 1.5, -0.5),
    ))
    generate_matlab_bundle(plant, tmp_path)
    assert (tmp_path / "softarm_legacy_affine_template.m").is_file()
    assert "softarm_legacy_step" in (
        tmp_path / "softarm_mass.m"
    ).read_text(encoding="utf-8")


def test_recursive_build_rejects_symbolic_only_options(tmp_path):
    config = _recursive_config()
    with pytest.raises(ValueError, match="backend='sympy'"):
        build_bundle(config, tmp_path / "wolfram", BuildOptions(backend="wolfram"))
    with pytest.raises(ValueError, match="--tex-appendix"):
        build_bundle(config, tmp_path / "tex", BuildOptions(tex_appendix=True))


def test_recursive_bundle_keeps_existing_actuator_interface(tmp_path):
    config = load_config(
        ROOT / "examples/config/euler_bernoulli_ritz_two_signed_pairs_n2.toml"
    )
    config = replace(config, dynamics=DynamicsConfig("recursive"))
    plant = derive(config)
    actuation = derive_actuation(plant)
    generate_matlab_bundle(plant, tmp_path, actuation=actuation)
    assert (tmp_path / "softarm_actuator_force.m").is_file()
    manifest = json.loads((tmp_path / "manifest.json").read_text(encoding="utf-8"))
    assert manifest["actuation"]["family"] == "tendon"


def test_recursive_constraint_request_fails_explicitly(tmp_path):
    config = load_config(
        ROOT
        / "examples/config/extensible_euler_bernoulli_pcs_flying_plane_contact_n1.toml"
    )
    config = replace(config, dynamics=DynamicsConfig("recursive"))
    with pytest.raises(BuildError, match="do not yet support constraints"):
        build_bundle(config, tmp_path)


def _parameter(plant, name):
    return next(item.symbol for item in plant.parameters if item.name == name)


def _fixed_base_recursive_expression(plant: RecursivePlant):
    definition = plant.definition
    kernel_builder = local_kernel_for(definition)
    ddq = sp.Matrix(sp.symbols(f"test_ddq1:{len(plant.q) + 1}", real=True))
    dof = len(plant.arm_q) // plant.config.segments
    gravity = sp.Matrix([0, 0, _parameter(plant, "gravity")])
    velocity = sp.zeros(3, 1)
    angular_velocity = sp.zeros(3, 1)
    acceleration = sp.zeros(3, 1)
    angular_acceleration = sp.zeros(3, 1)
    own_wrenches = []
    own_forces = []
    kernels = []
    for section in range(plant.config.segments):
        kernel = kernel_builder.build(definition, section)
        first = section * dof
        local_dq = plant.arm_dq[first : first + dof, 0]
        local_ddq = ddq[first : first + dof, 0]
        substitutions = dict(
            zip(kernel.acceleration_symbols, local_ddq, strict=True)
        )
        substitutions.update(
            zip(
                kernel.parent_velocity_symbols,
                velocity.col_join(angular_velocity),
                strict=True,
            )
        )
        substitutions.update(
            zip(
                kernel.parent_acceleration_symbols,
                acceleration.col_join(angular_acceleration),
                strict=True,
            )
        )
        substitutions.update(zip(kernel.gravity_symbols, gravity, strict=True))
        own_wrenches.append(kernel.own_wrench.subs(substitutions))
        own_forces.append(kernel.own_generalized_force.subs(substitutions))
        kernels.append(kernel)
        relative_v = kernel.end_linear_jacobian * local_dq
        relative_w = kernel.end_angular_jacobian * local_dq
        velocity_parent = (
            velocity
            + angular_velocity.cross(kernel.end_position)
            + relative_v
        )
        angular_velocity_parent = angular_velocity + relative_w
        acceleration_parent = (
            acceleration
            + angular_acceleration.cross(kernel.end_position)
            + angular_velocity.cross(
                angular_velocity.cross(kernel.end_position)
            )
            + 2 * angular_velocity.cross(relative_v)
            + kernel.end_linear_jacobian * local_ddq
            + kernel.end_linear_jacobian_rate * local_dq
        )
        angular_acceleration_parent = (
            angular_acceleration
            + angular_velocity.cross(relative_w)
            + kernel.end_angular_jacobian * local_ddq
            + kernel.end_angular_jacobian_rate * local_dq
        )
        inverse_rotation = kernel.end_rotation.T
        velocity = inverse_rotation * velocity_parent
        angular_velocity = inverse_rotation * angular_velocity_parent
        acceleration = inverse_rotation * acceleration_parent
        angular_acceleration = inverse_rotation * angular_acceleration_parent
        gravity = inverse_rotation * gravity

    tip_mass = _parameter(plant, "tip_mass")
    tip_inertia = sp.diag(
        _parameter(plant, "tip_Ixx"),
        _parameter(plant, "tip_Iyy"),
        _parameter(plant, "tip_Izz"),
    )
    tip_force = tip_mass * (acceleration - gravity)
    tip_moment = (
        tip_inertia * angular_acceleration
        + angular_velocity.cross(tip_inertia * angular_velocity)
    )
    wrench = tip_force.col_join(tip_moment)
    generalized = sp.zeros(len(plant.arm_q), 1)
    for section in reversed(range(plant.config.segments)):
        kernel = kernels[section]
        force = kernel.end_rotation * wrench[:3, 0]
        moment = kernel.end_rotation * wrench[3:, 0]
        first = section * dof
        generalized[first : first + dof, 0] = (
            own_forces[section]
            + kernel.end_linear_jacobian.T * force
            + kernel.end_angular_jacobian.T * moment
        )
        wrench = own_wrenches[section] + force.col_join(
            kernel.end_position.cross(force) + moment
        )
    damping = (
        sp.diag(*definition.damping)
        if isinstance(definition.damping, tuple)
        else definition.damping
    )
    generalized += damping * plant.arm_dq
    generalized += sp.Matrix(
        [sp.diff(definition.elastic, coordinate) for coordinate in plant.arm_q]
    )
    return generalized, ddq


@pytest.mark.parametrize(
    ("rod", "parameterization", "inertia", "integration", "ritz", "reference"),
    [
        ("euler_bernoulli", "ritz", "distributed", IntegrationConfig(), True, [0, 0]),
        ("euler_bernoulli", "pcs", "lumped", IntegrationConfig(), False, [0, 0]),
        (
            "extensible_euler_bernoulli",
            "ritz",
            "distributed",
            IntegrationConfig(),
            True,
            [0, 0, 0],
        ),
        (
            "extensible_euler_bernoulli",
            "pcs",
            "lumped",
            IntegrationConfig(),
            False,
            [0, 0, 0.5],
        ),
        ("euler_bernoulli", "pac", "lumped", IntegrationConfig(), False, [0, 0, 0]),
        (
            "extensible_euler_bernoulli",
            "pac",
            "lumped",
            IntegrationConfig(),
            False,
            [0, 0, 0, 0.5],
        ),
        ("cosserat", "pcs", "lumped", IntegrationConfig(), False, [0] * 6),
    ],
)
def test_recursive_inverse_dynamics_matches_symbolic_lagrange_at_reference(
    rod, parameterization, inertia, integration, ritz, reference
):
    common = dict(
        rod=rod,
        parameterization=parameterization,
        segments=1,
        inertia=inertia,
        integration=integration,
    )
    if ritz:
        common.update(
            ritz_x=(0.0, 0.0, 1.5, -0.5),
            ritz_y=(0.0, 0.0, 1.5, -0.5),
        )
        if rod == "extensible_euler_bernoulli":
            common["ritz_z"] = (0.0, 1.0)
    recursive = derive(ModelConfig(
        **common, dynamics=DynamicsConfig("recursive")
    ))
    symbolic = derive(ModelConfig(
        **common, dynamics=DynamicsConfig("symbolic_lagrange")
    ))
    expression, ddq_symbols = _fixed_base_recursive_expression(recursive)
    recursive_function = sp.lambdify(
        (recursive.q, recursive.dq, ddq_symbols, recursive.p),
        expression,
        [LAMBDA_MODULES, "numpy"],
    )
    symbolic_force = symbolic.mass * ddq_symbols
    if rod == "cosserat":
        symbolic_force += sp.Matrix(
            [sp.diff(symbolic.potential, coordinate) for coordinate in symbolic.q]
        ) + symbolic.damping * symbolic.dq
    else:
        symbolic_force += symbolic.bias
    symbolic_function = sp.lambdify(
        (symbolic.q, symbolic.dq, ddq_symbols, symbolic.p),
        symbolic_force,
        [LAMBDA_MODULES, "numpy"],
    )
    q = np.asarray(reference, dtype=float)
    dq = (
        np.zeros(len(q))
        if rod == "cosserat"
        else np.linspace(0.01, 0.02, len(q))
    )
    ddq = np.linspace(-0.03, 0.04, len(q))
    parameters = np.asarray([item.default for item in recursive.parameters])
    np.testing.assert_allclose(
        np.asarray(recursive_function(q, dq, ddq, parameters), dtype=float),
        np.asarray(symbolic_function(q, dq, ddq, parameters), dtype=float),
        rtol=2e-8,
        atol=2e-10,
    )
