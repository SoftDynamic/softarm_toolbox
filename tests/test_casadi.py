from __future__ import annotations

from dataclasses import replace
from pathlib import Path

import casadi as ca
import numpy as np
import pytest
import sympy as sp

from softarm import load_casadi_bundle
from softarm.actuation import derive_actuation
from softarm.codegen.casadi import SympyToCasadi, generate_casadi_bundle
from softarm.config import (
    ActuationConfig,
    ConstraintConfig,
    DynamicsConfig,
    IntegrationConfig,
    ModelConfig,
    TendonChannelConfig,
    TendonSpanConfig,
)
from softarm.constraints import derive_constraint
from softarm.derive import derive
from softarm.special import (
    AffineCosMoment,
    AffineSinMoment,
    CoscSqrt,
    CoscSqrtD,
    CoscSqrtDD,
    Sinc3Sqrt,
    Sinc3SqrtD,
    Sinc3SqrtDD,
    SincSqrt,
    SincSqrtD,
    SincSqrtDD,
    affine_cos_moment,
    affine_sin_moment,
)


def _ritz_config(formulation: str = "symbolic_lagrange") -> ModelConfig:
    return ModelConfig(
        rod="euler_bernoulli",
        parameterization="ritz",
        segments=1,
        dynamics=DynamicsConfig(formulation),
        integration=IntegrationConfig("analytic"),
        ritz_x=(0.0, 0.0, 1.5, -0.5),
        ritz_y=(0.0, 0.0, 1.5, -0.5),
    )


def _pcs_config(formulation: str) -> ModelConfig:
    return ModelConfig(
        rod="euler_bernoulli",
        parameterization="pcs",
        segments=2,
        dynamics=DynamicsConfig(formulation),
        inertia="lumped",
        integration=IntegrationConfig("analytic"),
    )


def test_special_functions_and_derivatives_are_finite_at_zero():
    z = sp.Symbol("z", real=True)
    expressions = sp.Matrix([
        SincSqrt(z),
        SincSqrtD(z),
        SincSqrtDD(z),
        CoscSqrt(z),
        CoscSqrtD(z),
        CoscSqrtDD(z),
        Sinc3Sqrt(z),
        Sinc3SqrtD(z),
        Sinc3SqrtDD(z),
    ])
    z_sx = ca.SX.sym("z")
    converted = SympyToCasadi({z: z_sx}, None).matrix(expressions)
    function = ca.Function("special", [z_sx], [converted, ca.jacobian(converted, z_sx)])
    values, derivatives = function(0.0)
    assert np.isfinite(np.asarray(values)).all()
    assert np.isfinite(np.asarray(derivatives)).all()
    np.testing.assert_allclose(
        np.asarray(values).ravel(),
        [1, -1 / 6, 1 / 60, 1 / 2, -1 / 24, 1 / 360, 1 / 6, -1 / 120, 1 / 2520],
        rtol=0,
        atol=1e-14,
    )


def test_affine_moment_mapaccum_matches_reference():
    c0_symbol, c1_symbol, xi_symbol = sp.symbols("c0 c1 xi", real=True)
    c0 = ca.SX.sym("c0")
    c1 = ca.SX.sym("c1")
    xi = ca.SX.sym("xi")
    converter = SympyToCasadi(
        {c0_symbol: c0, c1_symbol: c1, xi_symbol: xi}, 256
    )
    result = converter.matrix(sp.Matrix([
        AffineCosMoment(0, c0_symbol, c1_symbol, xi_symbol),
        AffineSinMoment(0, c0_symbol, c1_symbol, xi_symbol),
    ]))
    function = ca.Function("affine", [c0, c1, xi], [result])
    for values in ((0.2, -0.4, 1.0), (2 * np.pi, -2 * np.pi, 0.7)):
        actual = np.asarray(function(*values)).ravel()
        expected = [
            affine_cos_moment(0, *values),
            affine_sin_moment(0, *values),
        ]
        np.testing.assert_allclose(actual, expected, rtol=2e-13, atol=2e-13)


def test_converter_rejects_unknown_nodes_and_missing_pac_terms():
    x = sp.Symbol("x")
    x_sx = ca.SX.sym("x")
    with pytest.raises(TypeError, match="unsupported SymPy node"):
        SympyToCasadi({x: x_sx}, None).expression(sp.Function("Unknown")(x))
    with pytest.raises(ValueError, match="affine_terms"):
        SympyToCasadi({x: x_sx}, None).expression(AffineCosMoment(0, x, x, x))


def test_symbolic_bundle_round_trip_and_implicit_residual(tmp_path):
    plant = derive(_ritz_config())
    target = generate_casadi_bundle(plant, tmp_path)
    bundle = load_casadi_bundle(target)
    assert bundle.manifest["target"]["casadi_version"] == ca.__version__
    assert set(bundle.manifest["dynamics"]) == {"generalized_force"}
    explicit = bundle.function("dynamics.generalized_force.explicit")
    implicit = bundle.function("dynamics.generalized_force.implicit")
    x = np.zeros(explicit.size1_in("x"))
    u = np.zeros(explicit.size1_in("u"))
    p = np.asarray(bundle.parameters)
    xdot = np.asarray(explicit(x, u, p)).ravel()
    residual = np.asarray(implicit(xdot, x, u, np.zeros(0), p)).ravel()
    np.testing.assert_allclose(residual, 0.0, atol=1e-11)


def test_strict_tendon_bundle_exposes_kkt_solution_and_residual(tmp_path):
    config = _ritz_config()
    plant = derive(config)
    actuation = derive_actuation(
        plant,
        ActuationConfig(
            "tendon",
            "strict",
            (
                TendonChannelConfig(
                    "pair_x",
                    "signed",
                    (TendonSpanConfig(1, 0.02, 0.0),),
                ),
            ),
        ),
    )
    target = generate_casadi_bundle(plant, tmp_path, actuation=actuation)
    bundle = load_casadi_bundle(target)
    assert set(bundle.manifest["dynamics"]) == {
        "generalized_force",
        "tendon_force",
        "strict_tendon_acceleration",
    }
    solution = bundle.function("dynamics.strict_tendon_acceleration.solution")
    implicit = bundle.function("dynamics.strict_tendon_acceleration.implicit")
    q = np.zeros(len(plant.q))
    dq = np.zeros(len(plant.q))
    u = np.zeros(solution.size1_in("u"))
    p = np.asarray(bundle.parameters)
    ddq, tension = solution(q, dq, u, p)
    x = np.r_[q, dq]
    xdot = np.r_[dq, np.asarray(ddq).ravel()]
    residual = implicit(xdot, x, u, tension, p)
    np.testing.assert_allclose(np.asarray(residual), 0.0, atol=1e-10)


def test_active_constraint_bundle_exposes_kkt_solution_and_residual(tmp_path):
    config = replace(
        _ritz_config(),
        constraint=ConstraintConfig(
            "plane_point_contact",
            {
                "family": "plane_point_contact",
                "plane_normal": [1.0, 0.0, 0.0],
                "stabilization_frequency": 10.0,
                "stabilization_ratio": 1.0,
            },
        ),
    )
    plant = derive(config)
    constraint = derive_constraint(plant)
    bundle = load_casadi_bundle(
        generate_casadi_bundle(plant, tmp_path, constraint=constraint)
    )
    solution = bundle.function("dynamics.active_constraint.solution")
    implicit = bundle.function("dynamics.active_constraint.implicit")
    q = np.zeros(len(plant.q))
    dq = np.zeros(len(plant.q))
    u = np.zeros(solution.size1_in("u"))
    p = np.asarray(bundle.parameters)
    ddq, reaction = solution(q, dq, u, p)
    x = np.r_[q, dq]
    xdot = np.r_[dq, np.asarray(ddq).ravel()]
    residual = implicit(xdot, x, u, reaction, p)
    np.testing.assert_allclose(np.asarray(residual), 0.0, atol=1e-10)


def test_recursive_inverse_dynamics_matches_symbolic_lagrange(tmp_path):
    symbolic = derive(_pcs_config("symbolic_lagrange"))
    recursive = derive(_pcs_config("recursive"))
    symbolic_bundle = load_casadi_bundle(
        generate_casadi_bundle(symbolic, tmp_path / "symbolic")
    )
    recursive_bundle = load_casadi_bundle(
        generate_casadi_bundle(recursive, tmp_path / "recursive")
    )
    mass = symbolic_bundle.function("core.mass")
    bias = symbolic_bundle.function("core.bias")
    inverse = recursive_bundle.function("core.inverse_dynamics")
    rng = np.random.default_rng(4)
    q = 0.05 * rng.standard_normal(len(symbolic.q))
    dq = 0.04 * rng.standard_normal(len(symbolic.q))
    ddq = 0.03 * rng.standard_normal(len(symbolic.q))
    p = np.asarray(symbolic_bundle.parameters)
    expected = mass(q, p) @ ddq + bias(q, dq, p)
    np.testing.assert_allclose(inverse(q, dq, ddq, p), expected, rtol=2e-10, atol=2e-10)


def test_pac_terms_are_required_and_recorded(tmp_path):
    config = replace(_pcs_config("recursive"), parameterization="pac")
    # PAC uses three coordinates per section and requires Gauss integration.
    config = replace(
        config,
        inertia="distributed",
        integration=IntegrationConfig("gauss", 4),
    )
    plant = derive(config)
    with pytest.raises(ValueError, match="requires --casadi-affine-terms"):
        generate_casadi_bundle(plant, tmp_path / "missing")
    bundle = load_casadi_bundle(
        generate_casadi_bundle(plant, tmp_path / "pac", affine_terms=32)
    )
    assert bundle.manifest["target"]["affine_terms"] == 32


def test_generated_function_supports_c_code_generation(tmp_path, monkeypatch):
    plant = derive(_ritz_config())
    bundle = load_casadi_bundle(generate_casadi_bundle(plant, tmp_path / "bundle"))
    function = bundle.function("dynamics.generalized_force.explicit")
    output = tmp_path / "softarm_generated.c"
    monkeypatch.chdir(tmp_path)
    function.generate(output.name, {"with_header": True})
    assert output.is_file()
    assert output.with_suffix(".h").is_file()


def test_casadi_dynamics_use_linear_solves_without_explicit_inverse():
    source = (
        Path(__file__).parents[1] / "src" / "softarm" / "codegen" / "casadi.py"
    ).read_text(encoding="utf-8")
    assert "ca.solve(" in source
    assert ".inv(" not in source
    assert "sp.Inverse" not in source
