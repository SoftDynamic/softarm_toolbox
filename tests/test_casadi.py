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
    BaseConfig,
    ConstraintConfig,
    DynamicsConfig,
    IntegrationConfig,
    ModelConfig,
    TendonChannelConfig,
    TendonSpanConfig,
)
from softarm.constraints import derive_constraint, derive_constraint_definition
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


def _ritz_config(
    formulation: str = "symbolic_lagrange",
    *,
    rod: str = "euler_bernoulli",
    segments: int = 1,
) -> ModelConfig:
    return ModelConfig(
        rod=rod,
        parameterization="ritz",
        segments=segments,
        dynamics=DynamicsConfig(formulation),
        integration=IntegrationConfig("analytic"),
        ritz_x=(0.0, 0.0, 1.5, -0.5),
        ritz_y=(0.0, 0.0, 1.5, -0.5),
        ritz_z=(0.0, 1.0) if rod == "extensible_euler_bernoulli" else None,
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


def _assert_explicit_ad_matches_centered_difference(function, x, u, p):
    x_symbol = ca.MX.sym("x", len(x))
    u_symbol = ca.MX.sym("u", len(u))
    p_symbol = ca.MX.sym("p", len(p))
    expression = function(x_symbol, u_symbol, p_symbol)
    derivatives = ca.Function(
        "explicit_derivatives",
        [x_symbol, u_symbol, p_symbol],
        [ca.jacobian(expression, x_symbol), ca.jacobian(expression, u_symbol)],
    )
    ad_x, ad_u = derivatives(x, u, p)

    def centered(argument, evaluate):
        result = np.zeros((len(x), len(argument)))
        step = 2e-6
        for index in range(len(argument)):
            plus = np.array(argument, dtype=float)
            minus = np.array(argument, dtype=float)
            plus[index] += step
            minus[index] -= step
            result[:, index] = (
                np.asarray(evaluate(plus)).ravel()
                - np.asarray(evaluate(minus)).ravel()
            ) / (2 * step)
        return result

    finite_x = centered(x, lambda value: function(value, u, p))
    finite_u = centered(u, lambda value: function(x, value, p))
    np.testing.assert_allclose(ad_x, finite_x, rtol=2e-5, atol=2e-6)
    np.testing.assert_allclose(ad_u, finite_u, rtol=2e-5, atol=2e-6)


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
    _assert_explicit_ad_matches_centered_difference(
        bundle.function("dynamics.strict_tendon_acceleration.explicit"),
        x,
        u,
        p,
    )


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


@pytest.mark.parametrize(
    ("rod", "base_mode"),
    [
        ("euler_bernoulli", "fixed"),
        ("extensible_euler_bernoulli", "floating_rpy"),
    ],
)
def test_recursive_affine_ritz_matches_symbolic_lagrange(
    tmp_path, rod, base_mode
):
    base = BaseConfig(mode=base_mode)
    constraint_config = ConstraintConfig(
        "plane_point_contact",
        {"tool_offset": [0.01, -0.02, 0.03], "friction": 0.1},
    )
    symbolic_plant = derive(replace(
        _ritz_config("symbolic_lagrange", rod=rod, segments=2),
        base=base,
        constraint=constraint_config,
    ))
    recursive_plant = derive(replace(
        _ritz_config("recursive", rod=rod, segments=2),
        base=base,
        constraint=constraint_config,
    ))
    symbolic = load_casadi_bundle(generate_casadi_bundle(
        symbolic_plant,
        tmp_path / f"s-{base_mode}",
        constraint=derive_constraint(symbolic_plant),
    ))
    recursive = load_casadi_bundle(generate_casadi_bundle(
        recursive_plant,
        tmp_path / f"r-{base_mode}",
        constraint=constraint_config,
    ))
    rng = np.random.default_rng(17)
    q = 0.04 * rng.standard_normal(len(recursive_plant.q))
    if base_mode == "floating_rpy":
        q[3:6] = [0.12, -0.08, 0.16]
    dq = 0.03 * rng.standard_normal(len(q))
    ddq = 0.02 * rng.standard_normal(len(q))
    p = np.asarray(recursive.parameters)
    expected = (
        symbolic.function("core.mass")(q, p) @ ddq
        + symbolic.function("core.bias")(q, dq, p)
    )
    actual = recursive.function("core.inverse_dynamics")(q, dq, ddq, p)
    np.testing.assert_allclose(actual, expected, rtol=2e-8, atol=2e-9)
    np.testing.assert_allclose(
        recursive.function("core.kinematics")(q, p),
        symbolic.function("core.kinematics")(q, p),
        rtol=2e-10,
        atol=2e-10,
    )
    for logical, inputs in (
        ("core.constraint_value", (q, p)),
        ("core.constraint_jacobian", (q, p)),
        ("core.constraint_velocity_bias", (q, dq, p)),
        ("core.constraint_reaction_map", (q, dq, p)),
    ):
        np.testing.assert_allclose(
            recursive.function(logical)(*inputs),
            symbolic.function(logical)(*inputs),
            rtol=3e-8,
            atol=3e-9,
        )


def test_recursive_affine_ritz_derivative_structure(tmp_path):
    plant = derive(replace(
        _ritz_config("recursive", segments=2),
        base=BaseConfig(mode="floating_rpy"),
    ))
    bundle = load_casadi_bundle(generate_casadi_bundle(plant, tmp_path))
    kinematics = bundle.function("core.kinematics")
    q = ca.MX.sym("q", len(plant.q))
    p = ca.MX.sym("p", len(bundle.parameters))
    flat = ca.reshape(kinematics(q, p), -1, 1)
    first = ca.jacobian(flat, q)
    nbase = len(plant.base_q)
    second = ca.jacobian(first[:, nbase], q)
    arm_arm = second[:, nbase + 1]
    base_arm = second[:, 3]
    derivative = ca.Function(
        "affine_derivative_structure", [q, p], [arm_arm, base_arm]
    )
    state = np.zeros(len(plant.q))
    state[3:6] = [0.17, -0.11, 0.09]
    state[nbase:] = np.linspace(0.01, 0.04, len(plant.arm_q))
    arm_arm_value, base_arm_value = derivative(state, bundle.parameters)
    np.testing.assert_allclose(np.asarray(arm_arm_value), 0.0, atol=2e-12)
    assert np.max(np.abs(np.asarray(base_arm_value))) > 1e-5


def test_recursive_active_constraint_solution_satisfies_implicit_residual(tmp_path):
    config = replace(
        _pcs_config("recursive"),
        base=BaseConfig(mode="floating_rpy"),
        constraint=ConstraintConfig(
            "plane_point_contact",
            {
                "plane_normal": [1.0, 0.2, -0.1],
                "tool_offset": [0.02, -0.01, 0.03],
                "friction": 0.15,
            },
        ),
    )
    plant = derive(config)
    constraint = derive_constraint_definition(config.constraint)
    bundle = load_casadi_bundle(
        generate_casadi_bundle(plant, tmp_path, constraint=constraint)
    )
    solution = bundle.function("dynamics.active_constraint.solution")
    implicit = bundle.function("dynamics.active_constraint.implicit")
    rng = np.random.default_rng(31)
    q = 0.03 * rng.standard_normal(len(plant.q))
    q[3:6] = [0.10, -0.06, 0.13]
    dq = 0.02 * rng.standard_normal(len(plant.q))
    u = np.zeros(len(plant.arm_q) + 6 + 6 + 1)
    p = np.asarray(bundle.parameters)
    ddq, reaction = solution(q, dq, u, p)
    x = np.r_[q, dq]
    xdot = np.r_[dq, np.asarray(ddq).ravel()]
    residual = implicit(xdot, x, u, reaction, p)
    np.testing.assert_allclose(np.asarray(residual), 0.0, atol=2e-9)
    _assert_explicit_ad_matches_centered_difference(
        bundle.function("dynamics.active_constraint.explicit"), x, u, p
    )


def test_recursive_constraint_terms_match_symbolic_on_floating_base(tmp_path):
    constraint_config = ConstraintConfig(
        "plane_point_contact",
        {
            "plane_normal": [0.3, -0.2, -0.9],
            "plane_point": [0.02, -0.03, 0.04],
            "tool_offset": [0.015, -0.01, 0.025],
            "friction": 0.2,
        },
    )
    common = replace(
        _pcs_config("symbolic_lagrange"),
        segments=1,
        base=BaseConfig(mode="floating_rpy"),
        constraint=constraint_config,
    )
    symbolic_plant = derive(common)
    recursive_plant = derive(replace(
        common, dynamics=DynamicsConfig("recursive")
    ))
    symbolic = load_casadi_bundle(generate_casadi_bundle(
        symbolic_plant,
        tmp_path / "s",
        constraint=derive_constraint(symbolic_plant),
    ))
    recursive = load_casadi_bundle(generate_casadi_bundle(
        recursive_plant,
        tmp_path / "r",
        constraint=constraint_config,
    ))
    rng = np.random.default_rng(23)
    q = 0.04 * rng.standard_normal(len(recursive_plant.q))
    q[3:6] = [0.11, -0.07, 0.14]
    dq = 0.03 * rng.standard_normal(len(q))
    p = np.asarray(recursive.parameters)
    for logical, inputs in (
        ("core.constraint_value", (q, p)),
        ("core.constraint_jacobian", (q, p)),
        ("core.constraint_velocity_bias", (q, dq, p)),
        ("core.constraint_reaction_map", (q, dq, p)),
    ):
        np.testing.assert_allclose(
            recursive.function(logical)(*inputs),
            symbolic.function(logical)(*inputs),
            rtol=3e-9,
            atol=3e-10,
        )


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


def test_recursive_functions_support_mx_composition_and_c_generation(
    tmp_path, monkeypatch
):
    ritz_plant = derive(replace(
        _ritz_config("recursive", segments=2),
        base=BaseConfig(mode="floating_rpy"),
    ))
    contact_config = replace(
        _pcs_config("recursive"),
        base=BaseConfig(mode="floating_rpy"),
        constraint=ConstraintConfig("plane_point_contact", {"friction": 0.2}),
    )
    contact_plant = derive(contact_config)
    bundles = (
        (
            "ritz",
            load_casadi_bundle(generate_casadi_bundle(
                ritz_plant, tmp_path / "ritz"
            )),
            "dynamics.generalized_force.explicit",
        ),
        (
            "contact",
            load_casadi_bundle(generate_casadi_bundle(
                contact_plant,
                tmp_path / "contact",
                constraint=contact_config.constraint,
            )),
            "dynamics.active_constraint.explicit",
        ),
    )
    monkeypatch.chdir(tmp_path)
    for name, bundle, logical in bundles:
        function = bundle.function(logical)
        x = ca.MX.sym(f"x_{name}", function.size1_in("x"))
        u = ca.MX.sym(f"u_{name}", function.size1_in("u"))
        p = ca.MX.sym(f"p_{name}", function.size1_in("p"))
        expression = function(x, u, p)
        assert expression.shape == (function.size1_out(0), 1)
        output = tmp_path / f"softarm_{name}.c"
        function.generate(output.name, {"with_header": True})
        assert output.is_file()
        assert output.with_suffix(".h").is_file()
    ritz_document = (tmp_path / "ritz" / "softarm_model.tex").read_text(
        encoding="utf-8"
    )
    contact_document = (tmp_path / "contact" / "softarm_model.tex").read_text(
        encoding="utf-8"
    )
    assert "global first-order affine" in ritz_document
    assert "Exact-SE(3)" in contact_document


def test_casadi_dynamics_use_linear_solves_without_explicit_inverse():
    source = (
        Path(__file__).parents[1] / "src" / "softarm" / "codegen" / "casadi.py"
    ).read_text(encoding="utf-8")
    assert "ca.solve(" in source
    assert ".inv(" not in source
    assert "sp.Inverse" not in source
