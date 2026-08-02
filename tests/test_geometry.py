import numpy as np
import sympy as sp

from softarm.geometry import pac_transform, pcs_transform
from softarm.special import (
    LAMBDA_MODULES,
    AffineCosMoment,
    AffineSinMoment,
    Sinc3Sqrt,
    Sinc3SqrtD,
    Sinc3SqrtDD,
)


def test_pac_matches_independent_quadrature_and_preserves_arc_length():
    c0, c1, phi, length, xi = sp.symbols("c0 c1 phi length xi", real=True)
    transform = pac_transform(c0, c1, phi, length, xi)
    evaluate = sp.lambdify(
        (c0, c1, phi, length, xi), transform, [LAMBDA_MODULES, "numpy"]
    )
    nodes, weights = np.polynomial.legendre.leggauss(80)
    samples = 0.5 * (nodes + 1.0)
    for values in ((0.4, 0.3, 0.2, 0.7), (-0.2, -0.5, -0.7, 0.45)):
        c0_value, c1_value, phi_value, length_value = values
        actual = np.asarray(
            evaluate(c0_value, c1_value, phi_value, length_value, 1.0), dtype=float
        )
        alpha = c0_value * samples + 0.5 * c1_value * samples**2
        sine = 0.5 * np.dot(weights, np.sin(alpha))
        cosine = 0.5 * np.dot(weights, np.cos(alpha))
        expected = length_value * np.array([
            np.cos(phi_value) * sine, np.sin(phi_value) * sine, cosine
        ])
        np.testing.assert_allclose(actual[:3, 3], expected, rtol=2e-13, atol=2e-14)

    tangent = transform[:3, 3].diff(xi)
    assert sp.trigsimp(tangent.dot(tangent)) == length**2


def test_pac_straight_limit_and_first_two_derivatives_are_finite():
    c0, c1 = sp.symbols("c0 c1", real=True)
    expressions = [
        AffineCosMoment(0, c0, c1, 1), AffineSinMoment(0, c0, c1, 1)
    ]
    derivatives = expressions + [
        sp.diff(expression, variable, order)
        for expression in expressions
        for variable in (c0, c1)
        for order in (1, 2)
    ]
    values = sp.lambdify(
        (c0, c1), derivatives, [LAMBDA_MODULES, "numpy"]
    )(0.0, 0.0)
    assert np.isfinite(np.asarray(values, dtype=float)).all()
    np.testing.assert_allclose(values[:2], [1.0, 0.0], atol=1e-15)


def test_pcs_zero_curvature_is_finite_and_straight():
    bx, by, length = sp.symbols("bx by length", real=True)
    H = pcs_transform(
        sp.Matrix([-by / length, bx / length, 0]), sp.Matrix([0, 0, 1]), length
    )
    evaluate = sp.lambdify((bx, by, length), H, [LAMBDA_MODULES, "numpy"])
    actual = np.asarray(evaluate(0.0, 0.0, 0.5), dtype=float)
    expected = np.eye(4)
    expected[2, 3] = 0.5
    np.testing.assert_allclose(actual, expected, atol=1e-14)


def test_pcs_first_and_second_bending_derivatives_are_finite_at_zero():
    bx, by = sp.symbols("bx by", real=True)
    length = sp.Rational(1, 2)
    H = pcs_transform(
        sp.Matrix([-by / length, bx / length, 0]), sp.Matrix([0, 0, 1]), length
    )
    derivatives = list(H.diff(bx)) + list(H.diff(bx, 2)) + list(H.diff(bx, by))
    evaluate = sp.lambdify((bx, by), derivatives, [LAMBDA_MODULES, "numpy"])
    assert np.isfinite(np.asarray(evaluate(0.0, 0.0), dtype=float)).all()


def test_sinc3_analytic_continuation_and_derivatives_at_zero():
    z = sp.Symbol("z", real=True)
    values = sp.lambdify(
        z, [Sinc3Sqrt(z), Sinc3SqrtD(z), Sinc3SqrtDD(z)],
        [LAMBDA_MODULES, "numpy"],
    )(0.0)
    np.testing.assert_allclose(values, [1 / 6, -1 / 120, 1 / 2520], atol=1e-15)


def test_cosserat_reference_is_straight_and_derivatives_are_finite():
    kx, ky, kz, vx, vy, vz = sp.symbols("kx ky kz vx vy vz", real=True)
    H = pcs_transform(
        sp.Matrix([kx, ky, kz]), sp.Matrix([vx, vy, 1 + vz]), sp.Rational(1, 2)
    )
    evaluate = sp.lambdify(
        (kx, ky, kz, vx, vy, vz), H, [LAMBDA_MODULES, "numpy"]
    )
    actual = np.asarray(evaluate(0, 0, 0, 0, 0, 0), dtype=float)
    expected = np.eye(4)
    expected[2, 3] = 0.5
    np.testing.assert_allclose(actual, expected, atol=1e-14)
    derivatives = list(H.diff(kx)) + list(H.diff(kx, 2)) + list(H.diff(kx, ky))
    derivative_values = sp.lambdify(
        (kx, ky, kz, vx, vy, vz), derivatives, [LAMBDA_MODULES, "numpy"]
    )(0, 0, 0, 0, 0, 0)
    assert np.isfinite(np.asarray(derivative_values, dtype=float)).all()


def test_restricted_bend_stretch_is_one_pcs_strain_choice():
    bx, by, length = 0.17, -0.09, 0.53
    rest_length = 0.47
    kappa = sp.Matrix([-by / rest_length, bx / rest_length, 0.0])
    nu = sp.Matrix([0.0, 0.0, length / rest_length])
    actual = sp.lambdify(
        (), pcs_transform(kappa, nu, sp.Float(rest_length)),
        [LAMBDA_MODULES, "numpy"],
    )()
    theta = np.hypot(bx, by)
    sinc = np.sin(theta) / theta
    cosc = (1 - np.cos(theta)) / theta**2
    expected_position = [length * bx * cosc, length * by * cosc, length * sinc]
    np.testing.assert_allclose(actual[:3, 3], expected_position, rtol=1e-13, atol=1e-14)
