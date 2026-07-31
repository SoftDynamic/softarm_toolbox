import numpy as np
import sympy as sp

from softarm.geometry import cosserat_pcs_transform, pcc_transform
from softarm.special import (
    LAMBDA_MODULES, Sinc3Sqrt, Sinc3SqrtD, Sinc3SqrtDD,
)


def test_pcc_zero_curvature_is_finite_and_straight():
    bx, by, length = sp.symbols("bx by length", real=True)
    H = pcc_transform(bx, by, length)
    evaluate = sp.lambdify((bx, by, length), H, [LAMBDA_MODULES, "numpy"])
    actual = np.asarray(evaluate(0.0, 0.0, 0.5), dtype=float)
    expected = np.eye(4)
    expected[2, 3] = 0.5
    np.testing.assert_allclose(actual, expected, atol=1e-14)


def test_pcc_first_and_second_derivatives_are_finite_at_zero():
    bx, by = sp.symbols("bx by", real=True)
    H = pcc_transform(bx, by, sp.Rational(1, 2))
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
    H = cosserat_pcs_transform(
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


def test_cosserat_reduces_exactly_to_pcc_kinematics():
    bx, by, length = 0.17, -0.09, 0.53
    rest_length = 0.47
    pcc = sp.lambdify(
        (), pcc_transform(sp.Float(bx), sp.Float(by), sp.Float(length)),
        [LAMBDA_MODULES, "numpy"],
    )()
    kappa = sp.Matrix([-by / rest_length, bx / rest_length, 0.0])
    nu = sp.Matrix([0.0, 0.0, length / rest_length])
    cosserat = sp.lambdify(
        (), cosserat_pcs_transform(kappa, nu, sp.Float(rest_length)),
        [LAMBDA_MODULES, "numpy"],
    )()
    np.testing.assert_allclose(cosserat, pcc, rtol=1e-13, atol=1e-14)
