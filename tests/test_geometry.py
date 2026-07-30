import numpy as np
import sympy as sp

from softarm.geometry import pcc_transform
from softarm.special import LAMBDA_MODULES


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

