from __future__ import annotations

import math

import sympy as sp


class SincSqrt(sp.Function):
    nargs = 1

    def fdiff(self, argindex=1):
        if argindex != 1:
            raise sp.ArgumentIndexError(self, argindex)
        return SincSqrtD(self.args[0])


class SincSqrtD(sp.Function):
    nargs = 1

    def fdiff(self, argindex=1):
        if argindex != 1:
            raise sp.ArgumentIndexError(self, argindex)
        return SincSqrtDD(self.args[0])


class SincSqrtDD(sp.Function):
    nargs = 1


class CoscSqrt(sp.Function):
    nargs = 1

    def fdiff(self, argindex=1):
        if argindex != 1:
            raise sp.ArgumentIndexError(self, argindex)
        return CoscSqrtD(self.args[0])


class CoscSqrtD(sp.Function):
    nargs = 1

    def fdiff(self, argindex=1):
        if argindex != 1:
            raise sp.ArgumentIndexError(self, argindex)
        return CoscSqrtDD(self.args[0])


class CoscSqrtDD(sp.Function):
    nargs = 1


class Sinc3Sqrt(sp.Function):
    """Analytic continuation of (1 - sinc(sqrt(z))) / z."""

    nargs = 1

    def fdiff(self, argindex=1):
        if argindex != 1:
            raise sp.ArgumentIndexError(self, argindex)
        return Sinc3SqrtD(self.args[0])


class Sinc3SqrtD(sp.Function):
    nargs = 1

    def fdiff(self, argindex=1):
        if argindex != 1:
            raise sp.ArgumentIndexError(self, argindex)
        return Sinc3SqrtDD(self.args[0])


class Sinc3SqrtDD(sp.Function):
    nargs = 1


# This registry is the sole source used to teach generic CAS executors the
# symbolic derivative heads. It contains no backend-specific formula.
SPECIAL_DERIVATIVE_HEADS = (
    ("SincSqrt", "SincSqrtD"),
    ("SincSqrtD", "SincSqrtDD"),
    ("CoscSqrt", "CoscSqrtD"),
    ("CoscSqrtD", "CoscSqrtDD"),
    ("Sinc3Sqrt", "Sinc3SqrtD"),
    ("Sinc3SqrtD", "Sinc3SqrtDD"),
)


def _series_sinc(z: float) -> float:
    return 1.0 - z / 6.0 + z * z / 120.0 - z**3 / 5040.0 + z**4 / 362880.0


def _series_sinc_d(z: float) -> float:
    return -1.0 / 6.0 + z / 60.0 - z * z / 1680.0 + z**3 / 90720.0


def _series_sinc_dd(z: float) -> float:
    return 1.0 / 60.0 - z / 840.0 + z * z / 30240.0


def _series_cosc(z: float) -> float:
    return 0.5 - z / 24.0 + z * z / 720.0 - z**3 / 40320.0 + z**4 / 3628800.0


def _series_cosc_d(z: float) -> float:
    return -1.0 / 24.0 + z / 360.0 - z * z / 13440.0 + z**3 / 907200.0


def _series_cosc_dd(z: float) -> float:
    return 1.0 / 360.0 - z / 6720.0 + z * z / 302400.0


def _series_sinc3(z: float) -> float:
    return 1.0 / 6.0 - z / 120.0 + z * z / 5040.0 - z**3 / 362880.0 + z**4 / 39916800.0


def _series_sinc3_d(z: float) -> float:
    return -1.0 / 120.0 + z / 2520.0 - z * z / 120960.0 + z**3 / 9979200.0


def _series_sinc3_dd(z: float) -> float:
    return 1.0 / 2520.0 - z / 60480.0 + z * z / 3326400.0


def sinc_sqrt(z: float) -> float:
    if abs(z) < 1e-8:
        return _series_sinc(z)
    root = math.sqrt(z)
    return math.sin(root) / root


def sinc_sqrt_d(z: float) -> float:
    if abs(z) < 1e-8:
        return _series_sinc_d(z)
    root = math.sqrt(z)
    return (root * math.cos(root) - math.sin(root)) / (2.0 * root**3)


def sinc_sqrt_dd(z: float) -> float:
    if abs(z) < 1e-8:
        return _series_sinc_dd(z)
    root = math.sqrt(z)
    return ((3.0 - z) * math.sin(root) - 3.0 * root * math.cos(root)) / (4.0 * root**5)


def cosc_sqrt(z: float) -> float:
    if abs(z) < 1e-8:
        return _series_cosc(z)
    root = math.sqrt(z)
    return (1.0 - math.cos(root)) / z


def cosc_sqrt_d(z: float) -> float:
    if abs(z) < 1e-8:
        return _series_cosc_d(z)
    root = math.sqrt(z)
    return (root * math.sin(root) - 2.0 * (1.0 - math.cos(root))) / (2.0 * z**2)


def cosc_sqrt_dd(z: float) -> float:
    if abs(z) < 1e-8:
        return _series_cosc_dd(z)
    root = math.sqrt(z)
    return (z * math.cos(root) - 5.0 * root * math.sin(root) + 8.0 - 8.0 * math.cos(root)) / (4.0 * z**3)


def sinc3_sqrt(z: float) -> float:
    if abs(z) < 1e-8:
        return _series_sinc3(z)
    return (1.0 - sinc_sqrt(z)) / z


def sinc3_sqrt_d(z: float) -> float:
    if abs(z) < 1e-8:
        return _series_sinc3_d(z)
    sinc = sinc_sqrt(z)
    return (sinc - 1.0 - z * sinc_sqrt_d(z)) / z**2


def sinc3_sqrt_dd(z: float) -> float:
    if abs(z) < 1e-8:
        return _series_sinc3_dd(z)
    sinc = sinc_sqrt(z)
    return (
        2.0 - 2.0 * sinc + 2.0 * z * sinc_sqrt_d(z) - z**2 * sinc_sqrt_dd(z)
    ) / z**3


LAMBDA_MODULES = {
    "SincSqrt": sinc_sqrt,
    "SincSqrtD": sinc_sqrt_d,
    "SincSqrtDD": sinc_sqrt_dd,
    "CoscSqrt": cosc_sqrt,
    "CoscSqrtD": cosc_sqrt_d,
    "CoscSqrtDD": cosc_sqrt_dd,
    "Sinc3Sqrt": sinc3_sqrt,
    "Sinc3SqrtD": sinc3_sqrt_d,
    "Sinc3SqrtDD": sinc3_sqrt_dd,
}
