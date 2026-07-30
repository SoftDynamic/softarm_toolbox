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


LAMBDA_MODULES = {
    "SincSqrt": sinc_sqrt,
    "SincSqrtD": sinc_sqrt_d,
    "SincSqrtDD": sinc_sqrt_dd,
    "CoscSqrt": cosc_sqrt,
    "CoscSqrtD": cosc_sqrt_d,
    "CoscSqrtDD": cosc_sqrt_dd,
}

