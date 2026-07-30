from __future__ import annotations

import sympy as sp

from .special import CoscSqrt, SincSqrt


def homogeneous(rotation: sp.Matrix, position: sp.Matrix) -> sp.Matrix:
    return rotation.row_join(position).col_join(sp.Matrix([[0, 0, 0, 1]]))


def rotation_rpy(roll: sp.Expr, pitch: sp.Expr, yaw: sp.Expr) -> sp.Matrix:
    """ZYX roll-pitch-yaw rotation from a local frame to its parent frame."""
    cr, sr = sp.cos(roll), sp.sin(roll)
    cp, spitch = sp.cos(pitch), sp.sin(pitch)
    cy, sy = sp.cos(yaw), sp.sin(yaw)
    return sp.Matrix([
        [cy * cp, cy * spitch * sr - sy * cr, cy * spitch * cr + sy * sr],
        [sy * cp, sy * spitch * sr + cy * cr, sy * spitch * cr - cy * sr],
        [-spitch, cp * sr, cp * cr],
    ])


def transform_rpy(
    position: tuple[sp.Expr, sp.Expr, sp.Expr] | sp.Matrix,
    rpy: tuple[sp.Expr, sp.Expr, sp.Expr],
) -> sp.Matrix:
    return homogeneous(rotation_rpy(*rpy), sp.Matrix(position))


def pcc_transform(bx: sp.Expr, by: sp.Expr, length: sp.Expr, xi: sp.Expr = sp.S.One) -> sp.Matrix:
    x = xi * bx
    y = xi * by
    z = x * x + y * y
    a = SincSqrt(z)
    b = CoscSqrt(z)
    rotation = sp.Matrix([
        [1 - x * x * b, -x * y * b, x * a],
        [-x * y * b, 1 - y * y * b, y * a],
        [-x * a, -y * a, 1 - z * b],
    ])
    position = sp.Matrix([length * xi * x * b, length * xi * y * b, length * xi * a])
    return homogeneous(rotation, position)


def polynomial(coefficients: tuple[float, ...], xi: sp.Expr) -> sp.Expr:
    return sp.Add(*(sp.Rational(str(value)) * xi**index for index, value in enumerate(coefficients)))


def euler_ritz_transform(
    ax: sp.Expr,
    ay: sp.Expr,
    length: sp.Expr,
    xi: sp.Expr,
    psi_x: sp.Expr,
    psi_y: sp.Expr,
    dpsi_x: sp.Expr,
    dpsi_y: sp.Expr,
) -> sp.Matrix:
    slope_x = ax * dpsi_x / length
    slope_y = ay * dpsi_y / length
    # Consistent first-order Euler-Bernoulli cross-section rotation.
    rotation = sp.Matrix([
        [1, 0, slope_x],
        [0, 1, slope_y],
        [-slope_x, -slope_y, 1],
    ])
    position = sp.Matrix([ax * psi_x, ay * psi_y, length * xi])
    return homogeneous(rotation, position)


def vex(skew: sp.Matrix) -> sp.Matrix:
    return sp.Matrix([skew[2, 1], skew[0, 2], skew[1, 0]])


def angular_jacobian(rotation: sp.Matrix, q: sp.Matrix) -> sp.Matrix:
    columns = []
    for coordinate in q:
        rate = sp.diff(rotation, coordinate) * rotation.T
        skew = (rate - rate.T) / 2
        columns.append(vex(skew))
    return sp.Matrix.hstack(*columns)
