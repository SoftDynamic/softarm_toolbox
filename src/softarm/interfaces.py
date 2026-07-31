from __future__ import annotations

from typing import Protocol

import sympy as sp


class ActuatorMap(Protocol):
    """Optional physical actuator mapping to arm-only generalized force."""

    coordinates: sp.Matrix
    jacobian: sp.Matrix
    velocity_bias: sp.Matrix

    def generalized_force(self, tension: sp.Matrix) -> sp.Matrix:
        ...


class ConstraintMap(Protocol):
    """Acceleration-level constraint with a possibly non-ideal reaction map."""

    coordinates: sp.Matrix
    jacobian: sp.Matrix
    velocity_bias: sp.Matrix
    reaction_map: sp.Matrix
    stabilization_frequency: sp.Matrix
    stabilization_ratio: sp.Matrix
