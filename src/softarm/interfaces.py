from __future__ import annotations

from typing import Protocol
import sympy as sp


class ActuatorMap(Protocol):
    """Optional physical actuator mapping; the base plant consumes tau directly."""

    coordinates: sp.Matrix
    jacobian: sp.Matrix
    velocity_bias: sp.Matrix

    def generalized_force(self, tension: sp.Matrix) -> sp.Matrix:
        ...


class ConstraintMap(Protocol):
    """Optional acceleration-level constraint supplied by a future extension."""

    jacobian: sp.Matrix
    velocity_bias: sp.Matrix
