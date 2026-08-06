from __future__ import annotations

from abc import ABC, abstractmethod
from dataclasses import dataclass

import sympy as sp

from .config import ModelConfig


@dataclass(frozen=True)
class RuntimeParameter:
    name: str
    symbol: sp.Symbol
    default: float


class PlantModel(ABC):
    """Common build-time interface consumed by actuators and code generators."""

    config: ModelConfig
    q: sp.Matrix
    dq: sp.Matrix
    base_q: sp.Matrix
    base_dq: sp.Matrix
    arm_q: sp.Matrix
    arm_dq: sp.Matrix
    parameters: tuple[RuntimeParameter, ...]
    kinematics: sp.Matrix
    end_jacobian: sp.Matrix
    base_transform: sp.Matrix
    end_transform: sp.Matrix
    base_jacobian: sp.Matrix
    vehicle_wrench_map: sp.Matrix
    arm_force_map: sp.Matrix

    @property
    def p(self) -> sp.Matrix:
        return sp.Matrix([item.symbol for item in self.parameters])

    @property
    def coordinate_names(self) -> list[str]:
        return [str(item) for item in self.q]

    @property
    def base_coordinate_names(self) -> list[str]:
        return [str(item) for item in self.base_q]

    @property
    def arm_coordinate_names(self) -> list[str]:
        return [str(item) for item in self.arm_q]

    @property
    @abstractmethod
    def bias(self) -> sp.Matrix:
        """Return the generalized velocity, conservative, and damping force."""


@dataclass
class SymbolicPlant(PlantModel):
    config: ModelConfig
    q: sp.Matrix
    dq: sp.Matrix
    base_q: sp.Matrix
    base_dq: sp.Matrix
    arm_q: sp.Matrix
    arm_dq: sp.Matrix
    parameters: tuple[RuntimeParameter, ...]
    mass: sp.Matrix
    potential: sp.Expr
    damping: sp.Matrix
    kinematics: sp.Matrix
    end_jacobian: sp.Matrix
    base_transform: sp.Matrix
    end_transform: sp.Matrix
    base_jacobian: sp.Matrix
    vehicle_wrench_map: sp.Matrix
    arm_force_map: sp.Matrix
    _material_coordinate: sp.Symbol | None = None
    _material_kinematics: sp.Matrix | None = None
    _bias: sp.Matrix | None = None

    @property
    def bias(self) -> sp.Matrix:
        if self._bias is None:
            # Local import keeps the plant representation independent of a
            # particular dynamics assembly implementation.
            from .dynamics import SymbolicLagrangeAssembler

            self._bias = SymbolicLagrangeAssembler().assemble_bias(self)
        return self._bias
