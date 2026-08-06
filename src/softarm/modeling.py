from __future__ import annotations

from abc import ABC, abstractmethod
from collections.abc import Callable
from dataclasses import dataclass

import sympy as sp

from .config import ModelConfig
from .models import PlantModel, RuntimeParameter


class SectionKinematics(ABC):
    """Local section pose evaluated at a normalized material coordinate."""

    @abstractmethod
    def transform(self, section: int, xi: sp.Expr) -> sp.Matrix:
        """Return the local homogeneous transform for one section."""


@dataclass(frozen=True)
class FunctionalSectionKinematics(SectionKinematics):
    """Adapter for the existing symbolic local-transform definitions."""

    transform_function: Callable[[int, sp.Expr], sp.Matrix]

    def transform(self, section: int, xi: sp.Expr) -> sp.Matrix:
        return self.transform_function(section, xi)


@dataclass(frozen=True)
class ModelCoordinates:
    arm_q: sp.Matrix
    arm_dq: sp.Matrix


@dataclass(frozen=True)
class SectionProperties:
    lengths: tuple[sp.Symbol, ...]
    masses: tuple[sp.Symbol, ...]
    inertias: tuple[sp.Matrix, ...]


@dataclass(frozen=True)
class ModelDefinition:
    """Model-specific data required by a dynamics assembly strategy."""

    config: ModelConfig
    coordinates: ModelCoordinates
    parameters: tuple[RuntimeParameter, ...]
    kinematics: SectionKinematics
    sections: SectionProperties
    damping: tuple[sp.Expr, ...] | sp.Matrix
    elastic: sp.Expr
    distributed: bool
    linear_kinematics: bool = False


class DynamicsAssembler(ABC):
    """Strategy interface that converts a model definition into a plant."""

    @abstractmethod
    def assemble(self, definition: ModelDefinition) -> PlantModel:
        """Assemble one plant representation."""


class ModelDefinitionBuilder(ABC):
    """Template method for validating and assembling a registered model."""

    def build(
        self,
        config: ModelConfig,
        assembler: DynamicsAssembler,
    ) -> PlantModel:
        self.validate(config)
        return assembler.assemble(self.define(config))

    def validate(self, config: ModelConfig) -> None:
        """Validate model-specific assumptions before symbolic construction."""
        return None

    @abstractmethod
    def define(self, config: ModelConfig) -> ModelDefinition:
        """Create the assembly-neutral definition for a model combination."""
