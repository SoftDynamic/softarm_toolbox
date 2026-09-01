from __future__ import annotations

from collections.abc import Callable
from dataclasses import dataclass

import sympy as sp

from .config import ConstraintConfig
from .models import PlantModel, RuntimeParameter


@dataclass(frozen=True)
class ConstraintDefinition:
    """Target-independent point-jet constraint metadata and runtime parameters."""

    family: str
    channel_names: tuple[str, ...]
    channel_kinds: tuple[str, ...]
    parameters: tuple[RuntimeParameter, ...]
    tool_offset: sp.Matrix
    plane_point: sp.Matrix
    plane_normal: sp.Matrix
    friction: sp.Expr
    friction_velocity: sp.Expr
    stabilization_frequency: sp.Matrix
    stabilization_ratio: sp.Matrix

    @property
    def count(self) -> int:
        return len(self.channel_names)


@dataclass(frozen=True)
class ConstraintModel:
    """Acceleration-level constraints and their generalized reaction mapping."""

    family: str
    channel_names: tuple[str, ...]
    channel_kinds: tuple[str, ...]
    parameters: tuple[RuntimeParameter, ...]
    coordinates: sp.Matrix
    jacobian: sp.Matrix
    velocity_bias: sp.Matrix
    reaction_map: sp.Matrix
    stabilization_frequency: sp.Matrix
    stabilization_ratio: sp.Matrix

    @property
    def count(self) -> int:
        return len(self.channel_names)

    def generalized_force(self, reaction: sp.Matrix) -> sp.Matrix:
        if reaction.shape != (self.count, 1):
            raise ValueError(f"constraint reaction must have shape ({self.count}, 1)")
        return self.reaction_map * reaction

    def is_feasible(self, reaction: sp.Matrix) -> bool:
        if reaction.shape != (self.count, 1):
            raise ValueError(f"constraint reaction must have shape ({self.count}, 1)")
        return all(
            kind == "bilateral" or float(reaction[index]) >= 0.0
            for index, kind in enumerate(self.channel_kinds)
        )


ConstraintBuilder = Callable[[PlantModel, ConstraintConfig], ConstraintModel]


def _plane_point_contact_definition(config: ConstraintConfig) -> ConstraintDefinition:
    data = config.data
    parameters: list[RuntimeParameter] = []

    def scalar(name: str, default: float, positive: bool = False) -> sp.Symbol:
        symbol = sp.Symbol(f"constraint_{name}", real=True, positive=positive)
        parameters.append(RuntimeParameter(str(symbol), symbol, float(data.get(name, default))))
        return symbol

    def vector(name: str, default: tuple[float, float, float]) -> sp.Matrix:
        values = data.get(name, default)
        symbols = []
        for axis, value in zip("xyz", values, strict=True):
            symbol = sp.Symbol(f"constraint_{name}_{axis}", real=True)
            parameters.append(RuntimeParameter(str(symbol), symbol, float(value)))
            symbols.append(symbol)
        return sp.Matrix(symbols)

    tool_offset = vector("tool_offset", (0.0, 0.0, 0.0))
    plane_point = vector("plane_point", (0.0, 0.0, 0.0))
    normal_raw = vector("plane_normal", (0.0, 0.0, -1.0))
    normal_norm = sp.sqrt((normal_raw.T * normal_raw)[0])
    friction = scalar("friction", 0.0)
    friction_velocity = scalar("friction_velocity", 0.01, positive=True)
    frequency = scalar("stabilization_frequency", 20.0, positive=True)
    ratio = scalar("stabilization_ratio", 1.0)
    return ConstraintDefinition(
        family="plane_point_contact",
        channel_names=("normal_contact",),
        channel_kinds=("unilateral",),
        parameters=tuple(parameters),
        tool_offset=tool_offset,
        plane_point=plane_point,
        plane_normal=normal_raw / normal_norm,
        friction=friction,
        friction_velocity=friction_velocity,
        stabilization_frequency=sp.Matrix([frequency]),
        stabilization_ratio=sp.Matrix([ratio]),
    )


def derive_constraint_definition(
    config: ConstraintConfig | None,
) -> ConstraintDefinition | None:
    """Build constraint metadata without requiring global plant kinematics."""
    if config is None:
        return None
    if config.family != "plane_point_contact":
        raise ValueError(
            f"constraint family {config.family!r} has no recursive point-jet definition"
        )
    return _plane_point_contact_definition(config)


def point_constraint_expressions(
    definition: ConstraintDefinition,
    point_position: sp.Matrix,
    point_jacobian: sp.Matrix,
    point_velocity_bias: sp.Matrix,
    point_velocity: sp.Matrix,
) -> tuple[sp.Matrix, sp.Matrix, sp.Matrix, sp.Matrix]:
    """Evaluate the canonical plane-contact formula from a tool-point jet."""
    if definition.family != "plane_point_contact":
        raise ValueError(f"unsupported point constraint family: {definition.family}")
    normal = definition.plane_normal
    gap = sp.Matrix([(
        normal.T * (point_position - definition.plane_point)
    )[0]])
    jacobian = normal.T * point_jacobian
    velocity_bias = normal.T * point_velocity_bias
    tangential_velocity = (sp.eye(3) - normal * normal.T) * point_velocity
    regularized_speed = sp.sqrt(
        (tangential_velocity.T * tangential_velocity)[0]
        + definition.friction_velocity**2
    )
    contact_direction = (
        normal - definition.friction * tangential_velocity / regularized_speed
    )
    reaction_map = point_jacobian.T * contact_direction
    return gap, jacobian, velocity_bias, reaction_map


def _plane_point_contact_builder(
    plant: PlantModel,
    config: ConstraintConfig,
) -> ConstraintModel:
    definition = _plane_point_contact_definition(config)

    end_rotation = plant.end_transform[:3, :3]
    end_position = plant.end_transform[:3, 3]
    contact_position = end_position + end_rotation * definition.tool_offset
    point_jacobian = contact_position.jacobian(plant.q)
    point_velocity = point_jacobian * plant.dq
    point_velocity_bias = (
        point_velocity.jacobian(plant.q) * plant.dq
    )
    gap, jacobian, velocity_bias, reaction_map = point_constraint_expressions(
        definition,
        contact_position,
        point_jacobian,
        point_velocity_bias,
        point_velocity,
    )

    return ConstraintModel(
        family=definition.family,
        channel_names=definition.channel_names,
        channel_kinds=definition.channel_kinds,
        parameters=definition.parameters,
        coordinates=gap,
        jacobian=jacobian,
        velocity_bias=velocity_bias,
        reaction_map=reaction_map,
        stabilization_frequency=definition.stabilization_frequency,
        stabilization_ratio=definition.stabilization_ratio,
    )


_CONSTRAINT_BUILDERS: dict[str, ConstraintBuilder] = {
    "plane_point_contact": _plane_point_contact_builder,
}


def register_constraint(name: str, builder: ConstraintBuilder) -> None:
    if not name or name in _CONSTRAINT_BUILDERS:
        raise ValueError(f"constraint family {name!r} is already registered or invalid")
    _CONSTRAINT_BUILDERS[name] = builder


def derive_constraint(
    plant: PlantModel,
    config: ConstraintConfig | None = None,
) -> ConstraintModel | None:
    selected = config if config is not None else plant.config.constraint
    if selected is None:
        return None
    try:
        builder = _CONSTRAINT_BUILDERS[selected.family]
    except KeyError as error:
        raise ValueError(f"unregistered constraint family: {selected.family}") from error
    result = builder(plant, selected)
    m, nq = result.count, len(plant.q)
    expected = {
        "coordinates": ((m, 1), result.coordinates.shape),
        "jacobian": ((m, nq), result.jacobian.shape),
        "velocity bias": ((m, 1), result.velocity_bias.shape),
        "reaction map": ((nq, m), result.reaction_map.shape),
        "stabilization frequency": ((m, 1), result.stabilization_frequency.shape),
        "stabilization ratio": ((m, 1), result.stabilization_ratio.shape),
    }
    for label, (wanted, actual) in expected.items():
        if actual != wanted:
            raise ValueError(
                f"constraint builder returned invalid {label} shape {actual}; expected {wanted}"
            )
    if len(result.channel_kinds) != m or any(
        kind not in {"bilateral", "unilateral"} for kind in result.channel_kinds
    ):
        raise ValueError("constraint channel kinds must be bilateral or unilateral")
    parameter_names = [item.name for item in plant.parameters + result.parameters]
    if len(parameter_names) != len(set(parameter_names)):
        raise ValueError("constraint and plant parameter names must be unique")
    return result
