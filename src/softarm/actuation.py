from __future__ import annotations

from collections.abc import Callable
from dataclasses import dataclass

import sympy as sp

from .config import ActuationConfig
from .derive import RuntimeParameter, SymbolicPlant
from .geometry import polynomial


@dataclass(frozen=True)
class ActuationModel:
    """SymPy actuator coordinates and their virtual-work mapping."""

    family: str
    acceleration: str
    channel_names: tuple[str, ...]
    channel_kinds: tuple[str, ...]
    parameters: tuple[RuntimeParameter, ...]
    coordinates: sp.Matrix
    jacobian: sp.Matrix
    velocity_bias: sp.Matrix

    @property
    def count(self) -> int:
        return len(self.channel_names)

    def generalized_force(self, tension: sp.Matrix) -> sp.Matrix:
        if tension.shape != (self.count, 1):
            raise ValueError(f"actuator tension must have shape ({self.count}, 1)")
        return -self.jacobian.T * tension

    def is_feasible(self, tension: sp.Matrix) -> bool:
        if tension.shape != (self.count, 1):
            raise ValueError(f"actuator tension must have shape ({self.count}, 1)")
        return all(kind == "signed" or float(tension[index]) >= 0.0
                   for index, kind in enumerate(self.channel_kinds))


ActuatorBuilder = Callable[[SymbolicPlant, ActuationConfig], ActuationModel]


def _parameter_by_name(plant: SymbolicPlant) -> dict[str, sp.Symbol]:
    return {item.name: item.symbol for item in plant.parameters}


def _tendon_builder(
    plant: SymbolicPlant,
    config: ActuationConfig,
) -> ActuationModel:
    if not config.channels:
        raise ValueError("tendon actuation requires at least one channel")
    parameters: list[RuntimeParameter] = []
    coordinates: list[sp.Expr] = []
    plant_parameters = _parameter_by_name(plant)

    xi = sp.Symbol("xi", real=True)
    combination = (plant.config.rod, plant.config.parameterization)
    if plant.config.parameterization == "ritz":
        psi_x = polynomial(plant.config.ritz_x or (), xi)
        psi_y = polynomial(plant.config.ritz_y or (), xi)
        slope_x = sp.diff(psi_x, xi).subs(xi, sp.S.One)
        slope_y = sp.diff(psi_y, xi).subs(xi, sp.S.One)

    for channel in config.channels:
        coordinate = sp.S.Zero
        for span in channel.spans:
            section = span.section - 1
            parameter_name = f"act_{channel.name}_s{span.section}_radius"
            radius = sp.Symbol(parameter_name, positive=True, real=True)
            parameters.append(RuntimeParameter(parameter_name, radius, span.radius))
            cosine = sp.cos(sp.Float(str(span.angle)))
            sine = sp.sin(sp.Float(str(span.angle)))
            if combination == ("extensible_euler_bernoulli", "pcs"):
                offset = 3 * section
                bending = plant.arm_q[offset] * cosine + plant.arm_q[offset + 1] * sine
                if channel.kind == "unilateral":
                    coordinate += plant.arm_q[offset + 2]
                coordinate -= radius * bending
            elif combination == ("euler_bernoulli", "pcs"):
                offset = 2 * section
                length = plant_parameters[f"s{span.section}_length"]
                bending = plant.arm_q[offset] * cosine + plant.arm_q[offset + 1] * sine
                if channel.kind == "unilateral":
                    coordinate += length
                coordinate -= radius * bending
            elif combination == ("euler_bernoulli", "ritz"):
                offset = 2 * section
                length = plant_parameters[f"s{span.section}_length"]
                bending = (
                    plant.arm_q[offset] * cosine * slope_x
                    + plant.arm_q[offset + 1] * sine * slope_y
                ) / length
                if channel.kind == "unilateral":
                    coordinate += length
                coordinate -= radius * bending
            elif combination == ("extensible_euler_bernoulli", "ritz"):
                offset = 3 * section
                length = plant_parameters[f"s{span.section}_rest_length"]
                bending = (
                    plant.arm_q[offset] * cosine * slope_x
                    + plant.arm_q[offset + 1] * sine * slope_y
                ) / length
                if channel.kind == "unilateral":
                    coordinate += length + plant.arm_q[offset + 2]
                coordinate -= radius * bending
            elif combination == ("euler_bernoulli", "pac"):
                offset = 3 * section
                length = plant_parameters[f"s{span.section}_length"]
                total_bend = plant.arm_q[offset] + plant.arm_q[offset + 1] / 2
                bending = total_bend * sp.cos(
                    sp.Float(str(span.angle)) - plant.arm_q[offset + 2]
                )
                if channel.kind == "unilateral":
                    coordinate += length
                coordinate -= radius * bending
            elif combination == ("extensible_euler_bernoulli", "pac"):
                offset = 4 * section
                total_bend = plant.arm_q[offset] + plant.arm_q[offset + 1] / 2
                bending = total_bend * sp.cos(
                    sp.Float(str(span.angle)) - plant.arm_q[offset + 2]
                )
                if channel.kind == "unilateral":
                    coordinate += plant.arm_q[offset + 3]
                coordinate -= radius * bending
            elif combination == ("cosserat", "pcs"):
                offset = 6 * section
                length = plant_parameters[f"s{span.section}_length"]
                kappa = sp.Matrix([
                    plant_parameters[f"s{span.section}_kappa0_x"] + plant.arm_q[offset],
                    plant_parameters[f"s{span.section}_kappa0_y"] + plant.arm_q[offset + 1],
                    plant_parameters[f"s{span.section}_kappa0_z"] + plant.arm_q[offset + 2],
                ])
                nu = sp.Matrix([
                    plant_parameters[f"s{span.section}_nu0_x"] + plant.arm_q[offset + 3],
                    plant_parameters[f"s{span.section}_nu0_y"] + plant.arm_q[offset + 4],
                    plant_parameters[f"s{span.section}_nu0_z"] + plant.arm_q[offset + 5],
                ])
                routing_offset = sp.Matrix([radius * cosine, radius * sine, 0])
                plus_tangent = nu + kappa.cross(routing_offset)
                plus_length = length * sp.sqrt(plus_tangent.dot(plus_tangent))
                if channel.kind == "unilateral":
                    coordinate += plus_length
                else:
                    minus_tangent = nu - kappa.cross(routing_offset)
                    minus_length = length * sp.sqrt(minus_tangent.dot(minus_tangent))
                    coordinate += (plus_length - minus_length) / 2
            else:
                raise ValueError(
                    "built-in tendon routing does not support model combination "
                    f"{combination!r}"
                )
        coordinates.append(coordinate)

    coordinate_matrix = sp.Matrix(coordinates)
    jacobian = coordinate_matrix.jacobian(plant.arm_q)
    velocity_bias = (
        (jacobian * plant.arm_dq).jacobian(plant.arm_q)
        * plant.arm_dq
    )
    return ActuationModel(
        "tendon",
        config.acceleration,
        tuple(channel.name for channel in config.channels),
        tuple(channel.kind for channel in config.channels),
        tuple(parameters),
        coordinate_matrix,
        jacobian,
        velocity_bias,
    )


_ACTUATOR_BUILDERS: dict[str, ActuatorBuilder] = {"tendon": _tendon_builder}


def register_actuator(name: str, builder: ActuatorBuilder) -> None:
    """Register a custom SymPy actuator builder."""
    if not name or name in _ACTUATOR_BUILDERS:
        raise ValueError(f"actuator family {name!r} is already registered or invalid")
    _ACTUATOR_BUILDERS[name] = builder


def _validate_strict_rank(plant: SymbolicPlant, actuation: ActuationModel) -> None:
    if actuation.acceleration != "strict":
        return
    if actuation.count > len(plant.arm_q):
        raise ValueError(
            f"strict actuator acceleration has {actuation.count} constraints but only {len(plant.arm_q)} arm coordinates"
        )
    defaults = {
        item.symbol: item.default for item in plant.parameters + actuation.parameters
    }
    reference = {coordinate: 0.0 for coordinate in plant.arm_q}
    if (plant.config.rod, plant.config.parameterization) in {
        ("extensible_euler_bernoulli", "pcs"),
        ("extensible_euler_bernoulli", "pac"),
    }:
        stride = 3 if plant.config.parameterization == "pcs" else 4
        length_offset = 2 if plant.config.parameterization == "pcs" else 3
        for section in range(plant.config.segments):
            rest = next(
                item.default for item in plant.parameters
                if item.name == f"s{section + 1}_rest_length"
            )
            reference[plant.arm_q[stride * section + length_offset]] = rest
    nominal = actuation.jacobian.subs(defaults).subs(reference).evalf()

    def nearly_zero(value: sp.Expr) -> bool:
        return abs(complex(value)) <= 1e-10

    rank = nominal.rank(iszerofunc=nearly_zero)
    if rank != actuation.count:
        raise ValueError("strict actuator acceleration requires a full-row-rank nominal Jacobian")


def derive_actuation(
    plant: SymbolicPlant,
    config: ActuationConfig | None = None,
) -> ActuationModel | None:
    """Derive an optional actuator model without changing the base plant."""
    selected = config if config is not None else plant.config.actuation
    if selected is None:
        return None
    try:
        builder = _ACTUATOR_BUILDERS[selected.family]
    except KeyError as error:
        raise ValueError(f"unregistered actuator family: {selected.family}") from error
    result = builder(plant, selected)
    if result.coordinates.shape != (result.count, 1):
        raise ValueError("actuator builder returned inconsistent coordinate dimensions")
    if result.jacobian.shape != (result.count, len(plant.arm_q)):
        raise ValueError("actuator builder returned inconsistent Jacobian dimensions")
    if result.velocity_bias.shape != (result.count, 1):
        raise ValueError("actuator builder returned inconsistent velocity-bias dimensions")
    parameter_names = [item.name for item in plant.parameters + result.parameters]
    if len(parameter_names) != len(set(parameter_names)):
        raise ValueError("actuator and plant parameter names must be unique")
    _validate_strict_rank(plant, result)
    return result
