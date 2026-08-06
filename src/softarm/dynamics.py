from __future__ import annotations

from collections.abc import Sequence
from typing import Protocol

import sympy as sp

from .config import ModelConfig
from .geometry import angular_jacobian, transform_rpy
from .integration import integrate_unit, unit_gauss_rule
from .modeling import DynamicsAssembler, ModelDefinition
from .models import RuntimeParameter, SymbolicPlant


class BatchDifferentiator(Protocol):
    """Narrow symbolic capability required by dynamics assembly."""

    def differentiate(
        self,
        expressions: Sequence[sp.Expr],
        variables: Sequence[sp.Symbol],
    ) -> list[sp.Expr]: ...


class SympyBatchDifferentiator:
    """Reference implementation using the canonical SymPy expressions."""

    def differentiate(
        self,
        expressions: Sequence[sp.Expr],
        variables: Sequence[sp.Symbol],
    ) -> list[sp.Expr]:
        return [
            sp.diff(expression, variable)
            for expression in expressions
            for variable in variables
        ]


def _base_coordinates(config: ModelConfig) -> tuple[sp.Matrix, sp.Matrix]:
    if config.base.mode == "fixed":
        return sp.zeros(0, 1), sp.zeros(0, 1)
    q = sp.Matrix(
        sp.symbols(
            "base_x base_y base_z base_roll base_pitch base_yaw",
            real=True,
        )
    )
    dq = sp.Matrix(
        sp.symbols(
            "dbase_x dbase_y dbase_z dbase_roll dbase_pitch dbase_yaw",
            real=True,
        )
    )
    return q, dq


def _linearize_matrix(matrix: sp.Matrix, q: sp.Matrix) -> sp.Matrix:
    zero = {coordinate: 0 for coordinate in q}
    return matrix.applyfunc(
        lambda value: value.subs(zero)
        + sum(sp.diff(value, coordinate).subs(zero) * coordinate for coordinate in q)
    )


def _mixed_linear_angular_jacobian(
    rotation: sp.Matrix,
    q: sp.Matrix,
    linear_coordinates: sp.Matrix,
) -> sp.Matrix:
    """Spatial angular Jacobian exact in base attitude and first order in arm shape."""
    zero = {coordinate: 0 for coordinate in linear_coordinates}
    reference_rotation = rotation.subs(zero)
    columns = []
    for coordinate in q:
        rate = rotation.diff(coordinate).subs(zero) * reference_rotation.T
        skew = (rate - rate.T) / 2
        columns.append(sp.Matrix([skew[2, 1], skew[0, 2], skew[1, 0]]))
    return sp.Matrix.hstack(*columns)


def _point_terms(
    transform: sp.Matrix,
    q: sp.Matrix,
    mass: sp.Expr,
    inertia: sp.Matrix,
    linear_rotation: bool = False,
    linear_coordinates: sp.Matrix | None = None,
) -> tuple[sp.Matrix, sp.Expr]:
    position = transform[:3, 3]
    rotation = transform[:3, :3]
    linear_jacobian = position.jacobian(q)
    if linear_rotation:
        selected_coordinates = linear_coordinates if linear_coordinates is not None else q
        angular = _mixed_linear_angular_jacobian(rotation, q, selected_coordinates)
        rotation_for_inertia = rotation.subs(
            {coordinate: 0 for coordinate in selected_coordinates}
        )
    else:
        angular = angular_jacobian(rotation, q)
        rotation_for_inertia = rotation
    mass_term = (
        mass * (linear_jacobian.T * linear_jacobian)
        + angular.T
        * (rotation_for_inertia * inertia * rotation_for_inertia.T)
        * angular
    )
    return mass_term, position[2]


class SymbolicLagrangeAssembler(DynamicsAssembler):
    """Assemble the current global symbolic Euler-Lagrange plant."""

    def __init__(self, differentiator: BatchDifferentiator | None = None):
        self.differentiator = differentiator or SympyBatchDifferentiator()

    def assemble(self, definition: ModelDefinition) -> SymbolicPlant:
        config = definition.config
        arm_q = definition.coordinates.arm_q
        arm_dq = definition.coordinates.arm_dq
        parameters = list(definition.parameters)

        def global_parameter(name: str, default: float) -> sp.Symbol:
            value = float(config.parameters.get(name, default))
            symbol = sp.Symbol(name, real=True)
            parameters.append(RuntimeParameter(name, symbol, value))
            return symbol

        gravity = global_parameter("gravity", 9.81)
        tip_mass = global_parameter("tip_mass", 0.1)
        tip_ixx = global_parameter("tip_Ixx", 0.01)
        tip_iyy = global_parameter("tip_Iyy", 0.01)
        tip_izz = global_parameter("tip_Izz", 0.01)

        base_q, base_dq = _base_coordinates(config)
        q = base_q.col_join(arm_q)
        dq = base_dq.col_join(arm_dq)
        nq = len(q)
        mass_matrix = sp.zeros(nq)
        potential = definition.elastic
        if config.base.mode == "floating_rpy":
            vehicle_mass = global_parameter("vehicle_mass", 1.5)
            vehicle_ixx = global_parameter("vehicle_Ixx", 0.03)
            vehicle_iyy = global_parameter("vehicle_Iyy", 0.03)
            vehicle_izz = global_parameter("vehicle_Izz", 0.05)
            base_transform = transform_rpy(base_q[:3, 0], tuple(base_q[3:6, 0]))
            vehicle_inertia = sp.diag(vehicle_ixx, vehicle_iyy, vehicle_izz)
            vehicle_term, vehicle_height = _point_terms(
                base_transform, q, vehicle_mass, vehicle_inertia
            )
            mass_matrix += vehicle_term
            potential -= vehicle_mass * gravity * vehicle_height
        else:
            base_transform = sp.eye(4)

        mount = transform_rpy(
            tuple(sp.Float(str(value)) for value in config.base.mount_xyz),
            tuple(sp.Float(str(value)) for value in config.base.mount_rpy),
        )
        base = base_transform * mount
        transforms: list[sp.Matrix] = []
        material_transforms: list[sp.Matrix] = []
        xi = sp.Symbol("xi", real=True, nonnegative=True)

        for section in range(config.segments):
            local_end = definition.kinematics.transform(section, sp.S.One)
            end = base * local_end
            if definition.linear_kinematics:
                end = _linearize_matrix(end, arm_q)
            transforms.append(end)

            material_transform = base * definition.kinematics.transform(section, xi)
            if definition.linear_kinematics:
                material_transform = _linearize_matrix(material_transform, arm_q)
            material_transforms.append(material_transform)

            mass = definition.sections.masses[section]
            inertia = definition.sections.inertias[section]
            if definition.distributed:
                if config.integration.method == "gauss":
                    for node, weight in unit_gauss_rule(config.integration):
                        material = base * definition.kinematics.transform(section, node)
                        if definition.linear_kinematics:
                            material = _linearize_matrix(material, arm_q)
                        point_mass, height = _point_terms(
                            material,
                            q,
                            mass,
                            inertia,
                            definition.linear_kinematics,
                            arm_q if definition.linear_kinematics else None,
                        )
                        mass_matrix += weight * point_mass
                        potential -= mass * gravity * weight * height
                else:
                    material = base * definition.kinematics.transform(section, xi)
                    if definition.linear_kinematics:
                        material = _linearize_matrix(material, arm_q)
                    point_mass, height = _point_terms(
                        material,
                        q,
                        mass,
                        inertia,
                        definition.linear_kinematics,
                        arm_q if definition.linear_kinematics else None,
                    )
                    mass_matrix += integrate_unit(
                        point_mass,
                        xi,
                        config.integration,
                        f"section {section + 1} inertia",
                    )
                    potential -= mass * gravity * integrate_unit(
                        height,
                        xi,
                        config.integration,
                        f"section {section + 1} gravity",
                    )
            else:
                midpoint = base * definition.kinematics.transform(
                    section, sp.Rational(1, 2)
                )
                if definition.linear_kinematics:
                    midpoint = _linearize_matrix(midpoint, arm_q)
                point_mass, height = _point_terms(
                    midpoint,
                    q,
                    mass,
                    inertia,
                    definition.linear_kinematics,
                    arm_q if definition.linear_kinematics else None,
                )
                mass_matrix += point_mass
                potential -= mass * gravity * height
            base = end

        end_position = base[:3, 3]
        end_rotation = base[:3, :3]
        linear_end_jacobian = end_position.jacobian(q)
        angular_end_jacobian = (
            _mixed_linear_angular_jacobian(end_rotation, q, arm_q)
            if definition.linear_kinematics
            else angular_jacobian(end_rotation, q)
        )
        tip_inertia = sp.diag(tip_ixx, tip_iyy, tip_izz)
        tip_rotation = (
            end_rotation.subs({coordinate: 0 for coordinate in arm_q})
            if definition.linear_kinematics
            else end_rotation
        )
        mass_matrix += (
            tip_mass * (linear_end_jacobian.T * linear_end_jacobian)
            + angular_end_jacobian.T
            * (tip_rotation * tip_inertia * tip_rotation.T)
            * angular_end_jacobian
        )
        potential -= tip_mass * gravity * end_position[2]

        kinematics = sp.Matrix.hstack(*transforms)
        material_kinematics = sp.Matrix.hstack(*material_transforms)
        end_jacobian = linear_end_jacobian.col_join(angular_end_jacobian)
        base_position = base_transform[:3, 3]
        base_rotation = base_transform[:3, :3]
        if len(base_q):
            base_linear_jacobian = base_position.jacobian(q)
            base_angular_jacobian = angular_jacobian(base_rotation, q)
            base_jacobian = base_linear_jacobian.col_join(base_angular_jacobian)
            vehicle_wrench_map = base_jacobian.T * sp.diag(base_rotation, base_rotation)
        else:
            base_jacobian = sp.zeros(6, nq)
            vehicle_wrench_map = sp.zeros(nq, 6)
        arm_force_map = sp.zeros(nq, len(arm_q))
        if len(arm_q):
            arm_force_map[len(base_q) :, :] = sp.eye(len(arm_q))
        full_damping = sp.zeros(nq)
        arm_damping = (
            sp.diag(*definition.damping)
            if isinstance(definition.damping, tuple)
            else definition.damping
        )
        if arm_damping.shape != (len(arm_q), len(arm_q)):
            raise ValueError("arm damping matrix has inconsistent dimensions")
        if len(arm_q):
            full_damping[len(base_q) :, len(base_q) :] = arm_damping
        return SymbolicPlant(
            config,
            q,
            dq,
            base_q,
            base_dq,
            arm_q,
            arm_dq,
            tuple(parameters),
            mass_matrix,
            potential,
            full_damping,
            kinematics,
            end_jacobian,
            base_transform,
            base,
            base_jacobian,
            vehicle_wrench_map,
            arm_force_map,
            _material_coordinate=xi,
            _material_kinematics=material_kinematics,
        )

    def assemble_bias(self, plant: SymbolicPlant) -> sp.Matrix:
        return assemble_bias(plant, self.differentiator)


def assemble_bias(
    plant: SymbolicPlant,
    differentiator: BatchDifferentiator,
) -> sp.Matrix:
    """Assemble Coriolis, conservative, and damping terms from one formula."""
    nq = len(plant.q)
    upper = [
        (row, column)
        for row in range(nq)
        for column in range(row, nq)
    ]
    derivatives = differentiator.differentiate(
        [plant.mass[row, column] for row, column in upper] + [plant.potential],
        list(plant.q),
    )
    expected = (len(upper) + 1) * nq
    if len(derivatives) != expected:
        raise ValueError(
            "batch differentiator returned "
            f"{len(derivatives)} values; expected {expected}"
        )

    mass_derivatives: dict[tuple[int, int, int], sp.Expr] = {}
    for pair_index, (row, column) in enumerate(upper):
        for coordinate in range(nq):
            value = derivatives[pair_index * nq + coordinate]
            mass_derivatives[coordinate, row, column] = value
            mass_derivatives[coordinate, column, row] = value

    coriolis = sp.Matrix([
        sp.Add(*(
            mass_derivatives[k, i, j] * plant.dq[j] * plant.dq[k]
            - sp.Rational(1, 2)
            * mass_derivatives[i, j, k]
            * plant.dq[j]
            * plant.dq[k]
            for j in range(nq)
            for k in range(nq)
        ))
        for i in range(nq)
    ])
    conservative = sp.Matrix(derivatives[-nq:])
    return coriolis + conservative + plant.damping * plant.dq
