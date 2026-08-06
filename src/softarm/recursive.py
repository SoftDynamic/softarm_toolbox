from __future__ import annotations

from abc import ABC, abstractmethod
from dataclasses import dataclass

import sympy as sp

from .dynamics import _base_coordinates, _mixed_linear_angular_jacobian
from .geometry import angular_jacobian, transform_rpy
from .integration import integrate_unit, unit_gauss_rule
from .modeling import DynamicsAssembler, ModelDefinition
from .models import RecursivePlant, RuntimeParameter


def _cross(left: sp.Matrix, right: sp.Matrix) -> sp.Matrix:
    return left.cross(right)


def _directional_derivative(
    matrix: sp.Matrix,
    coordinates: sp.Matrix,
    velocities: sp.Matrix,
) -> sp.Matrix:
    result = sp.zeros(*matrix.shape)
    for coordinate, velocity in zip(coordinates, velocities, strict=True):
        result += matrix.diff(coordinate) * velocity
    return result


@dataclass(frozen=True)
class SectionVariationalKernel:
    section: int
    coordinate_offset: int
    dof: int
    acceleration_symbols: sp.Matrix
    parent_velocity_symbols: sp.Matrix
    parent_acceleration_symbols: sp.Matrix
    gravity_symbols: sp.Matrix
    end_rotation: sp.Matrix
    end_position: sp.Matrix
    end_linear_jacobian: sp.Matrix
    end_angular_jacobian: sp.Matrix
    end_linear_jacobian_rate: sp.Matrix
    end_angular_jacobian_rate: sp.Matrix
    own_wrench: sp.Matrix
    own_generalized_force: sp.Matrix


class LocalVariationalKernel(ABC):
    """Build one section's kinematic jet and inverse-dynamics contribution."""

    @abstractmethod
    def angular_jacobian(
        self,
        rotation: sp.Matrix,
        coordinates: sp.Matrix,
    ) -> sp.Matrix:
        """Return the local spatial angular Jacobian in the parent frame."""

    @abstractmethod
    def inertia_rotation(
        self,
        rotation: sp.Matrix,
        coordinates: sp.Matrix,
    ) -> sp.Matrix:
        """Return the rotation used to express cross-section inertia."""

    def build(
        self,
        definition: ModelDefinition,
        section: int,
    ) -> SectionVariationalKernel:
        arm_q = definition.coordinates.arm_q
        arm_dq = definition.coordinates.arm_dq
        dof = len(arm_q) // definition.config.segments
        offset = dof * section
        local_q = arm_q[offset : offset + dof, 0]
        local_dq = arm_dq[offset : offset + dof, 0]
        local_ddq = sp.Matrix(
            sp.symbols(f"section_ddq1:{dof + 1}", real=True)
        )
        parent_velocity = sp.Matrix(
            sp.symbols("parent_vx parent_vy parent_vz parent_wx parent_wy parent_wz", real=True)
        )
        parent_acceleration = sp.Matrix(
            sp.symbols("parent_ax parent_ay parent_az parent_alphax parent_alphay parent_alphaz", real=True)
        )
        gravity = sp.Matrix(sp.symbols("local_gx local_gy local_gz", real=True))
        parent_w = parent_velocity[3:, 0]
        parent_a = parent_acceleration[:3, 0]
        parent_alpha = parent_acceleration[3:, 0]
        xi = sp.Symbol("xi", real=True, nonnegative=True)

        def jet(material_coordinate: sp.Expr):
            transform = definition.kinematics.transform(section, material_coordinate)
            rotation = transform[:3, :3]
            position = transform[:3, 3]
            linear = position.jacobian(local_q)
            angular = self.angular_jacobian(rotation, local_q)
            linear_rate = _directional_derivative(linear, local_q, local_dq)
            angular_rate = _directional_derivative(angular, local_q, local_dq)
            return rotation, position, linear, angular, linear_rate, angular_rate

        mass = definition.sections.masses[section]
        inertia = definition.sections.inertias[section]

        def integrand(material_coordinate: sp.Expr) -> tuple[sp.Matrix, sp.Matrix]:
            rotation, position, linear, angular, linear_rate, angular_rate = jet(
                material_coordinate
            )
            relative_v = linear * local_dq
            relative_w = angular * local_dq
            point_w = parent_w + relative_w
            point_a = (
                parent_a
                + _cross(parent_alpha, position)
                + _cross(parent_w, _cross(parent_w, position))
                + 2 * _cross(parent_w, relative_v)
                + linear * local_ddq
                + linear_rate * local_dq
            )
            point_alpha = (
                parent_alpha
                + _cross(parent_w, relative_w)
                + angular * local_ddq
                + angular_rate * local_dq
            )
            inertia_rotation = self.inertia_rotation(rotation, local_q)
            expressed_inertia = inertia_rotation * inertia * inertia_rotation.T
            force = mass * (point_a - gravity)
            torque = expressed_inertia * point_alpha + _cross(
                point_w, expressed_inertia * point_w
            )
            wrench = force.col_join(_cross(position, force) + torque)
            generalized = linear.T * force + angular.T * torque
            return wrench, generalized

        if definition.distributed:
            if definition.config.integration.method == "gauss":
                own_wrench = sp.zeros(6, 1)
                own_generalized = sp.zeros(dof, 1)
                for node, weight in unit_gauss_rule(definition.config.integration):
                    node_wrench, node_generalized = integrand(node)
                    own_wrench += weight * node_wrench
                    own_generalized += weight * node_generalized
            else:
                wrench, generalized = integrand(xi)
                own_wrench = integrate_unit(
                    wrench,
                    xi,
                    definition.config.integration,
                    f"recursive section {section + 1} wrench",
                )
                own_generalized = integrate_unit(
                    generalized,
                    xi,
                    definition.config.integration,
                    f"recursive section {section + 1} generalized force",
                )
        else:
            own_wrench, own_generalized = integrand(sp.Rational(1, 2))

        end = jet(sp.S.One)
        return SectionVariationalKernel(
            section=section,
            coordinate_offset=offset,
            dof=dof,
            acceleration_symbols=local_ddq,
            parent_velocity_symbols=parent_velocity,
            parent_acceleration_symbols=parent_acceleration,
            gravity_symbols=gravity,
            end_rotation=end[0],
            end_position=end[1],
            end_linear_jacobian=end[2],
            end_angular_jacobian=end[3],
            end_linear_jacobian_rate=end[4],
            end_angular_jacobian_rate=end[5],
            own_wrench=own_wrench,
            own_generalized_force=own_generalized,
        )


class ExactSE3Kernel(LocalVariationalKernel):
    """Local kernel for PCS, PAC, and Cosserat transformations in SE(3)."""

    def angular_jacobian(
        self,
        rotation: sp.Matrix,
        coordinates: sp.Matrix,
    ) -> sp.Matrix:
        return angular_jacobian(rotation, coordinates)

    def inertia_rotation(
        self,
        rotation: sp.Matrix,
        coordinates: sp.Matrix,
    ) -> sp.Matrix:
        return rotation


class LegacyRitzKernel(LocalVariationalKernel):
    """First-order adapter for the existing displacement-Ritz kinematics."""

    def angular_jacobian(
        self,
        rotation: sp.Matrix,
        coordinates: sp.Matrix,
    ) -> sp.Matrix:
        return _mixed_linear_angular_jacobian(rotation, coordinates, coordinates)

    def inertia_rotation(
        self,
        rotation: sp.Matrix,
        coordinates: sp.Matrix,
    ) -> sp.Matrix:
        return rotation.subs({coordinate: 0 for coordinate in coordinates})


def local_kernel_for(definition: ModelDefinition) -> LocalVariationalKernel:
    return LegacyRitzKernel() if definition.linear_kinematics else ExactSE3Kernel()


class RecursiveInverseDynamicsAssembler(DynamicsAssembler):
    """Create an algorithmic plant without assembling global symbolic M or h."""

    def assemble(self, definition: ModelDefinition) -> RecursivePlant:
        config = definition.config
        parameters = list(definition.parameters)

        def global_parameter(name: str, default: float) -> sp.Symbol:
            symbol = sp.Symbol(name, real=True)
            parameters.append(
                RuntimeParameter(name, symbol, float(config.parameters.get(name, default)))
            )
            return symbol

        global_parameter("gravity", 9.81)
        global_parameter("tip_mass", 0.1)
        global_parameter("tip_Ixx", 0.01)
        global_parameter("tip_Iyy", 0.01)
        global_parameter("tip_Izz", 0.01)
        base_q, base_dq = _base_coordinates(config)
        if config.base.mode == "floating_rpy":
            global_parameter("vehicle_mass", 1.5)
            global_parameter("vehicle_Ixx", 0.03)
            global_parameter("vehicle_Iyy", 0.03)
            global_parameter("vehicle_Izz", 0.05)
            base_transform = transform_rpy(base_q[:3, 0], tuple(base_q[3:6, 0]))
        else:
            base_transform = sp.eye(4)
        arm_q = definition.coordinates.arm_q
        arm_dq = definition.coordinates.arm_dq
        q = base_q.col_join(arm_q)
        dq = base_dq.col_join(arm_dq)
        nq = len(q)
        arm_force_map = sp.zeros(nq, len(arm_q))
        arm_force_map[len(base_q) :, :] = sp.eye(len(arm_q))
        return RecursivePlant(
            config=config,
            q=q,
            dq=dq,
            base_q=base_q,
            base_dq=base_dq,
            arm_q=arm_q,
            arm_dq=arm_dq,
            parameters=tuple(parameters),
            kinematics=sp.zeros(0, 0),
            end_jacobian=sp.zeros(0, nq),
            base_transform=base_transform,
            end_transform=sp.zeros(0, 0),
            base_jacobian=sp.zeros(0, nq),
            vehicle_wrench_map=sp.zeros(nq, 6),
            arm_force_map=arm_force_map,
            definition=definition,
        )
