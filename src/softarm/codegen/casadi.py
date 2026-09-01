from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import sympy as sp

from ..actuation import ActuationModel
from ..constraints import ConstraintModel
from ..geometry import angular_jacobian, transform_rpy
from ..models import PlantModel, RecursivePlant, RuntimeParameter, SymbolicPlant
from ..recursive import local_kernel_for
from ..special import (
    AffineCosMoment,
    AffineSinMoment,
    CoscSqrt,
    CoscSqrtD,
    CoscSqrtDD,
    Sinc3Sqrt,
    Sinc3SqrtD,
    Sinc3SqrtDD,
    SincSqrt,
    SincSqrtD,
    SincSqrtDD,
)


def _import_casadi():
    try:
        import casadi as ca
    except ImportError as error:
        raise RuntimeError(
            "CasADi target requires the optional dependency; "
            "install softarm-toolbox[casadi]"
        ) from error
    return ca


_PAC_HEADS = (AffineCosMoment, AffineSinMoment)


class SympyToCasadi:
    """Translate the supported SymPy DAG to scalar CasADi SX expressions."""

    def __init__(self, symbols: dict[sp.Symbol, Any], affine_terms: int | None):
        self.ca = _import_casadi()
        self.symbols = symbols
        self.affine_terms = affine_terms
        self.cache: dict[sp.Expr, Any] = {}
        self._affine_cache: dict[tuple[int, sp.Expr, sp.Expr, sp.Expr], tuple[Any, Any]] = {}

    def matrix(self, matrix: sp.MatrixBase):
        values = [
            self.expression(matrix[row, column])
            for column in range(matrix.cols)
            for row in range(matrix.rows)
        ]
        if not values:
            return self.ca.SX.zeros(matrix.rows, matrix.cols)
        return self.ca.reshape(self.ca.vertcat(*values), matrix.rows, matrix.cols)

    def expression(self, expression: sp.Expr):
        cached = self.cache.get(expression)
        if cached is not None:
            return cached
        result = self._expression(expression)
        self.cache[expression] = result
        return result

    def _expression(self, expression: sp.Expr):
        ca = self.ca
        if expression in self.symbols:
            return self.symbols[expression]
        if expression is sp.S.Pi:
            return ca.pi
        if expression is sp.S.Exp1:
            return ca.exp(1)
        if isinstance(expression, sp.Integer):
            return int(expression)
        if isinstance(expression, sp.Rational):
            return int(expression.p) / int(expression.q)
        if isinstance(expression, sp.Float):
            return float(expression)
        if isinstance(expression, sp.Add):
            terms = [self.expression(item) for item in expression.args]
            return sum(terms[1:], terms[0])
        if isinstance(expression, sp.Mul):
            factors = [self.expression(item) for item in expression.args]
            result = factors[0]
            for factor in factors[1:]:
                result *= factor
            return result
        if isinstance(expression, sp.Pow):
            return self.expression(expression.base) ** self.expression(expression.exp)

        unary = {
            sp.sin: ca.sin,
            sp.cos: ca.cos,
            sp.tan: ca.tan,
            sp.asin: ca.asin,
            sp.acos: ca.acos,
            sp.atan: ca.atan,
            sp.sinh: ca.sinh,
            sp.cosh: ca.cosh,
            sp.tanh: ca.tanh,
            sp.exp: ca.exp,
            sp.log: ca.log,
            sp.Abs: ca.fabs,
            sp.sign: ca.sign,
        }
        operation = unary.get(expression.func)
        if operation is not None:
            return operation(self.expression(expression.args[0]))
        if expression.func is sp.atan2:
            return ca.atan2(*(self.expression(item) for item in expression.args))

        z = self.expression(expression.args[0]) if len(expression.args) == 1 else None
        if isinstance(expression, SincSqrt):
            return ca.if_else(
                ca.fabs(z) < 1e-8,
                1 - z / 6 + z**2 / 120 - z**3 / 5040 + z**4 / 362880,
                ca.sin(ca.sqrt(z)) / ca.sqrt(z),
            )
        if isinstance(expression, SincSqrtD):
            root = ca.sqrt(z)
            return ca.if_else(
                ca.fabs(z) < 1e-8,
                -1 / 6 + z / 60 - z**2 / 1680 + z**3 / 90720,
                (root * ca.cos(root) - ca.sin(root)) / (2 * root**3),
            )
        if isinstance(expression, SincSqrtDD):
            root = ca.sqrt(z)
            return ca.if_else(
                ca.fabs(z) < 1e-8,
                1 / 60 - z / 840 + z**2 / 30240,
                ((3 - z) * ca.sin(root) - 3 * root * ca.cos(root)) / (4 * root**5),
            )
        if isinstance(expression, CoscSqrt):
            return ca.if_else(
                ca.fabs(z) < 1e-8,
                1 / 2 - z / 24 + z**2 / 720 - z**3 / 40320 + z**4 / 3628800,
                (1 - ca.cos(ca.sqrt(z))) / z,
            )
        if isinstance(expression, CoscSqrtD):
            root = ca.sqrt(z)
            return ca.if_else(
                ca.fabs(z) < 1e-8,
                -1 / 24 + z / 360 - z**2 / 13440 + z**3 / 907200,
                (root * ca.sin(root) - 2 * (1 - ca.cos(root))) / (2 * z**2),
            )
        if isinstance(expression, CoscSqrtDD):
            root = ca.sqrt(z)
            return ca.if_else(
                ca.fabs(z) < 1e-8,
                1 / 360 - z / 6720 + z**2 / 302400,
                (z * ca.cos(root) - 5 * root * ca.sin(root) + 8 - 8 * ca.cos(root))
                / (4 * z**3),
            )
        if isinstance(expression, Sinc3Sqrt):
            sinc = self.expression(SincSqrt(expression.args[0]))
            return ca.if_else(
                ca.fabs(z) < 1e-8,
                1 / 6 - z / 120 + z**2 / 5040 - z**3 / 362880 + z**4 / 39916800,
                (1 - sinc) / z,
            )
        if isinstance(expression, Sinc3SqrtD):
            sinc = self.expression(SincSqrt(expression.args[0]))
            sinc_d = self.expression(SincSqrtD(expression.args[0]))
            return ca.if_else(
                ca.fabs(z) < 1e-8,
                -1 / 120 + z / 2520 - z**2 / 120960 + z**3 / 9979200,
                (sinc - 1 - z * sinc_d) / z**2,
            )
        if isinstance(expression, Sinc3SqrtDD):
            sinc = self.expression(SincSqrt(expression.args[0]))
            sinc_d = self.expression(SincSqrtD(expression.args[0]))
            sinc_dd = self.expression(SincSqrtDD(expression.args[0]))
            return ca.if_else(
                ca.fabs(z) < 1e-8,
                1 / 2520 - z / 60480 + z**2 / 3326400,
                (2 - 2 * sinc + 2 * z * sinc_d - z**2 * sinc_dd) / z**3,
            )
        if isinstance(expression, _PAC_HEADS):
            order, c0, c1, xi = expression.args
            real, imag = self._affine_moment(
                int(order), c0, c1, xi
            )
            return real if isinstance(expression, AffineCosMoment) else imag
        raise TypeError(
            f"unsupported SymPy node {type(expression).__name__}: {expression}"
        )

    def _affine_moment(
        self,
        order: int,
        c0_expr: sp.Expr,
        c1_expr: sp.Expr,
        xi_expr: sp.Expr,
    ) -> tuple[Any, Any]:
        if self.affine_terms is None:
            raise ValueError(
                "PAC CasADi export requires an explicit affine_terms value"
            )
        key = order, c0_expr, c1_expr, xi_expr
        cached = self._affine_cache.get(key)
        if cached is not None:
            return cached
        ca = self.ca
        c0 = self.expression(c0_expr)
        c1 = self.expression(c1_expr)
        xi = self.expression(xi_expr)
        state = ca.SX.sym("affine_state", 6)
        data = ca.SX.sym("affine_data", 5)
        degree, step_c0, step_c1, step_xi, step_order = ca.vertsplit(data)
        pp_real, pp_imag, p_real, p_imag, sum_real, sum_imag = ca.vertsplit(state)
        coefficient_real = (-step_c0 * p_imag - step_c1 * pp_imag) / degree
        coefficient_imag = (step_c0 * p_real + step_c1 * pp_real) / degree
        scale = step_xi ** (degree + step_order + 1) / (degree + step_order + 1)
        next_state = ca.vertcat(
            p_real,
            p_imag,
            coefficient_real,
            coefficient_imag,
            sum_real + coefficient_real * scale,
            sum_imag + coefficient_imag * scale,
        )
        step = ca.Function("softarm_affine_step", [state, data], [next_state])
        accumulated = step.mapaccum(
            f"softarm_affine_accum_{self.affine_terms}", self.affine_terms
        )
        initial = ca.vertcat(0, 0, 1, 0, xi ** (order + 1) / (order + 1), 0)
        columns = [
            ca.vertcat(index, c0, c1, xi, order)
            for index in range(1, self.affine_terms + 1)
        ]
        result = accumulated(initial, ca.horzcat(*columns))
        final = result[:, -1]
        values = final[4], final[5]
        self._affine_cache[key] = values
        return values


@dataclass
class _FunctionRegistry:
    target: Path
    functions: dict[str, dict[str, Any]]

    def add(self, logical_name: str, function) -> None:
        relative = Path("functions", logical_name.replace(".", "_") + ".casadi")
        path = self.target / relative
        function.save(str(path))
        self.functions[logical_name] = {
            "file": relative.as_posix(),
            "inputs": [
                {
                    "name": function.name_in(index),
                    "shape": [function.size1_in(index), function.size2_in(index)],
                }
                for index in range(function.n_in())
            ],
            "outputs": [
                {
                    "name": function.name_out(index),
                    "shape": [function.size1_out(index), function.size2_out(index)],
                }
                for index in range(function.n_out())
            ],
        }


def _sympy_function(
    name: str,
    inputs: list[tuple[str, sp.Matrix]],
    outputs: list[tuple[str, sp.Matrix]],
    affine_terms: int | None,
):
    ca = _import_casadi()
    symbolic_inputs = [ca.SX.sym(input_name, len(symbols)) for input_name, symbols in inputs]
    symbol_map = {
        symbol: symbolic_inputs[input_index][symbol_index]
        for input_index, (_, symbols) in enumerate(inputs)
        for symbol_index, symbol in enumerate(symbols)
    }
    converter = SympyToCasadi(symbol_map, affine_terms)
    expressions = [converter.matrix(matrix) for _, matrix in outputs]
    return ca.Function(
        name,
        symbolic_inputs,
        expressions,
        [item[0] for item in inputs],
        [item[0] for item in outputs],
    )


def _runtime_parameters(
    plant: PlantModel,
    actuation: ActuationModel | None,
    constraint: ConstraintModel | None,
) -> tuple[RuntimeParameter, ...]:
    return (
        plant.parameters
        + (() if actuation is None else actuation.parameters)
        + (() if constraint is None else constraint.parameters)
    )


def _directional_sympy(matrix: sp.Matrix, q: sp.Matrix, dq: sp.Matrix) -> sp.Matrix:
    result = sp.zeros(*matrix.shape)
    for coordinate, velocity in zip(q, dq, strict=True):
        result += matrix.diff(coordinate) * velocity
    return result


def _parameter(plant: RecursivePlant, name: str) -> sp.Symbol:
    return next(item.symbol for item in plant.parameters if item.name == name)


def _symbolic_core(
    plant: SymbolicPlant,
    actuation: ActuationModel | None,
    constraint: ConstraintModel | None,
    parameters: tuple[RuntimeParameter, ...],
    affine_terms: int | None,
) -> dict[str, Any]:
    p = sp.Matrix([item.symbol for item in parameters])
    q = plant.q
    dq = plant.dq
    functions = {
        "mass": _sympy_function(
            "softarm_mass", [("q", q), ("p", p)], [("M", plant.mass)], affine_terms
        ),
        "bias": _sympy_function(
            "softarm_bias",
            [("q", q), ("dq", dq), ("p", p)],
            [("h", plant.bias)],
            affine_terms,
        ),
        "kinematics": _sympy_function(
            "softarm_kinematics",
            [("q", q), ("p", p)],
            [("H", plant.kinematics)],
            affine_terms,
        ),
        "end_jacobian": _sympy_function(
            "softarm_end_jacobian",
            [("q", q), ("p", p)],
            [("J", plant.end_jacobian)],
            affine_terms,
        ),
        "vehicle_wrench_map": _sympy_function(
            "softarm_vehicle_wrench_map",
            [("q", q), ("p", p)],
            [("Bv", plant.vehicle_wrench_map)],
            affine_terms,
        ),
    }
    if plant._material_coordinate is not None and plant._material_kinematics is not None:
        xi = sp.Matrix([plant._material_coordinate])
        positions = sp.Matrix.hstack(*(
            plant._material_kinematics[:3, 4 * section + 3]
            for section in range(plant.config.segments)
        ))
        functions["centerline_at"] = _sympy_function(
            "softarm_centerline_at",
            [("q", q), ("p", p), ("xi", xi)],
            [("position", positions)],
            affine_terms,
        )
    if actuation is not None:
        functions.update({
            "actuator_coordinates": _sympy_function(
                "softarm_actuator_coordinates",
                [("q", q), ("p", p)],
                [("y", actuation.coordinates)],
                affine_terms,
            ),
            "actuator_jacobian": _sympy_function(
                "softarm_actuator_jacobian",
                [("q", q), ("p", p)],
                [("Ja", actuation.jacobian)],
                affine_terms,
            ),
            "actuator_velocity_bias": _sympy_function(
                "softarm_actuator_velocity_bias",
                [("q", q), ("dq", dq), ("p", p)],
                [("gamma", actuation.velocity_bias)],
                affine_terms,
            ),
        })
    if constraint is not None:
        functions.update({
            "constraint_value": _sympy_function(
                "softarm_constraint_value",
                [("q", q), ("p", p)],
                [("phi", constraint.coordinates)],
                affine_terms,
            ),
            "constraint_jacobian": _sympy_function(
                "softarm_constraint_jacobian",
                [("q", q), ("p", p)],
                [("A", constraint.jacobian)],
                affine_terms,
            ),
            "constraint_velocity_bias": _sympy_function(
                "softarm_constraint_velocity_bias",
                [("q", q), ("dq", dq), ("p", p)],
                [("gamma", constraint.velocity_bias)],
                affine_terms,
            ),
            "constraint_reaction_map": _sympy_function(
                "softarm_constraint_reaction_map",
                [("q", q), ("dq", dq), ("p", p)],
                [("G", constraint.reaction_map)],
                affine_terms,
            ),
            "constraint_stabilization": _sympy_function(
                "softarm_constraint_stabilization",
                [("p", p)],
                [(
                    "gains",
                    constraint.stabilization_frequency.row_join(
                        constraint.stabilization_ratio
                    ),
                )],
                affine_terms,
            ),
        })
    return functions


def _recursive_core(
    plant: RecursivePlant,
    actuation: ActuationModel | None,
    parameters: tuple[RuntimeParameter, ...],
    affine_terms: int | None,
) -> dict[str, Any]:
    ca = _import_casadi()
    definition = plant.definition
    nq = len(plant.q)
    nbase = len(plant.base_q)
    narm = len(plant.arm_q)
    segments = plant.config.segments
    dof = narm // segments
    p_symbols = sp.Matrix([item.symbol for item in parameters])

    base_ddq = (
        sp.Matrix(sp.symbols(f"base_ddq1:{nbase + 1}", real=True))
        if nbase
        else sp.zeros(0, 1)
    )
    gravity_world = sp.Matrix([0, 0, _parameter(plant, "gravity")])
    mount = transform_rpy(
        tuple(sp.Float(str(value)) for value in plant.config.base.mount_xyz),
        tuple(sp.Float(str(value)) for value in plant.config.base.mount_rpy),
    )
    mount_transform = plant.base_transform * mount
    mount_rotation = mount_transform[:3, :3]
    mount_position = mount_transform[:3, 3]
    if nbase:
        mount_linear = mount_position.jacobian(plant.base_q)
        mount_angular = angular_jacobian(mount_rotation, plant.base_q)
        mount_linear_rate = _directional_sympy(
            mount_linear, plant.base_q, plant.base_dq
        )
        mount_angular_rate = _directional_sympy(
            mount_angular, plant.base_q, plant.base_dq
        )
        mount_v_world = mount_linear * plant.base_dq
        mount_w_world = mount_angular * plant.base_dq
        mount_a_world = mount_linear * base_ddq + mount_linear_rate * plant.base_dq
        mount_alpha_world = (
            mount_angular * base_ddq + mount_angular_rate * plant.base_dq
        )
        parent_v = mount_rotation.T * mount_v_world
        parent_w = mount_rotation.T * mount_w_world
        parent_a = mount_rotation.T * mount_a_world
        parent_alpha = mount_rotation.T * mount_alpha_world
        local_gravity = mount_rotation.T * gravity_world
        mount_map = mount_linear.col_join(mount_angular).T * sp.diag(
            mount_rotation, mount_rotation
        )

        base_position = plant.base_transform[:3, 3]
        base_rotation = plant.base_transform[:3, :3]
        base_linear = base_position.jacobian(plant.base_q)
        base_angular = angular_jacobian(base_rotation, plant.base_q)
        base_linear_rate = _directional_sympy(
            base_linear, plant.base_q, plant.base_dq
        )
        base_angular_rate = _directional_sympy(
            base_angular, plant.base_q, plant.base_dq
        )
        base_w = base_angular * plant.base_dq
        base_a = base_linear * base_ddq + base_linear_rate * plant.base_dq
        base_alpha = base_angular * base_ddq + base_angular_rate * plant.base_dq
        vehicle_inertia = sp.diag(
            _parameter(plant, "vehicle_Ixx"),
            _parameter(plant, "vehicle_Iyy"),
            _parameter(plant, "vehicle_Izz"),
        )
        inertia_world = base_rotation * vehicle_inertia * base_rotation.T
        vehicle_force = _parameter(plant, "vehicle_mass") * (base_a - gravity_world)
        vehicle_torque = inertia_world * base_alpha + base_w.cross(
            inertia_world * base_w
        )
        base_tau = base_linear.T * vehicle_force + base_angular.T * vehicle_torque
        vehicle_map_base = base_linear.col_join(base_angular).T * sp.diag(
            base_rotation, base_rotation
        )
    else:
        parent_v = parent_w = parent_a = parent_alpha = sp.zeros(3, 1)
        local_gravity = mount_rotation.T * gravity_world
        mount_map = sp.zeros(0, 6)
        base_tau = sp.zeros(0, 1)
        vehicle_map_base = sp.zeros(0, 6)

    base_function = _sympy_function(
        "softarm_recursive_base",
        [("q", plant.q), ("dq", plant.dq), ("ddq_base", base_ddq), ("p", p_symbols)],
        [
            ("v", parent_v),
            ("w", parent_w),
            ("a", parent_a),
            ("alpha", parent_alpha),
            ("g", local_gravity),
            ("base_tau", base_tau),
            ("mount_map", mount_map),
        ],
        affine_terms,
    )
    mount_function = _sympy_function(
        "softarm_mount_transform",
        [("q", plant.q), ("p", p_symbols)],
        [("H", mount_transform)],
        affine_terms,
    )
    full_vehicle_map = sp.zeros(nq, 6)
    if nbase:
        full_vehicle_map[:nbase, :] = vehicle_map_base
    vehicle_map_function = _sympy_function(
        "softarm_vehicle_wrench_map",
        [("q", plant.q), ("p", p_symbols)],
        [("Bv", full_vehicle_map)],
        affine_terms,
    )

    kernel = local_kernel_for(definition).build(definition, 0)
    local_q_source = plant.arm_q[:dof, 0]
    local_dq_source = plant.arm_dq[:dof, 0]
    section_parameters = [
        item for item in plant.parameters if item.name.startswith("s1_")
    ]
    generic_q = sp.Matrix(sp.symbols(f"local_q1:{dof + 1}", real=True))
    generic_dq = sp.Matrix(sp.symbols(f"local_dq1:{dof + 1}", real=True))
    generic_p = sp.Matrix(
        sp.symbols(f"local_p1:{len(section_parameters) + 1}", real=True)
    )
    substitutions = dict(zip(local_q_source, generic_q, strict=True))
    substitutions.update(zip(local_dq_source, generic_dq, strict=True))
    substitutions.update(
        zip((item.symbol for item in section_parameters), generic_p, strict=True)
    )
    section_function = _sympy_function(
        "softarm_recursive_section_template",
        [
            ("q_local", generic_q),
            ("dq_local", generic_dq),
            ("ddq_local", kernel.acceleration_symbols),
            ("p_local", generic_p),
            ("v0", kernel.parent_velocity_symbols[:3, 0]),
            ("w0", kernel.parent_velocity_symbols[3:, 0]),
            ("a0", kernel.parent_acceleration_symbols[:3, 0]),
            ("alpha0", kernel.parent_acceleration_symbols[3:, 0]),
            ("g0", kernel.gravity_symbols),
        ],
        [
            ("own_wrench", kernel.own_wrench.xreplace(substitutions)),
            ("own_tau", kernel.own_generalized_force.xreplace(substitutions)),
        ],
        affine_terms,
    )
    end_function = _sympy_function(
        "softarm_recursive_end_template",
        [("q_local", generic_q), ("dq_local", generic_dq), ("p_local", generic_p)],
        [
            ("rotation", kernel.end_rotation.xreplace(substitutions)),
            ("position", kernel.end_position.xreplace(substitutions)),
            ("Jv", kernel.end_linear_jacobian.xreplace(substitutions)),
            ("Jw", kernel.end_angular_jacobian.xreplace(substitutions)),
            ("Jvd", kernel.end_linear_jacobian_rate.xreplace(substitutions)),
            ("Jwd", kernel.end_angular_jacobian_rate.xreplace(substitutions)),
        ],
        affine_terms,
    )
    xi_symbol = sp.Symbol("xi", real=True, nonnegative=True)
    pose = definition.kinematics.transform(0, xi_symbol).xreplace(substitutions)
    pose_function = _sympy_function(
        "softarm_recursive_pose_template",
        [("q_local", generic_q), ("p_local", generic_p), ("xi", sp.Matrix([xi_symbol]))],
        [("H", pose)],
        affine_terms,
    )

    arm_internal = sp.Matrix(
        [sp.diff(definition.elastic, coordinate) for coordinate in plant.arm_q]
    )
    damping = (
        sp.diag(*definition.damping)
        if isinstance(definition.damping, tuple)
        else definition.damping
    )
    arm_internal += damping * plant.arm_dq
    internal = sp.zeros(nq, 1)
    internal[nbase:, 0] = arm_internal
    internal_function = _sympy_function(
        "softarm_recursive_internal_force",
        [("q", plant.q), ("dq", plant.dq), ("p", p_symbols)],
        [("tau", internal)],
        affine_terms,
    )
    tip_velocity = sp.Matrix(
        sp.symbols("tip_vx tip_vy tip_vz tip_wx tip_wy tip_wz", real=True)
    )
    tip_acceleration = sp.Matrix(
        sp.symbols("tip_ax tip_ay tip_az tip_alphax tip_alphay tip_alphaz", real=True)
    )
    tip_gravity = sp.Matrix(sp.symbols("tip_gx tip_gy tip_gz", real=True))
    tip_inertia = sp.diag(
        _parameter(plant, "tip_Ixx"),
        _parameter(plant, "tip_Iyy"),
        _parameter(plant, "tip_Izz"),
    )
    tip_force = _parameter(plant, "tip_mass") * (
        tip_acceleration[:3, 0] - tip_gravity
    )
    tip_torque = tip_inertia * tip_acceleration[3:, 0] + tip_velocity[3:, 0].cross(
        tip_inertia * tip_velocity[3:, 0]
    )
    tip_function = _sympy_function(
        "softarm_recursive_tip",
        [
            ("velocity", tip_velocity),
            ("acceleration", tip_acceleration),
            ("gravity", tip_gravity),
            ("p", p_symbols),
        ],
        [("wrench", tip_force.col_join(tip_torque))],
        affine_terms,
    )

    parameter_indices = {item.name: index for index, item in enumerate(parameters)}
    suffixes = [item.name.removeprefix("s1_") for item in section_parameters]

    def local_values(q, dq, section: int):
        first = nbase + section * dof
        indices = [
            parameter_indices[f"s{section + 1}_{suffix}"] for suffix in suffixes
        ]
        return q[first : first + dof], dq[first : first + dof], indices

    q_mx = ca.MX.sym("q", nq)
    dq_mx = ca.MX.sym("dq", nq)
    ddq_mx = ca.MX.sym("ddq", nq)
    p_mx = ca.MX.sym("p", len(parameters))
    v, w, a, alpha, gravity, base_tau_mx, mount_map_mx = base_function(
        q_mx, dq_mx, ddq_mx[:nbase], p_mx
    )
    forward_state = ca.MX.sym("forward_state", 15)
    q_local_step = ca.MX.sym("q_local", dof)
    dq_local_step = ca.MX.sym("dq_local", dof)
    ddq_local_step = ca.MX.sym("ddq_local", dof)
    p_local_step = ca.MX.sym("p_local", len(section_parameters))
    step_v = forward_state[0:3]
    step_w = forward_state[3:6]
    step_a = forward_state[6:9]
    step_alpha = forward_state[9:12]
    step_gravity = forward_state[12:15]
    rotation, position, jv, jw, jvd, jwd = end_function(
        q_local_step, dq_local_step, p_local_step
    )
    own_wrench, own_tau = section_function(
        q_local_step,
        dq_local_step,
        ddq_local_step,
        p_local_step,
        step_v,
        step_w,
        step_a,
        step_alpha,
        step_gravity,
    )
    relative_v = jv @ dq_local_step
    relative_w = jw @ dq_local_step
    v_parent = step_v + ca.cross(step_w, position) + relative_v
    w_parent = step_w + relative_w
    a_parent = (
        step_a
        + ca.cross(step_alpha, position)
        + ca.cross(step_w, ca.cross(step_w, position))
        + 2 * ca.cross(step_w, relative_v)
        + jv @ ddq_local_step
        + jvd @ dq_local_step
    )
    alpha_parent = (
        step_alpha
        + ca.cross(step_w, relative_w)
        + jw @ ddq_local_step
        + jwd @ dq_local_step
    )
    next_state = ca.vertcat(
        rotation.T @ v_parent,
        rotation.T @ w_parent,
        rotation.T @ a_parent,
        rotation.T @ alpha_parent,
        rotation.T @ step_gravity,
    )
    section_data = ca.vertcat(
        ca.reshape(rotation, 9, 1),
        position,
        ca.reshape(jv, 3 * dof, 1),
        ca.reshape(jw, 3 * dof, 1),
        own_wrench,
        own_tau,
    )
    forward_step = ca.Function(
        "softarm_recursive_forward_step",
        [forward_state, q_local_step, dq_local_step, ddq_local_step, p_local_step],
        [next_state, section_data],
    )
    forward_map = forward_step.mapaccum(
        "softarm_recursive_forward", segments, {"base": 10}
    )
    q_columns = []
    dq_columns = []
    ddq_columns = []
    p_columns = []
    for section in range(segments):
        q_local, dq_local, indices = local_values(q_mx, dq_mx, section)
        q_columns.append(q_local)
        dq_columns.append(dq_local)
        ddq_columns.append(
            ddq_mx[nbase + section * dof : nbase + (section + 1) * dof]
        )
        p_columns.append(p_mx[indices])
    state_history, data_history = forward_map(
        ca.vertcat(v, w, a, alpha, gravity),
        ca.horzcat(*q_columns),
        ca.horzcat(*dq_columns),
        ca.horzcat(*ddq_columns),
        ca.horzcat(*p_columns),
    )
    final_state = state_history[:, -1]
    wrench = tip_function(
        final_state[:6], final_state[6:12], final_state[12:15], p_mx
    )

    backward_wrench = ca.MX.sym("backward_wrench", 6)
    backward_data = ca.MX.sym("backward_data", section_data.size1())
    cursor = 0
    backward_rotation = ca.reshape(backward_data[cursor : cursor + 9], 3, 3)
    cursor += 9
    backward_position = backward_data[cursor : cursor + 3]
    cursor += 3
    backward_jv = ca.reshape(backward_data[cursor : cursor + 3 * dof], 3, dof)
    cursor += 3 * dof
    backward_jw = ca.reshape(backward_data[cursor : cursor + 3 * dof], 3, dof)
    cursor += 3 * dof
    backward_own_wrench = backward_data[cursor : cursor + 6]
    cursor += 6
    backward_own_tau = backward_data[cursor : cursor + dof]
    force = backward_rotation @ backward_wrench[:3]
    moment = backward_rotation @ backward_wrench[3:]
    backward_tau = backward_own_tau + backward_jv.T @ force + backward_jw.T @ moment
    next_wrench = backward_own_wrench + ca.vertcat(
        force, ca.cross(backward_position, force) + moment
    )
    backward_step = ca.Function(
        "softarm_recursive_backward_step",
        [backward_wrench, backward_data],
        [next_wrench, backward_tau],
    )
    backward_map = backward_step.mapaccum(
        "softarm_recursive_backward", segments, {"base": 10}
    )
    reversed_data = ca.horzcat(
        *(data_history[:, section] for section in reversed(range(segments)))
    )
    wrench_history, reversed_tau = backward_map(wrench, reversed_data)
    root_wrench = wrench_history[:, -1]
    tau_sections = ca.vertcat(
        *(reversed_tau[:, section] for section in reversed(range(segments)))
    )
    inverse = ca.vertcat(base_tau_mx + mount_map_mx @ root_wrench, tau_sections)
    inverse += internal_function(q_mx, dq_mx, p_mx)
    inverse_function = ca.Function(
        "softarm_inverse_dynamics",
        [q_mx, dq_mx, ddq_mx, p_mx],
        [inverse],
        ["q", "dq", "ddq", "p"],
        ["tau"],
    )

    zero_dq = ca.MX.zeros(nq, 1)
    probe_ddq = ca.MX.sym("mass_ddq", nq)
    inverse_probe = inverse_function(q_mx, zero_dq, probe_ddq, p_mx)
    mass = ca.jacobian(inverse_probe, probe_ddq)
    mass = ca.substitute(mass, probe_ddq, ca.MX.zeros(nq, 1))
    mass = (mass + mass.T) / 2
    mass_function = ca.Function(
        "softarm_mass", [q_mx, p_mx], [mass], ["q", "p"], ["M"]
    )
    bias_function = ca.Function(
        "softarm_bias",
        [q_mx, dq_mx, p_mx],
        [inverse_function(q_mx, dq_mx, ca.MX.zeros(nq, 1), p_mx)],
        ["q", "dq", "p"],
        ["h"],
    )

    transforms = []
    current = mount_function(q_mx, p_mx)
    for section in range(segments):
        q_local, _, indices = local_values(q_mx, zero_dq, section)
        current = current @ pose_function(q_local, p_mx[indices], 1)
        transforms.append(current)
    kinematics_function = ca.Function(
        "softarm_kinematics",
        [q_mx, p_mx],
        [ca.horzcat(*transforms)],
        ["q", "p"],
        ["H"],
    )
    xi_mx = ca.MX.sym("xi")
    points = []
    current = mount_function(q_mx, p_mx)
    for section in range(segments):
        q_local, _, indices = local_values(q_mx, zero_dq, section)
        material = current @ pose_function(q_local, p_mx[indices], xi_mx)
        points.append(material[:3, 3])
        current = current @ pose_function(q_local, p_mx[indices], 1)
    centerline_function = ca.Function(
        "softarm_centerline_at",
        [q_mx, p_mx, xi_mx],
        [ca.horzcat(*points)],
        ["q", "p", "xi"],
        ["position"],
    )

    end_transform = transforms[-1]
    end_position_world = end_transform[:3, 3]
    end_rotation_world = end_transform[:3, :3]
    linear_jacobian = ca.jacobian(end_position_world, q_mx)
    rotation_jacobian = ca.jacobian(ca.reshape(end_rotation_world, 9, 1), q_mx)
    angular_columns = []
    for index in range(nq):
        rate = ca.reshape(rotation_jacobian[:, index], 3, 3) @ end_rotation_world.T
        skew = (rate - rate.T) / 2
        angular_columns.append(ca.vertcat(skew[2, 1], skew[0, 2], skew[1, 0]))
    end_jacobian_function = ca.Function(
        "softarm_end_jacobian",
        [q_mx, p_mx],
        [ca.vertcat(linear_jacobian, ca.horzcat(*angular_columns))],
        ["q", "p"],
        ["J"],
    )

    core = {
        "inverse_dynamics": inverse_function,
        "mass": mass_function,
        "bias": bias_function,
        "kinematics": kinematics_function,
        "centerline_at": centerline_function,
        "end_jacobian": end_jacobian_function,
        "vehicle_wrench_map": vehicle_map_function,
    }
    if actuation is not None:
        core.update({
            "actuator_coordinates": _sympy_function(
                "softarm_actuator_coordinates",
                [("q", plant.q), ("p", p_symbols)],
                [("y", actuation.coordinates)],
                affine_terms,
            ),
            "actuator_jacobian": _sympy_function(
                "softarm_actuator_jacobian",
                [("q", plant.q), ("p", p_symbols)],
                [("Ja", actuation.jacobian)],
                affine_terms,
            ),
            "actuator_velocity_bias": _sympy_function(
                "softarm_actuator_velocity_bias",
                [("q", plant.q), ("dq", plant.dq), ("p", p_symbols)],
                [("gamma", actuation.velocity_bias)],
                affine_terms,
            ),
        })
    return core


def _applied_force_function(plant: PlantModel, core: dict[str, Any], parameter_count: int):
    ca = _import_casadi()
    q = ca.MX.sym("q", len(plant.q))
    tau_arm = ca.MX.sym("tau_arm", len(plant.arm_q))
    w_vehicle = ca.MX.sym("w_vehicle", 6)
    w_tip = ca.MX.sym("w_tip", 6)
    p = ca.MX.sym("p", parameter_count)
    arm_force = ca.vertcat(ca.MX.zeros(len(plant.base_q), 1), tau_arm)
    force = (
        arm_force
        + core["vehicle_wrench_map"](q, p)
        @ w_vehicle
        + core["end_jacobian"](q, p).T
        @ w_tip
    )
    return ca.Function(
        "softarm_applied_force",
        [q, tau_arm, w_vehicle, w_tip, p],
        [force],
        ["q", "tau_arm", "w_vehicle", "w_tip", "p"],
        ["Q"],
    )


def _dynamics_functions(
    plant: PlantModel,
    core: dict[str, Any],
    actuation: ActuationModel | None,
    constraint: ConstraintModel | None,
    parameter_count: int,
) -> tuple[dict[str, Any], dict[str, Any]]:
    ca = _import_casadi()
    nq = len(plant.q)
    narm = len(plant.arm_q)
    nbase = len(plant.base_q)
    applied = core["applied_force"]

    def mechanical_terms(q, dq, ddq, p):
        if "inverse_dynamics" in core:
            return core["inverse_dynamics"](q, dq, ddq, p)
        return core["mass"](q, p) @ ddq + core["bias"](q, dq, p)

    functions: dict[str, Any] = {}
    variants: dict[str, Any] = {}

    def add_variant(
        logical: str,
        u_blocks: list[tuple[str, int]],
        explicit_ddq,
        implicit_mechanics,
        z_blocks: list[tuple[str, int]] | None = None,
    ) -> None:
        x = ca.MX.sym("x", 2 * nq)
        xdot = ca.MX.sym("xdot", 2 * nq)
        u = ca.MX.sym("u", sum(size for _, size in u_blocks))
        p = ca.MX.sym("p", parameter_count)
        q, dq = x[:nq], x[nq:]
        ddq = explicit_ddq(q, dq, u, p)
        functions[f"dynamics.{logical}.explicit"] = ca.Function(
            f"softarm_{logical}_f_expl",
            [x, u, p],
            [ca.vertcat(dq, ddq)],
            ["x", "u", "p"],
            ["f_expl"],
        )
        z_size = sum(size for _, size in (z_blocks or []))
        z = ca.MX.sym("z", z_size)
        residual = ca.vertcat(
            xdot[:nq] - dq,
            implicit_mechanics(q, dq, xdot[nq:], u, z, p),
        )
        functions[f"dynamics.{logical}.implicit"] = ca.Function(
            f"softarm_{logical}_f_impl",
            [xdot, x, u, z, p],
            [residual],
            ["xdot", "x", "u", "z", "p"],
            ["f_impl"],
        )
        variants[logical] = {
            "state": {"q": nq, "dq": nq},
            "state_derivative": {"qdot": nq, "ddq": nq},
            "control": [{"name": name, "size": size} for name, size in u_blocks],
            "algebraic": [
                {"name": name, "size": size} for name, size in (z_blocks or [])
            ],
            "explicit": f"dynamics.{logical}.explicit",
            "implicit": f"dynamics.{logical}.implicit",
        }

    generalized_blocks = [("tau_arm", narm), ("w_vehicle", 6), ("w_tip", 6)]

    def generalized_inputs(q, u, p):
        cursor = 0
        tau = u[cursor : cursor + narm]
        cursor += narm
        vehicle = u[cursor : cursor + 6]
        tip = u[cursor + 6 : cursor + 12]
        return applied(q, tau, vehicle, tip, p)

    def generalized_explicit(q, dq, u, p):
        rhs = generalized_inputs(q, u, p) - core["bias"](q, dq, p)
        return ca.solve(core["mass"](q, p), rhs)

    def generalized_implicit(q, dq, ddq, u, z, p):
        del z
        return mechanical_terms(q, dq, ddq, p) - generalized_inputs(q, u, p)

    add_variant(
        "generalized_force",
        generalized_blocks,
        generalized_explicit,
        generalized_implicit,
    )

    if actuation is not None:
        tendon_blocks = [
            ("tension", actuation.count),
            ("tau_arm_external", narm),
            ("w_vehicle", 6),
            ("w_tip", 6),
        ]

        def tendon_force_inputs(q, u, p):
            cursor = actuation.count
            tension = u[:cursor]
            external = u[cursor : cursor + narm]
            cursor += narm
            vehicle = u[cursor : cursor + 6]
            tip = u[cursor + 6 : cursor + 12]
            tau = external - core["actuator_jacobian"](q, p).T @ tension
            return applied(q, tau, vehicle, tip, p)

        def tendon_explicit(q, dq, u, p):
            rhs = tendon_force_inputs(q, u, p) - core["bias"](q, dq, p)
            return ca.solve(core["mass"](q, p), rhs)

        def tendon_implicit(q, dq, ddq, u, z, p):
            del z
            return mechanical_terms(q, dq, ddq, p) - tendon_force_inputs(q, u, p)

        add_variant("tendon_force", tendon_blocks, tendon_explicit, tendon_implicit)

        if actuation.acceleration == "strict":
            strict_blocks = [
                ("coordinate_acceleration", actuation.count),
                ("tau_arm_external", narm),
                ("w_vehicle", 6),
                ("w_tip", 6),
            ]

            def strict_parts(q, dq, u, p):
                cursor = actuation.count
                command = u[:cursor]
                external = u[cursor : cursor + narm]
                cursor += narm
                vehicle = u[cursor : cursor + 6]
                tip = u[cursor + 6 : cursor + 12]
                ja_arm = core["actuator_jacobian"](q, p)
                ja = ca.horzcat(ca.MX.zeros(actuation.count, nbase), ja_arm)
                gamma = core["actuator_velocity_bias"](q, dq, p)
                force = applied(q, external, vehicle, tip, p)
                return command, ja, gamma, force

            def strict_solution(q, dq, u, p):
                command, ja, gamma, force = strict_parts(q, dq, u, p)
                mass = core["mass"](q, p)
                bias = core["bias"](q, dq, p)
                kkt = ca.vertcat(
                    ca.horzcat(mass, ja.T),
                    ca.horzcat(ja, ca.MX.zeros(actuation.count, actuation.count)),
                )
                return ca.solve(kkt, ca.vertcat(force - bias, command - gamma))

            def strict_explicit(q, dq, u, p):
                return strict_solution(q, dq, u, p)[:nq]

            def strict_implicit(q, dq, ddq, u, z, p):
                command, ja, gamma, force = strict_parts(q, dq, u, p)
                return ca.vertcat(
                    mechanical_terms(q, dq, ddq, p) - force + ja.T @ z,
                    ja @ ddq + gamma - command,
                )

            add_variant(
                "strict_tendon_acceleration",
                strict_blocks,
                strict_explicit,
                strict_implicit,
                [("tension", actuation.count)],
            )
            q = ca.MX.sym("q", nq)
            dq = ca.MX.sym("dq", nq)
            u = ca.MX.sym("u", sum(size for _, size in strict_blocks))
            p = ca.MX.sym("p", parameter_count)
            solution = strict_solution(q, dq, u, p)
            functions["dynamics.strict_tendon_acceleration.solution"] = ca.Function(
                "softarm_strict_tendon_acceleration_solution",
                [q, dq, u, p],
                [solution[:nq], solution[nq:]],
                ["q", "dq", "u", "p"],
                ["ddq", "tension"],
            )
            variants["strict_tendon_acceleration"]["solution"] = (
                "dynamics.strict_tendon_acceleration.solution"
            )

    if constraint is not None:
        constraint_blocks = [
            ("tau_arm", narm),
            ("w_vehicle", 6),
            ("w_tip", 6),
            ("constraint_acceleration", constraint.count),
        ]

        def constraint_parts(q, dq, u, p):
            cursor = 0
            tau = u[:narm]
            cursor += narm
            vehicle = u[cursor : cursor + 6]
            cursor += 6
            tip = u[cursor : cursor + 6]
            command = u[cursor + 6 :]
            force = applied(q, tau, vehicle, tip, p)
            phi = core["constraint_value"](q, p)
            jacobian = core["constraint_jacobian"](q, p)
            gamma = core["constraint_velocity_bias"](q, dq, p)
            reaction_map = core["constraint_reaction_map"](q, dq, p)
            gains = core["constraint_stabilization"](p)
            omega, zeta = gains[:, 0], gains[:, 1]
            target = (
                command
                - gamma
                - 2 * zeta * omega * (jacobian @ dq)
                - omega**2 * phi
            )
            return force, jacobian, reaction_map, target

        def constraint_solution(q, dq, u, p):
            force, jacobian, reaction_map, target = constraint_parts(q, dq, u, p)
            mass = core["mass"](q, p)
            bias = core["bias"](q, dq, p)
            kkt = ca.vertcat(
                ca.horzcat(mass, -reaction_map),
                ca.horzcat(
                    jacobian, ca.MX.zeros(constraint.count, constraint.count)
                ),
            )
            return ca.solve(kkt, ca.vertcat(force - bias, target))

        def constraint_explicit(q, dq, u, p):
            return constraint_solution(q, dq, u, p)[:nq]

        def constraint_implicit(q, dq, ddq, u, z, p):
            force, jacobian, reaction_map, target = constraint_parts(q, dq, u, p)
            return ca.vertcat(
                mechanical_terms(q, dq, ddq, p) - force - reaction_map @ z,
                jacobian @ ddq - target,
            )

        add_variant(
            "active_constraint",
            constraint_blocks,
            constraint_explicit,
            constraint_implicit,
            [("reaction", constraint.count)],
        )
        q = ca.MX.sym("q", nq)
        dq = ca.MX.sym("dq", nq)
        u = ca.MX.sym("u", sum(size for _, size in constraint_blocks))
        p = ca.MX.sym("p", parameter_count)
        solution = constraint_solution(q, dq, u, p)
        functions["dynamics.active_constraint.solution"] = ca.Function(
            "softarm_active_constraint_solution",
            [q, dq, u, p],
            [solution[:nq], solution[nq:]],
            ["q", "dq", "u", "p"],
            ["ddq", "reaction"],
        )
        variants["active_constraint"]["solution"] = "dynamics.active_constraint.solution"

    return functions, variants


def _write_model_document(
    plant: PlantModel,
    target: Path,
    actuation: ActuationModel | None,
    constraint: ConstraintModel | None,
) -> None:
    if isinstance(plant, SymbolicPlant):
        from .latex import generate_latex_document

        generate_latex_document(
            plant,
            target,
            actuation=actuation,
            constraint=constraint,
            include_appendix=False,
        )
        return
    (target / "softarm_model.tex").write_text(
        "\\documentclass{article}\n\\usepackage{amsmath}\n\\begin{document}\n"
        "\\section*{Recursive Soft-Arm Dynamics}\n"
        "The CasADi model evaluates recursive inverse dynamics directly in the "
        "implicit residual. Explicit forward dynamics obtains the mass matrix "
        "by algorithmic differentiation with respect to generalized acceleration "
        "and solves the resulting linear system without forming an inverse.\n"
        "\\end{document}\n",
        encoding="utf-8",
    )


def generate_casadi_bundle(
    plant: PlantModel,
    output: str | Path,
    actuation: ActuationModel | None = None,
    constraint: ConstraintModel | None = None,
    *,
    affine_terms: int | None = None,
) -> Path:
    """Generate serialized CasADi Functions from the canonical SymPy model."""
    ca = _import_casadi()
    if affine_terms is not None and not 1 <= affine_terms <= 256:
        raise ValueError("affine_terms must lie in [1, 256]")
    if plant.config.parameterization == "pac" and affine_terms is None:
        raise ValueError("PAC CasADi export requires --casadi-affine-terms")
    if plant.config.parameterization != "pac" and affine_terms is not None:
        raise ValueError("--casadi-affine-terms is only valid for PAC models")
    if isinstance(plant, RecursivePlant) and constraint is not None:
        raise ValueError("recursive dynamics do not support constraints")

    target = Path(output).resolve()
    function_dir = target / "functions"
    function_dir.mkdir(parents=True, exist_ok=True)
    for stale in function_dir.glob("*.casadi"):
        stale.unlink()
    parameters = _runtime_parameters(plant, actuation, constraint)
    if isinstance(plant, SymbolicPlant):
        core = _symbolic_core(
            plant, actuation, constraint, parameters, affine_terms
        )
        dynamics_formulation = "symbolic_lagrange"
    elif isinstance(plant, RecursivePlant):
        core = _recursive_core(plant, actuation, parameters, affine_terms)
        dynamics_formulation = "recursive"
    else:
        raise TypeError("unsupported plant representation")
    if actuation is not None:
        q = ca.MX.sym("q", len(plant.q))
        tension = ca.MX.sym("tension", actuation.count)
        p = ca.MX.sym("p", len(parameters))
        core["actuator_force"] = ca.Function(
            "softarm_actuator_force",
            [q, tension, p],
            [-core["actuator_jacobian"](q, p).T @ tension],
            ["q", "tension", "p"],
            ["tau_arm"],
        )
    core["applied_force"] = _applied_force_function(plant, core, len(parameters))
    dynamics, variants = _dynamics_functions(
        plant, core, actuation, constraint, len(parameters)
    )

    registry = _FunctionRegistry(target, {})
    for name, function in core.items():
        registry.add(f"core.{name}", function)
    for name, function in dynamics.items():
        registry.add(name, function)

    actuation_manifest = None if actuation is None else {
        "family": actuation.family,
        "acceleration": actuation.acceleration,
        "channels": [
            {"name": name, "kind": kind}
            for name, kind in zip(
                actuation.channel_names, actuation.channel_kinds, strict=True
            )
        ],
    }
    constraint_manifest = None if constraint is None else {
        "family": constraint.family,
        "channels": [
            {"name": name, "kind": kind}
            for name, kind in zip(
                constraint.channel_names, constraint.channel_kinds, strict=True
            )
        ],
    }
    manifest = {
        "target": {
            "name": "casadi",
            "schema_version": 1,
            "casadi_version": ca.__version__,
            "symbolic_strategy": "SX kernels with MX composition",
            "affine_terms": affine_terms,
        },
        "model": {
            "rod": plant.config.rod,
            "parameterization": plant.config.parameterization,
            "segments": plant.config.segments,
            "base_mode": plant.config.base.mode,
            "mount_xyz": plant.config.base.mount_xyz,
            "mount_rpy": plant.config.base.mount_rpy,
            "centerline": "centerline_at" in core,
            "dynamics_formulation": dynamics_formulation,
        },
        "coordinates": {
            "base": plant.base_coordinate_names,
            "arm": plant.arm_coordinate_names,
        },
        "parameters": [
            {"name": item.name, "default": item.default} for item in parameters
        ],
        "actuation": actuation_manifest,
        "constraint": constraint_manifest,
        "functions": registry.functions,
        "dynamics": variants,
    }
    (target / "manifest.json").write_text(
        json.dumps(manifest, indent=2), encoding="utf-8"
    )
    _write_model_document(plant, target, actuation, constraint)
    return target
