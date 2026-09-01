from __future__ import annotations

from pathlib import Path

import sympy as sp

from ..constraints import (
    ConstraintDefinition,
    ConstraintModel,
    point_constraint_expressions,
)
from ..models import RecursivePlant, RuntimeParameter, SymbolicPlant
from .matlab import render_function, symbol_loads
from .optimization import FunctionOptimizer

_OPTIONAL_FILES = (
    "softarm_constraint_value.m",
    "softarm_constraint_jacobian.m",
    "softarm_constraint_velocity_bias.m",
    "softarm_constraint_reaction_map.m",
    "softarm_constraint_stabilization.m",
    "softarm_constraint_acceleration.m",
    "softarm_constraint_terms.m",
    "softarm_point_constraint_kernel.m",
)


def clear_constraint_functions(output: str | Path) -> None:
    target = Path(output).resolve()
    for filename in _OPTIONAL_FILES:
        (target / filename).unlink(missing_ok=True)


def generate_constraint_matlab(
    plant: SymbolicPlant,
    constraint: ConstraintModel,
    output: str | Path,
    preceding_parameters: tuple[RuntimeParameter, ...] = (),
    optimizer: FunctionOptimizer | None = None,
) -> Path:
    function_optimizer = optimizer or FunctionOptimizer()
    target = Path(output).resolve()
    target.mkdir(parents=True, exist_ok=True)
    q_loads = symbol_loads(plant.q, "q")
    dq_loads = symbol_loads(plant.dq, "dq")
    all_parameters = plant.parameters + preceding_parameters + constraint.parameters
    p_loads = symbol_loads([item.symbol for item in all_parameters], "p")
    common = q_loads + p_loads

    for filename, function, output_name, matrix, inputs, loads in (
        ("softarm_constraint_value.m", "softarm_constraint_value", "phi", constraint.coordinates, ["q", "p"], common),
        ("softarm_constraint_jacobian.m", "softarm_constraint_jacobian", "A", constraint.jacobian, ["q", "p"], common),
        ("softarm_constraint_velocity_bias.m", "softarm_constraint_velocity_bias", "gamma", constraint.velocity_bias, ["q", "dq", "p"], q_loads + dq_loads + p_loads),
        ("softarm_constraint_reaction_map.m", "softarm_constraint_reaction_map", "G", constraint.reaction_map, ["q", "dq", "p"], q_loads + dq_loads + p_loads),
    ):
        render_function(
            target / filename, function, output_name, matrix, matrix.shape,
            inputs, loads, function_optimizer,
        )

    stabilization = constraint.stabilization_frequency.row_join(
        constraint.stabilization_ratio
    )
    render_function(
        target / "softarm_constraint_stabilization.m",
        "softarm_constraint_stabilization",
        "gains",
        stabilization,
        stabilization.shape,
        ["p"],
        p_loads,
        function_optimizer,
    )

    unilateral = [
        str(index) for index, kind in enumerate(constraint.channel_kinds, start=1)
        if kind == "unilateral"
    ]
    feasibility = (
        f"isFeasible=all(reaction([{ ' '.join(unilateral) }])>=-sqrt(eps));\n"
        if unilateral else "isFeasible=true;\n"
    )
    normal_indices = [
        index for index, parameter in enumerate(all_parameters, start=1)
        if parameter.name.startswith("constraint_plane_normal_")
    ]
    normal_check = ""
    if normal_indices:
        vector = " ".join(str(index) for index in normal_indices)
        normal_check = (
            f"assert(norm(p([{vector}]))>sqrt(eps),'softarm:InvalidPlaneNormal',"
            "'Plane normal must be nonzero.');\n"
        )
    (target / "softarm_constraint_acceleration.m").write_text(
        (
            "function [ddq,reaction,isFeasible,reciprocalCondition] = "
            "softarm_constraint_acceleration(q,dq,tauArm,wVehicle,wTip,constraintAcceleration,p)\n"
            "%SOFTARM_CONSTRAINT_ACCELERATION Solve the active acceleration constraint.\n"
            "%#codegen\n"
            f"assert(numel(constraintAcceleration)=={constraint.count});\n"
            f"{normal_check}"
            "constraintAcceleration=constraintAcceleration(:);\n"
            "M=softarm_mass(q,p); h=softarm_bias(q,dq,p);\n"
            "Q=softarm_applied_force(q,tauArm,wVehicle,wTip,p);\n"
            "phi=softarm_constraint_value(q,p); A=softarm_constraint_jacobian(q,p);\n"
            "gamma=softarm_constraint_velocity_bias(q,dq,p);\n"
            "G=softarm_constraint_reaction_map(q,dq,p);\n"
            "gains=softarm_constraint_stabilization(p); omega=gains(:,1); zeta=gains(:,2);\n"
            "target=constraintAcceleration-gamma-2*zeta.*omega.*(A*dq)-(omega.^2).*phi;\n"
            "n=numel(q); m=size(A,1); K=[M,-G;A,zeros(m)]; reciprocalCondition=rcond(K);\n"
            "assert(isfinite(reciprocalCondition) && reciprocalCondition>eps,"
            "'softarm:SingularConstraint','Constraint KKT system is singular or ill-conditioned.');\n"
            "solution=K\\[Q-h;target]; ddq=solution(1:n); reaction=solution(n+1:n+m);\n"
            f"{feasibility}"
            "end\n"
        ),
        encoding="utf-8",
    )
    return target


def generate_recursive_constraint_matlab(
    plant: RecursivePlant,
    constraint: ConstraintDefinition,
    output: str | Path,
    preceding_parameters: tuple[RuntimeParameter, ...] = (),
    optimizer: FunctionOptimizer | None = None,
) -> Path:
    """Generate the public constraint API from a recursive tool-point jet."""
    function_optimizer = optimizer or FunctionOptimizer()
    target = Path(output).resolve()
    target.mkdir(parents=True, exist_ok=True)
    nq = len(plant.q)
    position = sp.Matrix(sp.symbols("point_px point_py point_pz", real=True))
    jacobian_symbols = sp.Matrix(
        sp.symbols(f"point_j1:{3 * nq + 1}", real=True)
    )
    point_jacobian = sp.Matrix(3, nq, list(jacobian_symbols))
    point_bias = sp.Matrix(sp.symbols("point_gx point_gy point_gz", real=True))
    point_velocity = sp.Matrix(sp.symbols("point_vx point_vy point_vz", real=True))
    phi, jacobian, gamma, reaction_map = point_constraint_expressions(
        constraint, position, point_jacobian, point_bias, point_velocity
    )
    packed = phi.col_join(jacobian.T).col_join(gamma).col_join(reaction_map)
    all_parameters = plant.parameters + preceding_parameters + constraint.parameters
    loads = (
        symbol_loads(position, "pointPosition")
        + [
            f"{point_jacobian[row, column]} = pointJacobian({row + 1},{column + 1});"
            for row in range(3)
            for column in range(nq)
        ]
        + symbol_loads(point_bias, "pointVelocityBias")
        + symbol_loads(point_velocity, "pointVelocity")
        + symbol_loads([item.symbol for item in all_parameters], "p")
    )
    render_function(
        target / "softarm_point_constraint_kernel.m",
        "softarm_point_constraint_kernel",
        "y",
        packed,
        packed.shape,
        [
            "pointPosition",
            "pointJacobian",
            "pointVelocityBias",
            "pointVelocity",
            "p",
        ],
        loads,
        function_optimizer,
    )
    offset_indices = [
        next(
            index
            for index, item in enumerate(all_parameters, start=1)
            if item.symbol == symbol
        )
        for symbol in constraint.tool_offset
    ]
    offsets = " ".join(str(index) for index in offset_indices)
    packed_size = 2 * nq + 2
    (target / "softarm_constraint_terms.m").write_text(
        (
            "function [phi,A,gamma,G]=softarm_constraint_terms(q,dq,p)\n"
            "%#codegen\n"
            f"offset=p([{offsets}]);\n"
            "[pointPosition,pointJacobian,pointVelocityBias,pointVelocity]="
            "softarm_tool_point_jet(q,dq,offset,p);\n"
            "y=softarm_point_constraint_kernel(pointPosition,pointJacobian,"
            "pointVelocityBias,pointVelocity,p);\n"
            "cursor=1; phi=y(cursor); cursor=cursor+1;\n"
            f"A=reshape(y(cursor:cursor+{nq - 1}),1,{nq}); cursor=cursor+{nq};\n"
            "gamma=y(cursor); cursor=cursor+1;\n"
            f"G=reshape(y(cursor:{packed_size}),{nq},1);\n"
            "end\n"
        ),
        encoding="utf-8",
    )
    wrappers = {
        "softarm_constraint_value.m": (
            "function phi=softarm_constraint_value(q,p)\n"
            f"[phi,~,~,~]=softarm_constraint_terms(q,zeros({nq},1),p);\nend\n"
        ),
        "softarm_constraint_jacobian.m": (
            "function A=softarm_constraint_jacobian(q,p)\n"
            f"[~,A,~,~]=softarm_constraint_terms(q,zeros({nq},1),p);\nend\n"
        ),
        "softarm_constraint_velocity_bias.m": (
            "function gamma=softarm_constraint_velocity_bias(q,dq,p)\n"
            "[~,~,gamma,~]=softarm_constraint_terms(q,dq,p);\nend\n"
        ),
        "softarm_constraint_reaction_map.m": (
            "function G=softarm_constraint_reaction_map(q,dq,p)\n"
            "[~,~,~,G]=softarm_constraint_terms(q,dq,p);\nend\n"
        ),
    }
    for filename, content in wrappers.items():
        (target / filename).write_text(content, encoding="utf-8")
    p_loads = symbol_loads([item.symbol for item in all_parameters], "p")
    gains = constraint.stabilization_frequency.row_join(
        constraint.stabilization_ratio
    )
    render_function(
        target / "softarm_constraint_stabilization.m",
        "softarm_constraint_stabilization",
        "gains",
        gains,
        gains.shape,
        ["p"],
        p_loads,
        function_optimizer,
    )
    _write_constraint_solver(
        plant, constraint, target, all_parameters
    )
    return target


def _write_constraint_solver(
    plant: RecursivePlant,
    constraint: ConstraintDefinition,
    target: Path,
    all_parameters: tuple[RuntimeParameter, ...],
) -> None:
    unilateral = [
        str(index)
        for index, kind in enumerate(constraint.channel_kinds, start=1)
        if kind == "unilateral"
    ]
    feasibility = (
        f"isFeasible=all(reaction([{' '.join(unilateral)}])>=-sqrt(eps));\n"
        if unilateral
        else "isFeasible=true;\n"
    )
    normal_indices = [
        index
        for index, parameter in enumerate(all_parameters, start=1)
        if parameter.name.startswith("constraint_plane_normal_")
    ]
    normal_check = ""
    if normal_indices:
        vector = " ".join(str(index) for index in normal_indices)
        normal_check = (
            f"assert(norm(p([{vector}]))>sqrt(eps),'softarm:InvalidPlaneNormal',"
            "'Plane normal must be nonzero.');\n"
        )
    (target / "softarm_constraint_acceleration.m").write_text(
        (
            "function [ddq,reaction,isFeasible,reciprocalCondition] = "
            "softarm_constraint_acceleration(q,dq,tauArm,wVehicle,wTip,"
            "constraintAcceleration,p)\n%#codegen\n"
            f"assert(numel(constraintAcceleration)=={constraint.count});\n"
            f"{normal_check}constraintAcceleration=constraintAcceleration(:);\n"
            "M=softarm_mass(q,p); h=softarm_bias(q,dq,p);\n"
            "Q=softarm_applied_force(q,tauArm,wVehicle,wTip,p);\n"
            "[phi,A,gamma,G]=softarm_constraint_terms(q,dq,p);\n"
            "gains=softarm_constraint_stabilization(p); omega=gains(:,1); "
            "zeta=gains(:,2);\n"
            "targetAcceleration=constraintAcceleration-gamma-2*zeta.*omega.*"
            "(A*dq)-(omega.^2).*phi;\n"
            "n=numel(q); m=size(A,1); K=[M,-G;A,zeros(m)]; "
            "reciprocalCondition=rcond(K);\n"
            "assert(isfinite(reciprocalCondition)&&reciprocalCondition>eps,"
            "'softarm:SingularConstraint','Constraint KKT system is singular or ill-conditioned.');\n"
            "solution=K\\[Q-h;targetAcceleration]; ddq=solution(1:n); "
            "reaction=solution(n+1:n+m);\n"
            f"{feasibility}end\n"
        ),
        encoding="utf-8",
    )
