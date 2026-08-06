from __future__ import annotations

from pathlib import Path

from ..constraints import ConstraintModel
from ..models import RuntimeParameter, SymbolicPlant
from .matlab import render_function, symbol_loads
from .optimization import FunctionOptimizer

_OPTIONAL_FILES = (
    "softarm_constraint_value.m",
    "softarm_constraint_jacobian.m",
    "softarm_constraint_velocity_bias.m",
    "softarm_constraint_reaction_map.m",
    "softarm_constraint_stabilization.m",
    "softarm_constraint_acceleration.m",
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
