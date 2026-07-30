from __future__ import annotations

from pathlib import Path

from ..actuation import ActuationModel
from ..derive import SymbolicPlant
from .matlab import render_function, symbol_loads


_OPTIONAL_FILES = (
    "softarm_actuator_coordinates.m",
    "softarm_actuator_jacobian.m",
    "softarm_actuator_velocity_bias.m",
    "softarm_actuator_force.m",
    "softarm_actuator_info.m",
    "softarm_actuator_acceleration.m",
)


def clear_actuator_functions(output: str | Path) -> None:
    target = Path(output).resolve()
    for filename in _OPTIONAL_FILES:
        (target / filename).unlink(missing_ok=True)


def generate_actuator_matlab(
    plant: SymbolicPlant,
    actuation: ActuationModel,
    output: str | Path,
    backend: str = "sympy",
    wolfram_kernel: str | None = None,
) -> Path:
    target = Path(output).resolve()
    target.mkdir(parents=True, exist_ok=True)
    q_loads = symbol_loads(plant.q, "q")
    dq_loads = symbol_loads(plant.dq, "dq")
    all_parameters = plant.parameters + actuation.parameters
    p_loads = symbol_loads([item.symbol for item in all_parameters], "p")
    common = q_loads + p_loads

    render_function(
        target / "softarm_actuator_coordinates.m",
        "softarm_actuator_coordinates",
        "y",
        actuation.coordinates,
        actuation.coordinates.shape,
        ["q", "p"],
        common,
        backend,
        wolfram_kernel,
    )
    render_function(
        target / "softarm_actuator_jacobian.m",
        "softarm_actuator_jacobian",
        "Ja",
        actuation.jacobian,
        actuation.jacobian.shape,
        ["q", "p"],
        common,
        backend,
        wolfram_kernel,
    )
    render_function(
        target / "softarm_actuator_velocity_bias.m",
        "softarm_actuator_velocity_bias",
        "jdotdq",
        actuation.velocity_bias,
        actuation.velocity_bias.shape,
        ["q", "dq", "p"],
        q_loads + dq_loads + p_loads,
        backend,
        wolfram_kernel,
    )

    unilateral = [
        str(index) for index, kind in enumerate(actuation.channel_kinds, start=1)
        if kind == "unilateral"
    ]
    index_vector = " ".join(unilateral)
    feasibility = (
        f"isFeasible=all(tension([{index_vector}])>=0);\n"
        if unilateral else "isFeasible=true;\n"
    )
    (target / "softarm_actuator_force.m").write_text(
        (
            "function [tau,isFeasible] = softarm_actuator_force(q,tension,p)\n"
            "%SOFTARM_ACTUATOR_FORCE Map ordered channel tensions to generalized force.\n"
            "%#codegen\n"
            f"assert(numel(tension)=={actuation.count});\n"
            "tension=tension(:);\n"
            f"{feasibility}"
            "assert(isFeasible,'softarm:NegativeUnilateralTension',"
            "'Unilateral tendon tension must be nonnegative.');\n"
            "Ja=softarm_actuator_jacobian(q,p);\n"
            "tau=-Ja.'*tension;\n"
            "end\n"
        ),
        encoding="utf-8",
    )

    if actuation.acceleration == "strict":
        acceleration_feasibility = (
            f"isPullOnlyFeasible=all(tension([{index_vector}])>=-sqrt(eps));\n"
            if unilateral else "isPullOnlyFeasible=true;\n"
        )
        (target / "softarm_actuator_acceleration.m").write_text(
            (
                "function [ddq,tension,isPullOnlyFeasible,reciprocalCondition] = "
                "softarm_actuator_acceleration(q,dq,coordinateAcceleration,tauArmExternal,wVehicle,wTip,p)\n"
                "%SOFTARM_ACTUATOR_ACCELERATION Enforce independent actuator-coordinate accelerations.\n"
                "%#codegen\n"
                f"assert(numel(coordinateAcceleration)=={actuation.count});\n"
                "coordinateAcceleration=coordinateAcceleration(:);\n"
                "M=softarm_mass(q,p); h=softarm_bias(q,dq,p);\n"
                "Ja=softarm_actuator_jacobian(q,p);\n"
                f"JaFull=[zeros({actuation.count},{len(plant.base_q)}),Ja];\n"
                "jdotdq=softarm_actuator_velocity_bias(q,dq,p);\n"
                "S=JaFull*(M\\JaFull.'); reciprocalCondition=rcond(S);\n"
                "assert(isfinite(reciprocalCondition) && reciprocalCondition>sqrt(eps),"
                "'softarm:SingularActuatorConstraint','Actuator acceleration constraints are singular.');\n"
                "Q=softarm_applied_force(q,tauArmExternal,wVehicle,wTip,p);\n"
                "n=numel(q); m=size(JaFull,1); solution=[M,JaFull.';JaFull,zeros(m)]\\"
                "[Q-h;coordinateAcceleration-jdotdq];\n"
                "ddq=solution(1:n); tension=solution(n+1:n+m);\n"
                f"{acceleration_feasibility}"
                "end\n"
            ),
            encoding="utf-8",
        )
    return target
