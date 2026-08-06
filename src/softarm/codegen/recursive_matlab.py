from __future__ import annotations

import json
from pathlib import Path

import sympy as sp

from ..actuation import ActuationModel
from ..geometry import angular_jacobian, transform_rpy
from ..models import RecursivePlant
from ..recursive import local_kernel_for
from .matlab import _HELPERS, _flatten, render_function, symbol_loads
from .optimization import FunctionOptimizer


def _parameter(plant: RecursivePlant, name: str) -> sp.Symbol:
    return next(item.symbol for item in plant.parameters if item.name == name)


def _directional(matrix: sp.Matrix, q: sp.Matrix, dq: sp.Matrix) -> sp.Matrix:
    result = sp.zeros(*matrix.shape)
    for coordinate, velocity in zip(q, dq, strict=True):
        result += matrix.diff(coordinate) * velocity
    return result


def _pack(*matrices: sp.Matrix) -> sp.Matrix:
    return sp.Matrix([value for matrix in matrices for value in _flatten(matrix)])


def _state_loads(symbols: sp.Matrix, source: str) -> list[str]:
    return symbol_loads(list(symbols), source)


def _generate_base_kernel(
    plant: RecursivePlant,
    target: Path,
    optimizer: FunctionOptimizer,
) -> None:
    config = plant.config
    base_q = plant.base_q
    base_dq = plant.base_dq
    base_ddq = sp.Matrix(
        sp.symbols(f"base_ddq1:{len(base_q) + 1}", real=True)
    ) if len(base_q) else sp.zeros(0, 1)
    gravity_world = sp.Matrix([0, 0, _parameter(plant, "gravity")])
    mount = transform_rpy(
        tuple(sp.Float(str(value)) for value in config.base.mount_xyz),
        tuple(sp.Float(str(value)) for value in config.base.mount_rpy),
    )
    mount_transform = plant.base_transform * mount
    mount_rotation = mount_transform[:3, :3]
    mount_position = mount_transform[:3, 3]

    if len(base_q):
        mount_linear = mount_position.jacobian(base_q)
        mount_angular = angular_jacobian(mount_rotation, base_q)
        mount_linear_rate = _directional(mount_linear, base_q, base_dq)
        mount_angular_rate = _directional(mount_angular, base_q, base_dq)
        mount_v_world = mount_linear * base_dq
        mount_w_world = mount_angular * base_dq
        mount_a_world = mount_linear * base_ddq + mount_linear_rate * base_dq
        mount_alpha_world = mount_angular * base_ddq + mount_angular_rate * base_dq
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
        base_linear = base_position.jacobian(base_q)
        base_angular = angular_jacobian(base_rotation, base_q)
        base_linear_rate = _directional(base_linear, base_q, base_dq)
        base_angular_rate = _directional(base_angular, base_q, base_dq)
        base_w = base_angular * base_dq
        base_a = base_linear * base_ddq + base_linear_rate * base_dq
        base_alpha = base_angular * base_ddq + base_angular_rate * base_dq
        vehicle_mass = _parameter(plant, "vehicle_mass")
        vehicle_inertia = sp.diag(
            _parameter(plant, "vehicle_Ixx"),
            _parameter(plant, "vehicle_Iyy"),
            _parameter(plant, "vehicle_Izz"),
        )
        inertia_world = base_rotation * vehicle_inertia * base_rotation.T
        vehicle_force = vehicle_mass * (base_a - gravity_world)
        vehicle_torque = inertia_world * base_alpha + base_w.cross(
            inertia_world * base_w
        )
        base_tau = base_linear.T * vehicle_force + base_angular.T * vehicle_torque
        vehicle_map = base_linear.col_join(base_angular).T * sp.diag(
            base_rotation, base_rotation
        )
    else:
        parent_v = sp.zeros(3, 1)
        parent_w = sp.zeros(3, 1)
        parent_a = sp.zeros(3, 1)
        parent_alpha = sp.zeros(3, 1)
        local_gravity = mount_rotation.T * gravity_world
        mount_map = sp.zeros(0, 6)
        base_tau = sp.zeros(0, 1)
        vehicle_map = sp.zeros(0, 6)

    q_loads = symbol_loads(plant.q, "q")
    dq_loads = symbol_loads(plant.dq, "dq")
    p_loads = symbol_loads([item.symbol for item in plant.parameters], "p")
    ddq_loads = _state_loads(base_ddq, "ddq")
    packed = _pack(
        parent_v,
        parent_w,
        parent_a,
        parent_alpha,
        local_gravity,
        base_tau,
        mount_map,
    )
    render_function(
        target / "softarm_recursive_base.m",
        "softarm_recursive_base",
        "y",
        packed,
        packed.shape,
        ["q", "dq", "ddq", "p"],
        q_loads + dq_loads + ddq_loads + p_loads,
        optimizer,
    )
    render_function(
        target / "softarm_mount_transform.m",
        "softarm_mount_transform",
        "H",
        mount_transform,
        (4, 4),
        ["q", "p"],
        q_loads + p_loads,
        optimizer,
    )
    mount_jacobian = (
        mount_position.jacobian(plant.q).col_join(
            angular_jacobian(mount_rotation, plant.q)
        )
        if len(plant.q)
        else sp.zeros(6, 0)
    )
    render_function(
        target / "softarm_mount_jacobian.m",
        "softarm_mount_jacobian",
        "J",
        mount_jacobian,
        mount_jacobian.shape,
        ["q", "p"],
        q_loads + p_loads,
        optimizer,
    )
    full_vehicle_map = sp.zeros(len(plant.q), 6)
    if len(base_q):
        full_vehicle_map[: len(base_q), :] = vehicle_map
    render_function(
        target / "softarm_vehicle_wrench_map.m",
        "softarm_vehicle_wrench_map",
        "Bv",
        full_vehicle_map,
        full_vehicle_map.shape,
        ["q", "p"],
        q_loads + p_loads,
        optimizer,
    )


def _generate_section_functions(
    plant: RecursivePlant,
    target: Path,
    optimizer: FunctionOptimizer,
) -> int:
    definition = plant.definition
    local_kernel = local_kernel_for(definition)
    dof = len(plant.arm_q) // plant.config.segments
    base_count = len(plant.base_q)
    own_cases: list[str] = []
    end_cases: list[str] = []
    pose_cases: list[str] = []
    xi = sp.Symbol("xi", real=True, nonnegative=True)
    kernel = local_kernel.build(definition, 0)
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
        zip(
            (item.symbol for item in section_parameters),
            generic_p,
            strict=True,
        )
    )

    own = _pack(kernel.own_wrench, kernel.own_generalized_force).xreplace(
        substitutions
    )
    own_loads = (
        _state_loads(generic_q, "qLocal")
        + _state_loads(generic_dq, "dqLocal")
        + _state_loads(kernel.acceleration_symbols, "ddqSection")
        + _state_loads(generic_p, "pLocal")
        + _state_loads(kernel.parent_velocity_symbols[:3, 0], "v0")
        + _state_loads(kernel.parent_velocity_symbols[3:, 0], "w0")
        + _state_loads(kernel.parent_acceleration_symbols[:3, 0], "a0")
        + _state_loads(kernel.parent_acceleration_symbols[3:, 0], "alpha0")
        + _state_loads(kernel.gravity_symbols, "g0")
    )
    render_function(
        target / "softarm_recursive_section_template.m",
        "softarm_recursive_section_template",
        "y",
        own,
        own.shape,
        [
            "qLocal",
            "dqLocal",
            "ddqSection",
            "pLocal",
            "v0",
            "w0",
            "a0",
            "alpha0",
            "g0",
        ],
        own_loads,
        optimizer,
    )

    end = _pack(
        kernel.end_rotation,
        kernel.end_position,
        kernel.end_linear_jacobian,
        kernel.end_angular_jacobian,
        kernel.end_linear_jacobian_rate,
        kernel.end_angular_jacobian_rate,
    ).xreplace(substitutions)
    render_function(
        target / "softarm_recursive_end_template.m",
        "softarm_recursive_end_template",
        "y",
        end,
        end.shape,
        ["qLocal", "dqLocal", "pLocal"],
        _state_loads(generic_q, "qLocal")
        + _state_loads(generic_dq, "dqLocal")
        + _state_loads(generic_p, "pLocal"),
        optimizer,
    )

    pose = definition.kinematics.transform(0, xi).xreplace(substitutions)
    render_function(
        target / "softarm_recursive_pose_template.m",
        "softarm_recursive_pose_template",
        "H",
        pose,
        (4, 4),
        ["qLocal", "pLocal", "xi"],
        _state_loads(generic_q, "qLocal") + _state_loads(generic_p, "pLocal"),
        optimizer,
    )

    parameter_indices = {
        item.name: index
        for index, item in enumerate(plant.parameters, start=1)
    }
    suffixes = [item.name.removeprefix("s1_") for item in section_parameters]
    for section in range(plant.config.segments):
        first = base_count + section * dof + 1
        last = first + dof - 1
        indices = ",".join(
            str(parameter_indices[f"s{section + 1}_{suffix}"])
            for suffix in suffixes
        )
        local_arguments = f"q({first}:{last}),dq({first}:{last}),p([{indices}])"
        own_cases.append(
            f"    case {section + 1}, y=softarm_recursive_section_template("
            f"{local_arguments.split(',p([')[0]},ddqSection,p([{indices}]),"
            "v0,w0,a0,alpha0,g0);"
        )
        end_cases.append(
            f"    case {section + 1}, y=softarm_recursive_end_template({local_arguments});"
        )
        pose_cases.append(
            f"    case {section + 1}, H=softarm_recursive_pose_template("
            f"q({first}:{last}),p([{indices}]),xi);"
        )

    (target / "softarm_recursive_section.m").write_text(
        "\n".join([
            "function y=softarm_recursive_section(section,q,dq,ddqSection,p,v0,w0,a0,alpha0,g0)",
            "%#codegen",
            "switch section",
            *own_cases,
            "    otherwise, error('softarm:InvalidSection','Invalid section index.');",
            "end",
            "end",
            "",
        ]),
        encoding="utf-8",
    )
    (target / "softarm_recursive_end.m").write_text(
        "\n".join([
            "function y=softarm_recursive_end(section,q,dq,p)",
            "%#codegen",
            "switch section",
            *end_cases,
            "    otherwise, error('softarm:InvalidSection','Invalid section index.');",
            "end",
            "end",
            "",
        ]),
        encoding="utf-8",
    )
    (target / "softarm_recursive_pose.m").write_text(
        "\n".join([
            "function H=softarm_recursive_pose(section,q,p,xi)",
            "%#codegen",
            "switch section",
            *pose_cases,
            "    otherwise, error('softarm:InvalidSection','Invalid section index.');",
            "end",
            "end",
            "",
        ]),
        encoding="utf-8",
    )
    return dof


def _generate_internal_and_tip(
    plant: RecursivePlant,
    target: Path,
    optimizer: FunctionOptimizer,
) -> None:
    definition = plant.definition
    arm_internal = sp.Matrix(
        [sp.diff(definition.elastic, coordinate) for coordinate in plant.arm_q]
    )
    damping = (
        sp.diag(*definition.damping)
        if isinstance(definition.damping, tuple)
        else definition.damping
    )
    arm_internal += damping * plant.arm_dq
    internal = sp.zeros(len(plant.q), 1)
    internal[len(plant.base_q) :, 0] = arm_internal
    render_function(
        target / "softarm_recursive_internal_force.m",
        "softarm_recursive_internal_force",
        "tau",
        internal,
        internal.shape,
        ["q", "dq", "p"],
        symbol_loads(plant.q, "q")
        + symbol_loads(plant.dq, "dq")
        + symbol_loads([item.symbol for item in plant.parameters], "p"),
        optimizer,
    )

    velocity = sp.Matrix(sp.symbols("tip_vx tip_vy tip_vz tip_wx tip_wy tip_wz", real=True))
    acceleration = sp.Matrix(
        sp.symbols("tip_ax tip_ay tip_az tip_alphax tip_alphay tip_alphaz", real=True)
    )
    gravity = sp.Matrix(sp.symbols("tip_gx tip_gy tip_gz", real=True))
    mass = _parameter(plant, "tip_mass")
    inertia = sp.diag(
        _parameter(plant, "tip_Ixx"),
        _parameter(plant, "tip_Iyy"),
        _parameter(plant, "tip_Izz"),
    )
    force = mass * (acceleration[:3, 0] - gravity)
    torque = inertia * acceleration[3:, 0] + velocity[3:, 0].cross(
        inertia * velocity[3:, 0]
    )
    wrench = force.col_join(torque)
    loads = (
        _state_loads(velocity[:3, 0], "v")
        + _state_loads(velocity[3:, 0], "w")
        + _state_loads(acceleration[:3, 0], "a")
        + _state_loads(acceleration[3:, 0], "alpha")
        + _state_loads(gravity, "g")
        + symbol_loads([item.symbol for item in plant.parameters], "p")
    )
    render_function(
        target / "softarm_recursive_tip.m",
        "softarm_recursive_tip",
        "wrench",
        wrench,
        (6, 1),
        ["v", "w", "a", "alpha", "g", "p"],
        loads,
        optimizer,
    )


def _write_runtime(
    plant: RecursivePlant,
    target: Path,
    dof: int,
) -> None:
    segments = plant.config.segments
    base_count = len(plant.base_q)
    nq = len(plant.q)
    base_pack = 15 + base_count + 6 * base_count
    runtime = f"""function tau=softarm_inverse_dynamics(q,dq,ddq,p)
%SOFTARM_INVERSE_DYNAMICS Recursive section inverse dynamics.
%#codegen
assert(numel(q)=={nq}&&numel(dq)=={nq}&&numel(ddq)=={nq});
q=q(:); dq=dq(:); ddq=ddq(:);
base=softarm_recursive_base(q,dq,ddq(1:{base_count}),p);
assert(numel(base)=={base_pack});
v=base(1:3); w=base(4:6); a=base(7:9); alpha=base(10:12); g=base(13:15);
baseTau=base(16:15+{base_count});
mountMap=reshape(base(16+{base_count}:end),{base_count},6);
Rend=zeros(3,3,{segments}); rend=zeros(3,{segments});
Jv=zeros(3,{dof},{segments}); Jw=zeros(3,{dof},{segments});
ownWrench=zeros(6,{segments}); ownTau=zeros({dof},{segments});
for section=1:{segments}
    first={base_count}+(section-1)*{dof}+1; last=first+{dof}-1;
    jet=softarm_recursive_end(section,q,dq,p);
    cursor=1; Rend(:,:,section)=reshape(jet(cursor:cursor+8),3,3); cursor=cursor+9;
    rend(:,section)=jet(cursor:cursor+2); cursor=cursor+3;
    Jv(:,:,section)=reshape(jet(cursor:cursor+{3*dof-1}),3,{dof}); cursor=cursor+{3*dof};
    Jw(:,:,section)=reshape(jet(cursor:cursor+{3*dof-1}),3,{dof}); cursor=cursor+{3*dof};
    Jvd=reshape(jet(cursor:cursor+{3*dof-1}),3,{dof}); cursor=cursor+{3*dof};
    Jwd=reshape(jet(cursor:cursor+{3*dof-1}),3,{dof});
    local=softarm_recursive_section(section,q,dq,ddq(first:last),p,v,w,a,alpha,g);
    ownWrench(:,section)=local(1:6); ownTau(:,section)=local(7:end);
    relativeV=Jv(:,:,section)*dq(first:last); relativeW=Jw(:,:,section)*dq(first:last);
    vParent=v+cross(w,rend(:,section))+relativeV;
    wParent=w+relativeW;
    aParent=a+cross(alpha,rend(:,section))+cross(w,cross(w,rend(:,section)))+2*cross(w,relativeV)+Jv(:,:,section)*ddq(first:last)+Jvd*dq(first:last);
    alphaParent=alpha+cross(w,relativeW)+Jw(:,:,section)*ddq(first:last)+Jwd*dq(first:last);
    rotation=Rend(:,:,section);
    v=rotation.'*vParent; w=rotation.'*wParent;
    a=rotation.'*aParent; alpha=rotation.'*alphaParent; g=rotation.'*g;
end
wrench=softarm_recursive_tip(v,w,a,alpha,g,p);
tauArm=zeros({dof*segments},1);
for section={segments}:-1:1
    rotation=Rend(:,:,section); force=rotation*wrench(1:3); moment=rotation*wrench(4:6);
    localFirst=(section-1)*{dof}+1; localLast=localFirst+{dof}-1;
    tauArm(localFirst:localLast)=ownTau(:,section)+Jv(:,:,section).'*force+Jw(:,:,section).'*moment;
    wrench=ownWrench(:,section)+[force;cross(rend(:,section),force)+moment];
end
tau=[baseTau+mountMap*wrench;tauArm]+softarm_recursive_internal_force(q,dq,p);
end
"""
    (target / "softarm_inverse_dynamics.m").write_text(runtime, encoding="utf-8")
    (target / "softarm_mass.m").write_text(
        f"""function M=softarm_mass(q,p)
%#codegen
n={nq}; z=zeros(n,1); reference=softarm_inverse_dynamics(q,z,z,p); M=zeros(n,n);
for column=1:n
    acceleration=zeros(n,1); acceleration(column)=1;
    M(:,column)=softarm_inverse_dynamics(q,z,acceleration,p)-reference;
end
M=(M+M.')/2;
end
""",
        encoding="utf-8",
    )
    (target / "softarm_bias.m").write_text(
        f"function h=softarm_bias(q,dq,p)\n%#codegen\nh=softarm_inverse_dynamics(q,dq,zeros({nq},1),p);\nend\n",
        encoding="utf-8",
    )
    (target / "softarm_kinematics.m").write_text(
        f"""function H=softarm_kinematics(q,p)
%#codegen
H=zeros(4,4,{segments}); current=softarm_mount_transform(q,p);
for section=1:{segments}
    current=current*softarm_recursive_pose(section,q,p,1); H(:,:,section)=current;
end
end
""",
        encoding="utf-8",
    )
    (target / "softarm_centerline_at.m").write_text(
        f"""function P=softarm_centerline_at(q,p,xi)
%#codegen
P=zeros(3,{segments}); current=softarm_mount_transform(q,p);
for section=1:{segments}
    material=current*softarm_recursive_pose(section,q,p,xi); P(:,section)=material(1:3,4);
    current=current*softarm_recursive_pose(section,q,p,1);
end
end
""",
        encoding="utf-8",
    )
    (target / "softarm_centerline.m").write_text(
        f"""function P=softarm_centerline(q,p,xi)
%#codegen
xi=double(xi(:).'); assert(~isempty(xi)&&all(isfinite(xi))&&all(xi>=0)&&all(xi<=1));
P=zeros(3,{segments}*numel(xi));
for sample=1:numel(xi)
    points=softarm_centerline_at(q,p,xi(sample));
    for section=1:{segments}, P(:,(section-1)*numel(xi)+sample)=points(:,section); end
end
end
""",
        encoding="utf-8",
    )
    (target / "softarm_end_jacobian.m").write_text(
        f"""function J=softarm_end_jacobian(q,p)
%#codegen
J=softarm_mount_jacobian(q,p); T=softarm_mount_transform(q,p); zeroDq=zeros({nq},1);
for section=1:{segments}
    jet=softarm_recursive_end(section,q,zeroDq,p); cursor=1;
    rotation=reshape(jet(cursor:cursor+8),3,3); cursor=cursor+9;
    position=jet(cursor:cursor+2); cursor=cursor+3;
    localJv=reshape(jet(cursor:cursor+{3*dof-1}),3,{dof}); cursor=cursor+{3*dof};
    localJw=reshape(jet(cursor:cursor+{3*dof-1}),3,{dof});
    worldOffset=T(1:3,1:3)*position;
    J(1:3,:)=J(1:3,:)-softarm_skew(worldOffset)*J(4:6,:);
    first={base_count}+(section-1)*{dof}+1; last=first+{dof}-1;
    J(1:3,first:last)=J(1:3,first:last)+T(1:3,1:3)*localJv;
    J(4:6,first:last)=J(4:6,first:last)+T(1:3,1:3)*localJw;
    T=T*[rotation position;0 0 0 1];
end
end
""",
        encoding="utf-8",
    )
    (target / "softarm_skew.m").write_text(
        "function S=softarm_skew(v)\n%#codegen\nS=[0,-v(3),v(2);v(3),0,-v(1);-v(2),v(1),0];\nend\n",
        encoding="utf-8",
    )


def generate_recursive_matlab_bundle(
    plant: RecursivePlant,
    output: str | Path,
    actuation: ActuationModel | None = None,
    optimizer: FunctionOptimizer | None = None,
) -> Path:
    function_optimizer = optimizer or FunctionOptimizer()
    target = Path(output).resolve()
    target.mkdir(parents=True, exist_ok=True)
    if plant.definition.linear_kinematics:
        from .legacy_ritz_matlab import generate_legacy_ritz_core

        generate_legacy_ritz_core(plant, target, function_optimizer)
    else:
        _generate_base_kernel(plant, target, function_optimizer)
        dof = _generate_section_functions(plant, target, function_optimizer)
        _generate_internal_and_tip(plant, target, function_optimizer)
        _write_runtime(plant, target, dof)

    (target / "softarm_applied_force.m").write_text(
        f"""function Q=softarm_applied_force(q,tauArm,wVehicle,wTip,p)
%#codegen
assert(numel(tauArm)=={len(plant.arm_q)}); assert(numel(wVehicle)==6); assert(numel(wTip)==6);
Q=[zeros({len(plant.base_q)},1);tauArm(:)]+softarm_vehicle_wrench_map(q,p)*wVehicle(:)+softarm_end_jacobian(q,p).'*wTip(:);
end
""",
        encoding="utf-8",
    )
    (target / "softarm_forward_dynamics.m").write_text(
        "function ddq=softarm_forward_dynamics(q,dq,tauArm,wVehicle,wTip,p)\n%#codegen\nM=softarm_mass(q,p); h=softarm_bias(q,dq,p); Q=softarm_applied_force(q,tauArm,wVehicle,wTip,p); ddq=M\\(Q-h);\nend\n",
        encoding="utf-8",
    )
    (target / "softarm_state_rhs.m").write_text(
        "function dx=softarm_state_rhs(x,tauArm,wVehicle,wTip,p)\n%#codegen\nn=numel(x)/2; q=x(1:n); dq=x(n+1:end); dx=[dq;softarm_forward_dynamics(q,dq,tauArm,wVehicle,wTip,p)];\nend\n",
        encoding="utf-8",
    )
    for filename, content in _HELPERS.items():
        (target / filename).write_text(content, encoding="utf-8")

    from .actuator_matlab import clear_actuator_functions, generate_actuator_matlab

    clear_actuator_functions(target)
    if actuation is not None:
        generate_actuator_matlab(plant, actuation, target, function_optimizer)
    from .constraint_matlab import clear_constraint_functions

    clear_constraint_functions(target)
    manifest = {
        "model": {
            "rod": plant.config.rod,
            "parameterization": plant.config.parameterization,
            "segments": plant.config.segments,
            "base_mode": plant.config.base.mode,
            "mount_xyz": plant.config.base.mount_xyz,
            "mount_rpy": plant.config.base.mount_rpy,
            "centerline": True,
            "dynamics_formulation": "recursive",
        },
        "coordinates": {
            "base": plant.base_coordinate_names,
            "arm": plant.arm_coordinate_names,
        },
        "parameters": [
            {"name": item.name, "default": item.default}
            for item in plant.parameters
            + (() if actuation is None else actuation.parameters)
        ],
        "actuation": None if actuation is None else {
            "family": actuation.family,
            "acceleration": actuation.acceleration,
            "channels": [
                {"name": name, "kind": kind}
                for name, kind in zip(
                    actuation.channel_names, actuation.channel_kinds, strict=True
                )
            ],
        },
        "constraint": None,
    }
    (target / "manifest.json").write_text(
        json.dumps(manifest, indent=2), encoding="utf-8"
    )
    (target / "softarm_model.tex").write_text(
        "\\documentclass{article}\n\\usepackage{amsmath}\n\\begin{document}\n"
        "\\section*{Recursive Soft-Arm Dynamics}\n"
        "The generated inverse-dynamics routine propagates section boundary "
        "velocities and accelerations forward and pulls spatial wrenches backward. "
        "The mass matrix is obtained by unit-acceleration inverse-dynamics calls, "
        "and $h(q,\\dot q)$ by a zero-acceleration call.\n\\end{document}\n",
        encoding="utf-8",
    )
    return target
