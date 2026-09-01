from __future__ import annotations

from pathlib import Path

import sympy as sp

from ..config import IntegrationConfig
from ..geometry import transform_rpy
from ..integration import unit_gauss_rule
from ..models import RecursivePlant
from .matlab import _flatten, render_function, symbol_loads
from .optimization import FunctionOptimizer


def _pack(*matrices: sp.Matrix) -> sp.Matrix:
    return sp.Matrix([value for matrix in matrices for value in _flatten(matrix)])


def _parameter_index(plant: RecursivePlant, name: str) -> int:
    return next(
        index
        for index, item in enumerate(plant.parameters, start=1)
        if item.name == name
    )


def _generate_affine_template(
    plant: RecursivePlant,
    target: Path,
    optimizer: FunctionOptimizer,
) -> int:
    definition = plant.definition
    dof = len(plant.arm_q) // plant.config.segments
    local_q = plant.arm_q[:dof, 0]
    xi = sp.Symbol("xi", real=True, nonnegative=True)
    transform = definition.kinematics.transform(0, xi)
    zero = {coordinate: 0 for coordinate in local_q}
    reference = transform.subs(zero)
    derivatives = [transform.diff(coordinate).subs(zero) for coordinate in local_q]
    section_parameters = [
        item for item in plant.parameters if item.name.startswith("s1_")
    ]
    generic_p = sp.Matrix(
        sp.symbols(f"local_p1:{len(section_parameters) + 1}", real=True)
    )
    substitutions = dict(
        zip(
            (item.symbol for item in section_parameters),
            generic_p,
            strict=True,
        )
    )
    packed = _pack(reference, *derivatives).xreplace(substitutions)
    render_function(
        target / "softarm_legacy_affine_template.m",
        "softarm_legacy_affine_template",
        "y",
        packed,
        packed.shape,
        ["pLocal", "xi"],
        symbol_loads(generic_p, "pLocal"),
        optimizer,
    )
    parameter_indices = {
        item.name: index
        for index, item in enumerate(plant.parameters, start=1)
    }
    suffixes = [item.name.removeprefix("s1_") for item in section_parameters]
    cases = []
    for section in range(plant.config.segments):
        indices = ",".join(
            str(parameter_indices[f"s{section + 1}_{suffix}"])
            for suffix in suffixes
        )
        cases.append(
            f"    case {section + 1}, y=softarm_legacy_affine_template(p([{indices}]),xi);"
        )
    (target / "softarm_legacy_affine.m").write_text(
        "\n".join([
            "function [F0,dF]=softarm_legacy_affine(section,p,xi)",
            "%#codegen",
            "switch section",
            *cases,
            "    otherwise, error('softarm:InvalidSection','Invalid section index.');",
            "end",
            f"F0=reshape(y(1:16),4,4); dF=reshape(y(17:end),4,4,{dof});",
            "end",
            "",
        ]),
        encoding="utf-8",
    )
    return dof


def _generate_base_and_mount(
    plant: RecursivePlant,
    target: Path,
    optimizer: FunctionOptimizer,
) -> None:
    p_symbols = sp.Matrix([item.symbol for item in plant.parameters])
    mount = transform_rpy(
        tuple(sp.Float(str(value)) for value in plant.config.base.mount_xyz),
        tuple(sp.Float(str(value)) for value in plant.config.base.mount_rpy),
    )
    render_function(
        target / "softarm_affine_mount.m",
        "softarm_affine_mount",
        "H",
        mount,
        (4, 4),
        ["p"],
        symbol_loads(p_symbols, "p"),
        optimizer,
    )
    render_function(
        target / "softarm_mount_transform.m",
        "softarm_mount_transform",
        "H",
        plant.base_transform * mount,
        (4, 4),
        ["q", "p"],
        symbol_loads(plant.q, "q") + symbol_loads(p_symbols, "p"),
        optimizer,
    )
    nbase = len(plant.base_q)
    if not nbase:
        (target / "softarm_affine_base_jet.m").write_text(
            "function [B,dB,d2B]=softarm_affine_base_jet(q,p)\n"
            "%#codegen\nB=eye(4); dB=zeros(4,4,0); d2B=zeros(4,4,0,0);\n"
            "end\n",
            encoding="utf-8",
        )
        return
    derivatives = [
        plant.base_transform.diff(coordinate) for coordinate in plant.base_q
    ]
    second = [
        derivative.diff(coordinate)
        for derivative in derivatives
        for coordinate in plant.base_q
    ]
    packed = _pack(plant.base_transform, *derivatives, *second)
    render_function(
        target / "softarm_affine_base_jet_raw.m",
        "softarm_affine_base_jet_raw",
        "y",
        packed,
        packed.shape,
        ["q", "p"],
        symbol_loads(plant.q, "q") + symbol_loads(p_symbols, "p"),
        optimizer,
    )
    (target / "softarm_affine_base_jet.m").write_text(
        f"""function [B,dB,d2B]=softarm_affine_base_jet(q,p)
%#codegen
y=softarm_affine_base_jet_raw(q,p); cursor=1;
B=reshape(y(cursor:cursor+15),4,4); cursor=cursor+16;
dB=zeros(4,4,{nbase}); d2B=zeros(4,4,{nbase},{nbase});
for first=1:{nbase}
    dB(:,:,first)=reshape(y(cursor:cursor+15),4,4); cursor=cursor+16;
end
for first=1:{nbase}
    for second=1:{nbase}
        d2B(:,:,first,second)=reshape(y(cursor:cursor+15),4,4); cursor=cursor+16;
    end
end
end
""",
        encoding="utf-8",
    )


_AFFINE_HELPERS = {
    "softarm_affine_step.m": """function [next0,dNext]=softarm_affine_step(parent0,dParent,F0,dF,first,last)
%#codegen
n=size(dParent,3); dNext=zeros(4,4,n); next0=parent0*F0;
for coordinate=1:n, dNext(:,:,coordinate)=dParent(:,:,coordinate)*F0; end
for local=1:size(dF,3), dNext(:,:,first+local-1)=parent0*dF(:,:,local); end
end
""",
    "softarm_affine_kinematic_jet.m": """function [T,Jv,dJv,Jw,dJw,R0]=softarm_affine_kinematic_jet(B,dB,d2B,A0,dA,qa)
%#codegen
nb=size(dB,3); na=size(dA,3); n=nb+na; A=A0;
for arm=1:na, A=A+dA(:,:,arm)*qa(arm); end
T=B*A; R0=B(1:3,1:3)*A0(1:3,1:3); Jv=zeros(3,n); dJv=zeros(3,n,n);
D=zeros(3,3,n); dD=zeros(3,3,n,n); dR0=zeros(3,3,n);
for column=1:nb
    derivative=dB(:,:,column)*A; Jv(:,column)=derivative(1:3,4);
    D(:,:,column)=dB(1:3,1:3,column)*A0(1:3,1:3);
    dR0(:,:,column)=D(:,:,column);
end
for arm=1:na
    column=nb+arm; derivative=B*dA(:,:,arm); Jv(:,column)=derivative(1:3,4);
    D(:,:,column)=B(1:3,1:3)*dA(1:3,1:3,arm);
end
for direction=1:nb
    for column=1:nb
        derivative=d2B(:,:,column,direction)*A;
        dJv(:,column,direction)=derivative(1:3,4);
        dD(:,:,column,direction)=d2B(1:3,1:3,column,direction)*A0(1:3,1:3);
    end
    for arm=1:na
        column=nb+arm; derivative=dB(:,:,direction)*dA(:,:,arm);
        dJv(:,column,direction)=derivative(1:3,4);
        dD(:,:,column,direction)=dB(1:3,1:3,direction)*dA(1:3,1:3,arm);
    end
end
for direction=1:na
    for column=1:nb
        derivative=dB(:,:,column)*dA(:,:,direction);
        dJv(:,column,nb+direction)=derivative(1:3,4);
    end
end
Jw=zeros(3,n); dJw=zeros(3,n,n);
for column=1:n
    rate=D(:,:,column)*R0.'; skew=(rate-rate.')/2;
    Jw(:,column)=[skew(3,2);skew(1,3);skew(2,1)];
    for direction=1:nb
        rateDerivative=dD(:,:,column,direction)*R0.'+D(:,:,column)*dR0(:,:,direction).';
        skew=(rateDerivative-rateDerivative.')/2;
        dJw(:,column,direction)=[skew(3,2);skew(1,3);skew(2,1)];
    end
end
end
""",
    "softarm_affine_body_terms.m": """function [M,Mdot,gradK,gravityTau]=softarm_affine_body_terms(B,dB,d2B,A0,dA,qa,dq,mass,inertia,gravity)
%#codegen
[~,Jv,dJv,Jw,dJw,R0]=softarm_affine_kinematic_jet(B,dB,d2B,A0,dA,qa);
n=numel(dq); worldInertia=R0*inertia*R0.'; M=mass*(Jv.'*Jv)+Jw.'*worldInertia*Jw;
Mdot=zeros(n); gradK=zeros(n,1); nb=size(dB,3);
for direction=1:n
    dWorldInertia=zeros(3);
    if direction<=nb
        dR=dB(1:3,1:3,direction)*A0(1:3,1:3);
        dWorldInertia=dR*inertia*R0.'+R0*inertia*dR.';
    end
    dM=mass*(dJv(:,:,direction).'*Jv+Jv.'*dJv(:,:,direction))+dJw(:,:,direction).'*worldInertia*Jw+Jw.'*dWorldInertia*Jw+Jw.'*worldInertia*dJw(:,:,direction);
    Mdot=Mdot+dM*dq(direction); gradK(direction)=0.5*dq.'*dM*dq;
end
gravityTau=-mass*gravity*Jv(3,:).';
end
""",
}


def _quadrature(plant: RecursivePlant) -> tuple[list[str], list[str]]:
    config = plant.config
    degree = max(
        len(config.ritz_x or ()) - 1,
        len(config.ritz_y or ()) - 1,
        len(config.ritz_z or ()) - 1,
    )
    rule = unit_gauss_rule(IntegrationConfig("gauss", max(3, degree + 1)))
    nodes = [sp.octave_code(node) for node, _ in rule]
    weights = [sp.octave_code(weight) for _, weight in rule]
    return nodes, weights


def _write_affine_runtime(plant: RecursivePlant, target: Path, dof: int) -> None:
    segments = plant.config.segments
    nq = len(plant.q)
    nbase = len(plant.base_q)
    narm = len(plant.arm_q)
    nodes, weights = _quadrature(plant)
    node_text = ";".join(nodes)
    weight_text = ";".join(weights)
    parameter_indices = {
        item.name: index for index, item in enumerate(plant.parameters, start=1)
    }

    def section_indices(suffix: str) -> str:
        return ";".join(
            str(parameter_indices[f"s{section}_{suffix}"])
            for section in range(1, segments + 1)
        )

    mass_indices = section_indices("mass")
    ixx_indices = section_indices("Ixx")
    iyy_indices = section_indices("Iyy")
    izz_indices = section_indices("Izz")
    gravity_index = parameter_indices["gravity"]
    tip_mass_index = parameter_indices["tip_mass"]
    tip_ixx_index = parameter_indices["tip_Ixx"]
    tip_iyy_index = parameter_indices["tip_Iyy"]
    tip_izz_index = parameter_indices["tip_Izz"]
    vehicle = ""
    if nbase:
        vehicle = f"""
vehicleInertia=diag([p({parameter_indices['vehicle_Ixx']}),p({parameter_indices['vehicle_Iyy']}),p({parameter_indices['vehicle_Izz']})]);
[bodyM,bodyMdot,bodyGrad,bodyGravity]=softarm_affine_body_terms(B,dB,d2B,eye(4),zeros(4,4,{narm}),qa,dq,p({parameter_indices['vehicle_mass']}),vehicleInertia,p({gravity_index}));
M=M+bodyM; Mdot=Mdot+bodyMdot; gradK=gradK+bodyGrad; gravityTau=gravityTau+bodyGravity;
"""
    common = f"""q=q(:); dq=dq(:); qa=q({nbase + 1}:end);
[B,dB,d2B]=softarm_affine_base_jet(q,p); A0=softarm_affine_mount(p); dA=zeros(4,4,{narm});
nodes=[{node_text}]; weights=[{weight_text}];
massIndex=[{mass_indices}]; ixxIndex=[{ixx_indices}]; iyyIndex=[{iyy_indices}]; izzIndex=[{izz_indices}];
M=zeros({nq}); Mdot=zeros({nq}); gradK=zeros({nq},1); gravityTau=zeros({nq},1);
{vehicle}
for section=1:{segments}
    first=(section-1)*{dof}+1; last=first+{dof}-1;
    inertia=diag([p(ixxIndex(section)),p(iyyIndex(section)),p(izzIndex(section))]);
    for sample=1:numel(nodes)
        [F0,dF]=softarm_legacy_affine(section,p,nodes(sample));
        [material0,dMaterial]=softarm_affine_step(A0,dA,F0,dF,first,last);
        [bodyM,bodyMdot,bodyGrad,bodyGravity]=softarm_affine_body_terms(B,dB,d2B,material0,dMaterial,qa,dq,p(massIndex(section)),inertia,p({gravity_index}));
        M=M+weights(sample)*bodyM; Mdot=Mdot+weights(sample)*bodyMdot;
        gradK=gradK+weights(sample)*bodyGrad; gravityTau=gravityTau+weights(sample)*bodyGravity;
    end
    [F0,dF]=softarm_legacy_affine(section,p,1);
    [A0,dA]=softarm_affine_step(A0,dA,F0,dF,first,last);
end
tipInertia=diag([p({tip_ixx_index}),p({tip_iyy_index}),p({tip_izz_index})]);
[bodyM,bodyMdot,bodyGrad,bodyGravity]=softarm_affine_body_terms(B,dB,d2B,A0,dA,qa,dq,p({tip_mass_index}),tipInertia,p({gravity_index}));
M=M+bodyM; Mdot=Mdot+bodyMdot; gradK=gradK+bodyGrad; gravityTau=gravityTau+bodyGravity;
M=(M+M.')/2; h=Mdot*dq-gradK+gravityTau+softarm_recursive_internal_force(q,dq,p);
"""
    (target / "softarm_affine_components.m").write_text(
        "function [M,h]=softarm_affine_components(q,dq,p)\n%#codegen\n"
        + common
        + "end\n",
        encoding="utf-8",
    )
    (target / "softarm_inverse_dynamics.m").write_text(
        "function tau=softarm_inverse_dynamics(q,dq,ddq,p)\n%#codegen\n"
        "[M,h]=softarm_affine_components(q,dq,p); tau=M*ddq(:)+h;\nend\n",
        encoding="utf-8",
    )
    (target / "softarm_mass.m").write_text(
        f"function M=softarm_mass(q,p)\n%#codegen\n"
        f"n={nq}; z=zeros(n,1); reference=softarm_inverse_dynamics(q,z,z,p); M=zeros(n);\n"
        "for column=1:n, acceleration=zeros(n,1); acceleration(column)=1; "
        "M(:,column)=softarm_inverse_dynamics(q,z,acceleration,p)-reference; end\n"
        "M=(M+M.')/2;\nend\n",
        encoding="utf-8",
    )
    (target / "softarm_bias.m").write_text(
        f"function h=softarm_bias(q,dq,p)\n%#codegen\nh=softarm_inverse_dynamics(q,dq,zeros({nq},1),p);\nend\n",
        encoding="utf-8",
    )
    propagation = f"""[B,dB,d2B]=softarm_affine_base_jet(q,p); qa=q({nbase + 1}:end);
A0=softarm_affine_mount(p); dA=zeros(4,4,{narm});
"""
    (target / "softarm_kinematics.m").write_text(
        f"""function H=softarm_kinematics(q,p)
%#codegen
{propagation}H=zeros(4,4,{segments});
for section=1:{segments}
    first=(section-1)*{dof}+1; last=first+{dof}-1; [F0,dF]=softarm_legacy_affine(section,p,1);
    [A0,dA]=softarm_affine_step(A0,dA,F0,dF,first,last); A=A0;
    for arm=1:{narm}, A=A+dA(:,:,arm)*qa(arm); end
    H(:,:,section)=B*A;
end
end
""",
        encoding="utf-8",
    )
    (target / "softarm_centerline_at.m").write_text(
        f"""function P=softarm_centerline_at(q,p,xi)
%#codegen
{propagation}P=zeros(3,{segments});
for section=1:{segments}
    first=(section-1)*{dof}+1; last=first+{dof}-1; [F0,dF]=softarm_legacy_affine(section,p,xi);
    [material0,dMaterial]=softarm_affine_step(A0,dA,F0,dF,first,last); A=material0;
    for arm=1:{narm}, A=A+dMaterial(:,:,arm)*qa(arm); end
    T=B*A; P(:,section)=T(1:3,4); [F0,dF]=softarm_legacy_affine(section,p,1);
    [A0,dA]=softarm_affine_step(A0,dA,F0,dF,first,last);
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
    endpoint = f"""{propagation}
for section=1:{segments}
    first=(section-1)*{dof}+1; last=first+{dof}-1; [F0,dF]=softarm_legacy_affine(section,p,1);
    [A0,dA]=softarm_affine_step(A0,dA,F0,dF,first,last);
end
"""
    (target / "softarm_end_jacobian.m").write_text(
        "function J=softarm_end_jacobian(q,p)\n%#codegen\n"
        + endpoint
        + "[~,Jv,~,Jw,~,~]=softarm_affine_kinematic_jet(B,dB,d2B,A0,dA,qa); J=[Jv;Jw];\nend\n",
        encoding="utf-8",
    )
    (target / "softarm_tool_point_jet.m").write_text(
        "function [position,J,gamma,velocity]=softarm_tool_point_jet(q,dq,toolOffset,p)\n%#codegen\n"
        + endpoint
        + "tool=[eye(3),toolOffset(:);0 0 0 1]; A0=A0*tool;\n"
        + f"for arm=1:{narm}, dA(:,:,arm)=dA(:,:,arm)*tool; end\n"
        + "[T,J,dJ,~,~,~]=softarm_affine_kinematic_jet(B,dB,d2B,A0,dA,qa); "
        + f"Jdot=zeros(3,{nq}); for direction=1:{nq}, Jdot=Jdot+dJ(:,:,direction)*dq(direction); end\n"
        + "position=T(1:3,4); velocity=J*dq; gamma=Jdot*dq;\nend\n",
        encoding="utf-8",
    )
    if nbase:
        vehicle_map = (
            "[B,dB,d2B]=softarm_affine_base_jet(q,p); "
            f"[~,Jv,~,Jw,~,~]=softarm_affine_kinematic_jet(B,dB,d2B,eye(4),zeros(4,4,{narm}),q({nbase + 1}:end)); "
            "J=[Jv;Jw]; R=B(1:3,1:3); Bv=J.'*blkdiag(R,R);"
        )
    else:
        vehicle_map = f"Bv=zeros({nq},6);"
    (target / "softarm_vehicle_wrench_map.m").write_text(
        f"function Bv=softarm_vehicle_wrench_map(q,p)\n%#codegen\n{vehicle_map}\nend\n",
        encoding="utf-8",
    )
    for filename, content in _AFFINE_HELPERS.items():
        (target / filename).write_text(content, encoding="utf-8")


def generate_legacy_ritz_core(
    plant: RecursivePlant,
    output: str | Path,
    optimizer: FunctionOptimizer,
) -> Path:
    target = Path(output).resolve()
    target.mkdir(parents=True, exist_ok=True)
    dof = _generate_affine_template(plant, target, optimizer)
    _generate_base_and_mount(plant, target, optimizer)
    from . import recursive_matlab

    recursive_matlab._generate_internal_and_tip(plant, target, optimizer)
    _write_affine_runtime(plant, target, dof)
    return target
