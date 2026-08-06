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


def _write_runtime(plant: RecursivePlant, target: Path, dof: int) -> None:
    segments = plant.config.segments
    nq = len(plant.q)
    nodes, weights = _quadrature(plant)
    node_text = ";".join(nodes)
    weight_text = ";".join(weights)
    mass_indices = [
        _parameter_index(plant, f"s{section}_mass")
        for section in range(1, segments + 1)
    ]
    ixx_indices = [
        _parameter_index(plant, f"s{section}_Ixx")
        for section in range(1, segments + 1)
    ]
    iyy_indices = [
        _parameter_index(plant, f"s{section}_Iyy")
        for section in range(1, segments + 1)
    ]
    izz_indices = [
        _parameter_index(plant, f"s{section}_Izz")
        for section in range(1, segments + 1)
    ]
    gravity_index = _parameter_index(plant, "gravity")
    tip_mass_index = _parameter_index(plant, "tip_mass")
    tip_ixx_index = _parameter_index(plant, "tip_Ixx")
    tip_iyy_index = _parameter_index(plant, "tip_Iyy")
    tip_izz_index = _parameter_index(plant, "tip_Izz")
    def indices(values: list[int]) -> str:
        return ";".join(str(value) for value in values)
    common = f"""n={nq}; T0=softarm_mount_transform(q,p); dT=zeros(4,4,n);
nodes=[{node_text}]; weights=[{weight_text}];
massIndex=[{indices(mass_indices)}]; ixxIndex=[{indices(ixx_indices)}]; iyyIndex=[{indices(iyy_indices)}]; izzIndex=[{indices(izz_indices)}];
"""
    (target / "softarm_mass.m").write_text(
        f"""function M=softarm_mass(q,p)
%#codegen
{common}M=zeros(n,n);
for section=1:{segments}
    first=(section-1)*{dof}+1; last=first+{dof}-1;
    for sample=1:numel(nodes)
        [F0,dF]=softarm_legacy_affine(section,p,nodes(sample));
        [material0,dMaterial]=softarm_legacy_step(T0,dT,F0,dF,first,last);
        [Jv,Jw]=softarm_legacy_jacobian(material0,dMaterial);
        inertia=diag([p(ixxIndex(section)),p(iyyIndex(section)),p(izzIndex(section))]);
        worldInertia=material0(1:3,1:3)*inertia*material0(1:3,1:3).';
        M=M+weights(sample)*(p(massIndex(section))*(Jv.'*Jv)+Jw.'*worldInertia*Jw);
    end
    [F0,dF]=softarm_legacy_affine(section,p,1);
    [T0,dT]=softarm_legacy_step(T0,dT,F0,dF,first,last);
end
[Jv,Jw]=softarm_legacy_jacobian(T0,dT);
tipInertia=diag([p({tip_ixx_index}),p({tip_iyy_index}),p({tip_izz_index})]);
worldTipInertia=T0(1:3,1:3)*tipInertia*T0(1:3,1:3).';
M=M+p({tip_mass_index})*(Jv.'*Jv)+Jw.'*worldTipInertia*Jw; M=(M+M.')/2;
end
""",
        encoding="utf-8",
    )
    (target / "softarm_bias.m").write_text(
        f"""function h=softarm_bias(q,dq,p)
%#codegen
{common}h=softarm_recursive_internal_force(q,dq,p);
for section=1:{segments}
    first=(section-1)*{dof}+1; last=first+{dof}-1;
    for sample=1:numel(nodes)
        [F0,dF]=softarm_legacy_affine(section,p,nodes(sample));
        [material0,dMaterial]=softarm_legacy_step(T0,dT,F0,dF,first,last);
        [Jv,~]=softarm_legacy_jacobian(material0,dMaterial);
        h=h-weights(sample)*p(massIndex(section))*p({gravity_index})*Jv(3,:).';
    end
    [F0,dF]=softarm_legacy_affine(section,p,1);
    [T0,dT]=softarm_legacy_step(T0,dT,F0,dF,first,last);
end
[Jv,~]=softarm_legacy_jacobian(T0,dT);
h=h-p({tip_mass_index})*p({gravity_index})*Jv(3,:).';
end
""",
        encoding="utf-8",
    )
    (target / "softarm_inverse_dynamics.m").write_text(
        "function tau=softarm_inverse_dynamics(q,dq,ddq,p)\n%#codegen\ntau=softarm_mass(q,p)*ddq(:)+softarm_bias(q,dq,p);\nend\n",
        encoding="utf-8",
    )
    (target / "softarm_legacy_step.m").write_text(
        """function [next0,dNext]=softarm_legacy_step(parent0,dParent,F0,dF,first,last)
%#codegen
n=size(dParent,3); dNext=zeros(4,4,n); next0=parent0*F0;
for coordinate=1:n, dNext(:,:,coordinate)=dParent(:,:,coordinate)*F0; end
for local=1:size(dF,3), dNext(:,:,first+local-1)=dNext(:,:,first+local-1)+parent0*dF(:,:,local); end
end
""",
        encoding="utf-8",
    )
    (target / "softarm_legacy_jacobian.m").write_text(
        """function [Jv,Jw]=softarm_legacy_jacobian(reference,dTransform)
%#codegen
n=size(dTransform,3); Jv=zeros(3,n); Jw=zeros(3,n); R0=reference(1:3,1:3);
for coordinate=1:n
    Jv(:,coordinate)=dTransform(1:3,4,coordinate);
    rate=dTransform(1:3,1:3,coordinate)*R0.'; skewRate=(rate-rate.')/2;
    Jw(:,coordinate)=[skewRate(3,2);skewRate(1,3);skewRate(2,1)];
end
end
""",
        encoding="utf-8",
    )
    (target / "softarm_kinematics.m").write_text(
        f"""function H=softarm_kinematics(q,p)
%#codegen
n={nq}; T0=softarm_mount_transform(q,p); dT=zeros(4,4,n); H=zeros(4,4,{segments});
for section=1:{segments}
    first=(section-1)*{dof}+1; last=first+{dof}-1; [F0,dF]=softarm_legacy_affine(section,p,1);
    [T0,dT]=softarm_legacy_step(T0,dT,F0,dF,first,last); current=T0;
    for coordinate=1:n, current=current+dT(:,:,coordinate)*q(coordinate); end
    H(:,:,section)=current;
end
end
""",
        encoding="utf-8",
    )
    (target / "softarm_end_jacobian.m").write_text(
        f"""function J=softarm_end_jacobian(q,p)
%#codegen
n={nq}; T0=softarm_mount_transform(q,p); dT=zeros(4,4,n);
for section=1:{segments}
    first=(section-1)*{dof}+1; last=first+{dof}-1; [F0,dF]=softarm_legacy_affine(section,p,1);
    [T0,dT]=softarm_legacy_step(T0,dT,F0,dF,first,last);
end
[Jv,Jw]=softarm_legacy_jacobian(T0,dT); J=[Jv;Jw];
end
""",
        encoding="utf-8",
    )
    (target / "softarm_centerline_at.m").write_text(
        f"""function P=softarm_centerline_at(q,p,xi)
%#codegen
n={nq}; T0=softarm_mount_transform(q,p); dT=zeros(4,4,n); P=zeros(3,{segments});
for section=1:{segments}
    first=(section-1)*{dof}+1; last=first+{dof}-1; [F0,dF]=softarm_legacy_affine(section,p,xi);
    [material0,dMaterial]=softarm_legacy_step(T0,dT,F0,dF,first,last); material=material0;
    for coordinate=1:n, material=material+dMaterial(:,:,coordinate)*q(coordinate); end
    P(:,section)=material(1:3,4); [F0,dF]=softarm_legacy_affine(section,p,1);
    [T0,dT]=softarm_legacy_step(T0,dT,F0,dF,first,last);
end
end
""",
        encoding="utf-8",
    )


def generate_legacy_ritz_core(
    plant: RecursivePlant,
    output: str | Path,
    optimizer: FunctionOptimizer,
) -> Path:
    if plant.config.base.mode != "fixed":
        raise ValueError("recursive legacy Ritz currently requires base.mode='fixed'")
    target = Path(output).resolve()
    target.mkdir(parents=True, exist_ok=True)
    dof = _generate_affine_template(plant, target, optimizer)
    mount = transform_rpy(
        tuple(sp.Float(str(value)) for value in plant.config.base.mount_xyz),
        tuple(sp.Float(str(value)) for value in plant.config.base.mount_rpy),
    )
    render_function(
        target / "softarm_mount_transform.m",
        "softarm_mount_transform",
        "H",
        mount,
        (4, 4),
        ["q", "p"],
        symbol_loads(plant.q, "q")
        + symbol_loads([item.symbol for item in plant.parameters], "p"),
        optimizer,
    )
    from . import recursive_matlab

    recursive_matlab._generate_internal_and_tip(plant, target, optimizer)
    recursive_matlab._write_runtime(plant, target, dof)
    _write_runtime(plant, target, dof)
    (target / "softarm_vehicle_wrench_map.m").write_text(
        f"function Bv=softarm_vehicle_wrench_map(q,p)\n%#codegen\nBv=zeros({len(plant.q)},6);\nend\n",
        encoding="utf-8",
    )
    return target
