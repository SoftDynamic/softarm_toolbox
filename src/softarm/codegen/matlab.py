from __future__ import annotations

import json
from collections.abc import Iterable
from pathlib import Path

import sympy as sp
from sympy.printing.octave import OctaveCodePrinter

from ..actuation import ActuationModel
from ..constraints import ConstraintModel
from ..derive import SymbolicPlant
from .optimization import FunctionOptimizer


class _MatlabPrinter(OctaveCodePrinter):
    def _print_AffineCosMoment(self, expr):
        return "softarm_affine_cos_moment(" + ",".join(
            self._print(arg) for arg in expr.args
        ) + ")"

    def _print_AffineSinMoment(self, expr):
        return "softarm_affine_sin_moment(" + ",".join(
            self._print(arg) for arg in expr.args
        ) + ")"

    def _print_SincSqrt(self, expr):
        return f"softarm_sinc_sqrt({self._print(expr.args[0])})"

    def _print_SincSqrtD(self, expr):
        return f"softarm_sinc_sqrt_d({self._print(expr.args[0])})"

    def _print_SincSqrtDD(self, expr):
        return f"softarm_sinc_sqrt_dd({self._print(expr.args[0])})"

    def _print_CoscSqrt(self, expr):
        return f"softarm_cosc_sqrt({self._print(expr.args[0])})"

    def _print_CoscSqrtD(self, expr):
        return f"softarm_cosc_sqrt_d({self._print(expr.args[0])})"

    def _print_CoscSqrtDD(self, expr):
        return f"softarm_cosc_sqrt_dd({self._print(expr.args[0])})"

    def _print_Sinc3Sqrt(self, expr):
        return f"softarm_sinc3_sqrt({self._print(expr.args[0])})"

    def _print_Sinc3SqrtD(self, expr):
        return f"softarm_sinc3_sqrt_d({self._print(expr.args[0])})"

    def _print_Sinc3SqrtDD(self, expr):
        return f"softarm_sinc3_sqrt_dd({self._print(expr.args[0])})"


def _flatten(matrix: sp.Matrix) -> list[sp.Expr]:
    return [matrix[row, column] for column in range(matrix.cols) for row in range(matrix.rows)]


def symbol_loads(symbols: Iterable[sp.Symbol], source: str) -> list[str]:
    return [f"{symbol} = {source}({index});" for index, symbol in enumerate(symbols, start=1)]


def render_function(
    path: Path,
    name: str,
    output_name: str,
    matrix: sp.Matrix,
    shape: tuple[int, ...],
    inputs: list[str],
    loads: list[str],
    optimizer: FunctionOptimizer,
) -> None:
    optimized = optimizer.optimize(_flatten(matrix), shape)
    printer = _MatlabPrinter()
    lines = [f"function {output_name} = {name}({','.join(inputs)})", "% Generated from the SymPy model. Do not edit.", "%#codegen"]
    lines.extend(loads)
    lines.extend(
        f"{symbol} = {printer.doprint(expr)};"
        for symbol, expr in optimized.replacements
    )
    vector = ";".join(printer.doprint(expr) for expr in optimized.expressions)
    dimensions = ",".join(str(item) for item in shape)
    lines.append(f"{output_name} = reshape([{vector}],{dimensions});")
    lines.append("end")
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def generate_centerline_matlab(
    plant: SymbolicPlant,
    target: Path,
    loads: list[str],
    optimizer: FunctionOptimizer,
) -> bool:
    scalar_path = target / "softarm_centerline_at.m"
    wrapper_path = target / "softarm_centerline.m"
    if plant._material_coordinate is None or plant._material_kinematics is None:
        scalar_path.unlink(missing_ok=True)
        wrapper_path.unlink(missing_ok=True)
        return False

    positions = sp.Matrix.hstack(*(
        plant._material_kinematics[:3, 4 * section + 3]
        for section in range(plant.config.segments)
    ))
    render_function(
        scalar_path,
        "softarm_centerline_at",
        "P",
        positions,
        positions.shape,
        ["q", "p", str(plant._material_coordinate)],
        loads,
        optimizer,
    )
    segment_count = plant.config.segments
    wrapper_path.write_text(
        (
            "function P = softarm_centerline(q,p,xi)\n"
            "%SOFTARM_CENTERLINE Sample model-derived section centerlines.\n"
            "xi=double(xi(:).');\n"
            "assert(~isempty(xi)&&all(isfinite(xi))&&all(xi>=0)&&all(xi<=1),"
            "'softarm:InvalidMaterialCoordinate','xi must be finite and lie in [0,1].');\n"
            f"P=zeros(3,{segment_count}*numel(xi));\n"
            "for sample=1:numel(xi)\n"
            "    sectionPoints=softarm_centerline_at(q,p,xi(sample));\n"
            f"    for section=1:{segment_count}\n"
            "        P(:,(section-1)*numel(xi)+sample)=sectionPoints(:,section);\n"
            "    end\n"
            "end\n"
            "end\n"
        ),
        encoding="utf-8",
    )
    return True


_HELPERS = {
    "softarm_affine_cos_moment.m": """function y = softarm_affine_cos_moment(n,c0,c1,xi)\n%#codegen\n[y,~]=softarm_affine_moment_parts(n,c0,c1,xi);\nend\n""",
    "softarm_affine_sin_moment.m": """function y = softarm_affine_sin_moment(n,c0,c1,xi)\n%#codegen\n[~,y]=softarm_affine_moment_parts(n,c0,c1,xi);\nend\n""",
    "softarm_affine_moment_parts.m": """function [c,s] = softarm_affine_moment_parts(n,c0,c1,xi)\n%#codegen\nassert(n>=0&&n==floor(n));\nif xi==0, c=0; s=0; return; end\npreviousPreviousReal=0; previousPreviousImag=0;\npreviousReal=1; previousImag=0;\nc=xi^(n+1)/(n+1); s=0; smallTerms=0;\nfor degree=1:256\n    realCoefficient=(-c0*previousImag-c1*previousPreviousImag)/degree;\n    imagCoefficient=(c0*previousReal+c1*previousPreviousReal)/degree;\n    scale=xi^(degree+n+1)/(degree+n+1);\n    realTerm=realCoefficient*scale; imagTerm=imagCoefficient*scale;\n    c=c+realTerm; s=s+imagTerm;\n    if hypot(realTerm,imagTerm)<=2e-16*max(1,hypot(c,s))\n        smallTerms=smallTerms+1;\n        if smallTerms>=4, break; end\n    else\n        smallTerms=0;\n    end\n    previousPreviousReal=previousReal; previousPreviousImag=previousImag;\n    previousReal=realCoefficient; previousImag=imagCoefficient;\nend\nend\n""",
    "softarm_sinc_sqrt.m": """function y = softarm_sinc_sqrt(z)\n%#codegen\nif abs(z)<1e-8, y=1-z/6+z^2/120-z^3/5040+z^4/362880; else, s=sqrt(z); y=sin(s)/s; end\nend\n""",
    "softarm_sinc_sqrt_d.m": """function y = softarm_sinc_sqrt_d(z)\n%#codegen\nif abs(z)<1e-8, y=-1/6+z/60-z^2/1680+z^3/90720; else, s=sqrt(z); y=(s*cos(s)-sin(s))/(2*s^3); end\nend\n""",
    "softarm_sinc_sqrt_dd.m": """function y = softarm_sinc_sqrt_dd(z)\n%#codegen\nif abs(z)<1e-8, y=1/60-z/840+z^2/30240; else, s=sqrt(z); y=((3-z)*sin(s)-3*s*cos(s))/(4*s^5); end\nend\n""",
    "softarm_cosc_sqrt.m": """function y = softarm_cosc_sqrt(z)\n%#codegen\nif abs(z)<1e-8, y=1/2-z/24+z^2/720-z^3/40320+z^4/3628800; else, s=sqrt(z); y=(1-cos(s))/z; end\nend\n""",
    "softarm_cosc_sqrt_d.m": """function y = softarm_cosc_sqrt_d(z)\n%#codegen\nif abs(z)<1e-8, y=-1/24+z/360-z^2/13440+z^3/907200; else, s=sqrt(z); y=(s*sin(s)-2*(1-cos(s)))/(2*z^2); end\nend\n""",
    "softarm_cosc_sqrt_dd.m": """function y = softarm_cosc_sqrt_dd(z)\n%#codegen\nif abs(z)<1e-8, y=1/360-z/6720+z^2/302400; else, s=sqrt(z); y=(z*cos(s)-5*s*sin(s)+8-8*cos(s))/(4*z^3); end\nend\n""",
    "softarm_sinc3_sqrt.m": """function y = softarm_sinc3_sqrt(z)\n%#codegen\nif abs(z)<1e-8, y=1/6-z/120+z^2/5040-z^3/362880+z^4/39916800; else, y=(1-softarm_sinc_sqrt(z))/z; end\nend\n""",
    "softarm_sinc3_sqrt_d.m": """function y = softarm_sinc3_sqrt_d(z)\n%#codegen\nif abs(z)<1e-8, y=-1/120+z/2520-z^2/120960+z^3/9979200; else, a=softarm_sinc_sqrt(z); y=(a-1-z*softarm_sinc_sqrt_d(z))/z^2; end\nend\n""",
    "softarm_sinc3_sqrt_dd.m": """function y = softarm_sinc3_sqrt_dd(z)\n%#codegen\nif abs(z)<1e-8, y=1/2520-z/60480+z^2/3326400; else, a=softarm_sinc_sqrt(z); y=(2-2*a+2*z*softarm_sinc_sqrt_d(z)-z^2*softarm_sinc_sqrt_dd(z))/z^3; end\nend\n""",
}


def generate_matlab_bundle(
    plant: SymbolicPlant,
    output: str | Path,
    actuation: ActuationModel | None = None,
    constraint: ConstraintModel | None = None,
    tex_appendix: bool = False,
    optimizer: FunctionOptimizer | None = None,
) -> Path:
    function_optimizer = optimizer or FunctionOptimizer()
    target = Path(output).resolve()
    target.mkdir(parents=True, exist_ok=True)
    q_loads = symbol_loads(plant.q, "q")
    dq_loads = symbol_loads(plant.dq, "dq")
    p_loads = symbol_loads([item.symbol for item in plant.parameters], "p")
    common = q_loads + p_loads
    render_function(target / "softarm_mass.m", "softarm_mass", "M", plant.mass, plant.mass.shape, ["q", "p"], common, function_optimizer)
    render_function(target / "softarm_bias.m", "softarm_bias", "h", plant.bias, plant.bias.shape, ["q", "dq", "p"], q_loads + dq_loads + p_loads, function_optimizer)
    render_function(
        target / "softarm_kinematics.m", "softarm_kinematics", "H", plant.kinematics,
        (4, 4, plant.config.segments), ["q", "p"], common, function_optimizer,
    )
    has_centerline = generate_centerline_matlab(
        plant, target, common, function_optimizer,
    )
    render_function(
        target / "softarm_end_jacobian.m", "softarm_end_jacobian", "J", plant.end_jacobian,
        plant.end_jacobian.shape, ["q", "p"], common, function_optimizer,
    )
    render_function(
        target / "softarm_vehicle_wrench_map.m", "softarm_vehicle_wrench_map", "Bv",
        plant.vehicle_wrench_map, plant.vehicle_wrench_map.shape, ["q", "p"], common,
        function_optimizer,
    )
    (target / "softarm_applied_force.m").write_text(
        (
            "function Q = softarm_applied_force(q,tauArm,wVehicle,wTip,p)\n"
            "%SOFTARM_APPLIED_FORCE Assemble arm, vehicle-body, and world-tip loads.\n"
            "%#codegen\n"
            f"assert(numel(tauArm)=={len(plant.arm_q)}); assert(numel(wVehicle)==6); assert(numel(wTip)==6);\n"
            f"Q=[zeros({len(plant.base_q)},1);tauArm(:)]+softarm_vehicle_wrench_map(q,p)*wVehicle(:)+softarm_end_jacobian(q,p).'*wTip(:);\n"
            "end\n"
        ),
        encoding="utf-8",
    )
    (target / "softarm_forward_dynamics.m").write_text(
        "function ddq = softarm_forward_dynamics(q,dq,tauArm,wVehicle,wTip,p)\n%#codegen\nM=softarm_mass(q,p); h=softarm_bias(q,dq,p); Q=softarm_applied_force(q,tauArm,wVehicle,wTip,p); ddq=M\\(Q-h);\nend\n",
        encoding="utf-8",
    )
    (target / "softarm_state_rhs.m").write_text(
        "function dx = softarm_state_rhs(x,tauArm,wVehicle,wTip,p)\n%#codegen\nn=numel(x)/2; q=x(1:n); dq=x(n+1:end); dx=[dq;softarm_forward_dynamics(q,dq,tauArm,wVehicle,wTip,p)];\nend\n",
        encoding="utf-8",
    )
    for filename, content in _HELPERS.items():
        (target / filename).write_text(content, encoding="utf-8")
    from .actuator_matlab import clear_actuator_functions, generate_actuator_matlab

    clear_actuator_functions(target)
    if actuation is not None:
        generate_actuator_matlab(plant, actuation, target, function_optimizer)
    from .constraint_matlab import clear_constraint_functions, generate_constraint_matlab

    clear_constraint_functions(target)
    if constraint is not None:
        generate_constraint_matlab(
            plant, constraint, target,
            () if actuation is None else actuation.parameters,
            function_optimizer,
        )
    runtime_parameters = (
        plant.parameters
        + (() if actuation is None else actuation.parameters)
        + (() if constraint is None else constraint.parameters)
    )
    actuation_manifest = None if actuation is None else {
        "family": actuation.family,
        "acceleration": actuation.acceleration,
        "channels": [
            {"name": name, "kind": kind}
            for name, kind in zip(actuation.channel_names, actuation.channel_kinds, strict=True)
        ],
    }
    constraint_manifest = None if constraint is None else {
        "family": constraint.family,
        "channels": [
            {"name": name, "kind": kind}
            for name, kind in zip(constraint.channel_names, constraint.channel_kinds, strict=True)
        ],
    }
    manifest = {
        "model": {
            "rod": plant.config.rod,
            "parameterization": plant.config.parameterization,
            "segments": plant.config.segments,
            "base_mode": plant.config.base.mode,
            "mount_xyz": plant.config.base.mount_xyz,
            "mount_rpy": plant.config.base.mount_rpy,
            "centerline": has_centerline,
        },
        "coordinates": {
            "base": plant.base_coordinate_names,
            "arm": plant.arm_coordinate_names,
        },
        "parameters": [{"name": item.name, "default": item.default} for item in runtime_parameters],
        "actuation": actuation_manifest,
        "constraint": constraint_manifest,
    }
    (target / "manifest.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    from .latex import generate_latex_document

    generate_latex_document(
        plant,
        target,
        actuation=actuation,
        constraint=constraint,
        include_appendix=tex_appendix,
        optimizer=function_optimizer,
    )
    return target
