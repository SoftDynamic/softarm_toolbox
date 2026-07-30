from __future__ import annotations

import json
from pathlib import Path
from typing import Iterable

import sympy as sp
from sympy.printing.octave import OctaveCodePrinter

from ..backends import optimize
from ..actuation import ActuationModel
from ..derive import SymbolicPlant


class _MatlabPrinter(OctaveCodePrinter):
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
    backend: str,
    wolfram_kernel: str | None,
) -> None:
    expressions = optimize(_flatten(matrix), backend, wolfram_kernel)
    replacements, reduced = sp.cse(expressions, symbols=sp.numbered_symbols("t"), order="none")
    printer = _MatlabPrinter()
    lines = [f"function {output_name} = {name}({','.join(inputs)})", "% Generated from the SymPy model. Do not edit.", "%#codegen"]
    lines.extend(loads)
    lines.extend(f"{symbol} = {printer.doprint(expr)};" for symbol, expr in replacements)
    vector = ";".join(printer.doprint(expr) for expr in reduced)
    dimensions = ",".join(str(item) for item in shape)
    lines.append(f"{output_name} = reshape([{vector}],{dimensions});")
    lines.append("end")
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


_HELPERS = {
    "softarm_sinc_sqrt.m": """function y = softarm_sinc_sqrt(z)\n%#codegen\nif abs(z)<1e-8, y=1-z/6+z^2/120-z^3/5040+z^4/362880; else, s=sqrt(z); y=sin(s)/s; end\nend\n""",
    "softarm_sinc_sqrt_d.m": """function y = softarm_sinc_sqrt_d(z)\n%#codegen\nif abs(z)<1e-8, y=-1/6+z/60-z^2/1680+z^3/90720; else, s=sqrt(z); y=(s*cos(s)-sin(s))/(2*s^3); end\nend\n""",
    "softarm_sinc_sqrt_dd.m": """function y = softarm_sinc_sqrt_dd(z)\n%#codegen\nif abs(z)<1e-8, y=1/60-z/840+z^2/30240; else, s=sqrt(z); y=((3-z)*sin(s)-3*s*cos(s))/(4*s^5); end\nend\n""",
    "softarm_cosc_sqrt.m": """function y = softarm_cosc_sqrt(z)\n%#codegen\nif abs(z)<1e-8, y=1/2-z/24+z^2/720-z^3/40320+z^4/3628800; else, s=sqrt(z); y=(1-cos(s))/z; end\nend\n""",
    "softarm_cosc_sqrt_d.m": """function y = softarm_cosc_sqrt_d(z)\n%#codegen\nif abs(z)<1e-8, y=-1/24+z/360-z^2/13440+z^3/907200; else, s=sqrt(z); y=(s*sin(s)-2*(1-cos(s)))/(2*z^2); end\nend\n""",
    "softarm_cosc_sqrt_dd.m": """function y = softarm_cosc_sqrt_dd(z)\n%#codegen\nif abs(z)<1e-8, y=1/360-z/6720+z^2/302400; else, s=sqrt(z); y=(z*cos(s)-5*s*sin(s)+8-8*cos(s))/(4*z^3); end\nend\n""",
}


def generate_matlab_bundle(
    plant: SymbolicPlant,
    output: str | Path,
    backend: str = "sympy",
    wolfram_kernel: str | None = None,
    actuation: ActuationModel | None = None,
) -> Path:
    target = Path(output).resolve()
    target.mkdir(parents=True, exist_ok=True)
    q_loads = symbol_loads(plant.q, "q")
    dq_loads = symbol_loads(plant.dq, "dq")
    p_loads = symbol_loads([item.symbol for item in plant.parameters], "p")
    common = q_loads + p_loads
    render_function(target / "softarm_mass.m", "softarm_mass", "M", plant.mass, plant.mass.shape, ["q", "p"], common, backend, wolfram_kernel)
    render_function(target / "softarm_bias.m", "softarm_bias", "h", plant.bias, plant.bias.shape, ["q", "dq", "p"], q_loads + dq_loads + p_loads, backend, wolfram_kernel)
    render_function(
        target / "softarm_kinematics.m", "softarm_kinematics", "H", plant.kinematics,
        (4, 4, plant.config.segments), ["q", "p"], common, backend, wolfram_kernel,
    )
    render_function(
        target / "softarm_end_jacobian.m", "softarm_end_jacobian", "J", plant.end_jacobian,
        plant.end_jacobian.shape, ["q", "p"], common, backend, wolfram_kernel,
    )
    (target / "softarm_forward_dynamics.m").write_text(
        "function ddq = softarm_forward_dynamics(q,dq,tau,w,p)\n%#codegen\nM=softarm_mass(q,p); h=softarm_bias(q,dq,p); J=softarm_end_jacobian(q,p); ddq=M\\(tau+J.'*w-h);\nend\n",
        encoding="utf-8",
    )
    (target / "softarm_state_rhs.m").write_text(
        "function dx = softarm_state_rhs(x,tau,w,p)\n%#codegen\nn=numel(x)/2; q=x(1:n); dq=x(n+1:end); dx=[dq;softarm_forward_dynamics(q,dq,tau,w,p)];\nend\n",
        encoding="utf-8",
    )
    for filename, content in _HELPERS.items():
        (target / filename).write_text(content, encoding="utf-8")
    from .actuator_matlab import clear_actuator_functions, generate_actuator_matlab

    clear_actuator_functions(target)
    if actuation is not None:
        generate_actuator_matlab(plant, actuation, target, backend, wolfram_kernel)
    runtime_parameters = plant.parameters + (() if actuation is None else actuation.parameters)
    manifest = {
        "segments": plant.config.segments,
        "coordinates": plant.coordinate_names,
        "parameters": [{"name": item.name, "default": item.default} for item in runtime_parameters],
    }
    (target / "manifest.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    return target
