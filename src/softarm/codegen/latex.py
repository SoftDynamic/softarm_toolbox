from __future__ import annotations

from pathlib import Path
import re
from typing import Iterable

import sympy as sp
from sympy.printing.latex import LatexPrinter

from ..actuation import ActuationModel
from ..backends import optimize
from ..constraints import ConstraintModel
from ..derive import RuntimeParameter, SymbolicPlant
from ..geometry import euler_ritz_transform, pcc_transform, polynomial, transform_rpy


def _escape_text(value: str) -> str:
    replacements = {
        "\\": r"\textbackslash{}",
        "&": r"\&",
        "%": r"\%",
        "$": r"\$",
        "#": r"\#",
        "_": r"\_",
        "{": r"\{",
        "}": r"\}",
        "~": r"\textasciitilde{}",
        "^": r"\textasciicircum{}",
    }
    return "".join(replacements.get(character, character) for character in value)


def _symbol_name(name: str) -> str:
    base = {
        "base_x": r"x_B", "base_y": r"y_B", "base_z": r"z_B",
        "base_roll": r"\phi_B", "base_pitch": r"\theta_B", "base_yaw": r"\psi_B",
        "dbase_x": r"\dot{x}_B", "dbase_y": r"\dot{y}_B", "dbase_z": r"\dot{z}_B",
        "dbase_roll": r"\dot{\phi}_B", "dbase_pitch": r"\dot{\theta}_B",
        "dbase_yaw": r"\dot{\psi}_B", "gravity": "g",
        "vehicle_mass": r"m_B", "vehicle_Ixx": r"I_{B,xx}",
        "vehicle_Iyy": r"I_{B,yy}", "vehicle_Izz": r"I_{B,zz}",
        "tip_mass": r"m_e", "tip_Ixx": r"I_{e,xx}",
        "tip_Iyy": r"I_{e,yy}", "tip_Izz": r"I_{e,zz}",
    }
    if name in base:
        return base[name]
    match = re.fullmatch(r"(d?)(bx|by|ax|ay|kx|ky|kz|vx|vy|vz|l)(\d+)", name)
    if match:
        derivative, token, section = match.groups()
        stem, axis = {
            "bx": ("b", "x"), "by": ("b", "y"),
            "ax": ("a", "x"), "ay": ("a", "y"), "l": ("L", None),
            "kx": (r"\delta\kappa", "x"), "ky": (r"\delta\kappa", "y"),
            "kz": (r"\delta\kappa", "z"), "vx": (r"\delta\nu", "x"),
            "vy": (r"\delta\nu", "y"), "vz": (r"\delta\nu", "z"),
        }[token]
        value = rf"{stem}_{{{section}}}" if axis is None else rf"{stem}_{{{axis},{section}}}"
        return rf"\dot{{{value}}}" if derivative else value
    match = re.fullmatch(r"s(\d+)_(.+)", name)
    if match:
        section, token = match.groups()
        labels = {
            "rest_length": ("L", "0"), "length": ("L", None), "mass": ("m", None),
            "Ixx": ("I", "xx"), "Iyy": ("I", "yy"), "Izz": ("I", "zz"),
            "k_bx": ("k", "b_x"), "k_by": ("k", "b_y"), "k_l": ("k", "L"),
            "d_bx": ("d", "b_x"), "d_by": ("d", "b_y"), "d_l": ("d", "L"),
            "EI_x": ("EI", "x"), "EI_y": ("EI", "y"), "d_ax": ("d", "a_x"),
            "d_ay": ("d", "a_y"),
            "kappa0_x": (r"\kappa_0", "x"), "kappa0_y": (r"\kappa_0", "y"),
            "kappa0_z": (r"\kappa_0", "z"), "nu0_x": (r"\nu_0", "x"),
            "nu0_y": (r"\nu_0", "y"), "nu0_z": (r"\nu_0", "z"),
            "GJ": ("GJ", None), "GA_x": ("GA", "x"), "GA_y": ("GA", "y"),
            "EA": ("EA", None), "d_kx": ("d", r"\kappa_x"),
            "d_ky": ("d", r"\kappa_y"), "d_kz": ("d", r"\kappa_z"),
            "d_vx": ("d", r"\nu_x"), "d_vy": ("d", r"\nu_y"),
            "d_vz": ("d", r"\nu_z"),
        }
        if token in labels:
            stem, detail = labels[token]
            index = section if detail is None else f"{detail},{section}"
            return rf"{stem}_{{{index}}}"
        return rf"\mathrm{{{_escape_text(token)}}}_{{{section}}}"
    match = re.fullmatch(r"act_(.+)_s(\d+)_radius", name)
    if match:
        channel, section = match.groups()
        return rf"r_{{\mathrm{{{_escape_text(channel)}}},{section}}}"
    match = re.fullmatch(r"constraint_(.+)_([xyz])", name)
    if match:
        token, axis = match.groups()
        return rf"{axis}_{{\mathrm{{{_escape_text(token)}}}}}"
    if name == "constraint_friction":
        return r"\mu"
    if name == "constraint_friction_velocity":
        return r"v_s"
    if name == "constraint_stabilization_frequency":
        return r"\omega_c"
    if name == "constraint_stabilization_ratio":
        return r"\zeta_c"
    match = re.fullmatch(r"cse_([a-z0-9]+)_(\d+)", name)
    if match:
        block, index = match.groups()
        return rf"c^{{\mathrm{{{block}}}}}_{{{index}}}"
    common = {
        "xi": r"\xi", "phi": r"\phi", "theta": r"\theta", "psi": r"\psi",
        "bx": r"b_x", "by": r"b_y", "ax": r"a_x", "ay": r"a_y", "L": "L",
    }
    if name in common:
        return common[name]
    return rf"\mathtt{{{_escape_text(name)}}}"


class SoftArmLatexPrinter(LatexPrinter):
    def _print_Symbol(self, expr: sp.Symbol) -> str:
        return _symbol_name(str(expr))

    def _special(self, name: str, argument: sp.Expr, derivative: str = "") -> str:
        return rf"\{name}{derivative}\!\left({self._print(argument)}\right)"

    def _print_SincSqrt(self, expr):
        return self._special("Sfun", expr.args[0])

    def _print_SincSqrtD(self, expr):
        return self._special("Sfun", expr.args[0], "'")

    def _print_SincSqrtDD(self, expr):
        return self._special("Sfun", expr.args[0], "''")

    def _print_CoscSqrt(self, expr):
        return self._special("Cfun", expr.args[0])

    def _print_CoscSqrtD(self, expr):
        return self._special("Cfun", expr.args[0], "'")

    def _print_CoscSqrtDD(self, expr):
        return self._special("Cfun", expr.args[0], "''")

    def _print_Sinc3Sqrt(self, expr):
        return self._special("Tfun", expr.args[0])

    def _print_Sinc3SqrtD(self, expr):
        return self._special("Tfun", expr.args[0], "'")

    def _print_Sinc3SqrtDD(self, expr):
        return self._special("Tfun", expr.args[0], "''")


def _math(expr: sp.Expr | sp.MatrixBase, printer: SoftArmLatexPrinter) -> str:
    return printer.doprint(expr)


def _equation(body: str, label: str | None = None) -> str:
    suffix = "" if label is None else rf"\label{{{label}}}"
    return "\\begin{equation}\n" + body + suffix + "\n\\end{equation}\n"


def _parameters(
    plant: SymbolicPlant,
    actuation: ActuationModel | None,
    constraint: ConstraintModel | None,
) -> tuple[RuntimeParameter, ...]:
    return (
        plant.parameters
        + (() if actuation is None else actuation.parameters)
        + (() if constraint is None else constraint.parameters)
    )


def _parameter_table(parameters: Iterable[RuntimeParameter], printer: SoftArmLatexPrinter) -> str:
    lines = [
        r"\begin{longtable}{lll}",
        r"\hline",
        r"Symbol & Code name & Default \\",
        r"\hline",
        r"\endfirsthead",
        r"\hline Symbol & Code name & Default \\",
        r"\hline",
        r"\endhead",
    ]
    for parameter in parameters:
        lines.append(
            rf"${_math(parameter.symbol, printer)}$ & \texttt{{{_escape_text(parameter.name)}}} & "
            rf"${parameter.default:.12g}$ \\"
        )
    lines.extend([r"\hline", r"\end{longtable}"])
    return "\n".join(lines)


def _channel_table(names: Iterable[str], kinds: Iterable[str]) -> str:
    lines = [
        r"\begin{center}",
        r"\begin{tabular}{ll}",
        r"\hline",
        r"Channel & Type \\",
        r"\hline",
    ]
    for name, kind in zip(names, kinds, strict=True):
        lines.append(
            rf"\texttt{{{_escape_text(name)}}} & \texttt{{{_escape_text(kind)}}} \\"
        )
    lines.extend([r"\hline", r"\end{tabular}", r"\end{center}"])
    return "\n".join(lines)


def _model_summary(plant: SymbolicPlant) -> str:
    source = "programmatic configuration"
    if plant.config.source is not None:
        source = plant.config.source.name
    integration = plant.config.integration.method
    if plant.config.integration.order is not None:
        integration += f" (order {plant.config.integration.order})"
    return "\n".join([
        r"\section{Model Summary}",
        r"This document is generated from the same SymPy model used for numerical code generation.",
        r"\begin{description}",
        rf"\item[Source configuration] \texttt{{{_escape_text(source)}}}",
        rf"\item[Arm family] \texttt{{{_escape_text(plant.config.family)}}}",
        rf"\item[Number of sections] {plant.config.segments}",
        rf"\item[Base mode] \texttt{{{_escape_text(plant.config.base.mode)}}}",
        rf"\item[Arm inertia model] \texttt{{{_escape_text(plant.config.inertia)}}}",
        rf"\item[Material integration] \texttt{{{_escape_text(integration)}}}",
        r"\end{description}",
    ])


def _notation(
    plant: SymbolicPlant,
    actuation: ActuationModel | None,
    constraint: ConstraintModel | None,
    printer: SoftArmLatexPrinter,
) -> str:
    lines = [r"\section{Notation and Parameters}"]
    if len(plant.base_q):
        lines.append(_equation(r"\boldsymbol q_B=" + _math(plant.base_q, printer)))
        lines.append(_equation(r"\dot{\boldsymbol q}_B=" + _math(plant.base_dq, printer)))
    else:
        lines.append(r"The arm base is fixed, hence $\boldsymbol q_B$ is empty.")
    lines.append(_equation(r"\boldsymbol q_a=" + _math(plant.arm_q, printer)))
    lines.append(_equation(r"\dot{\boldsymbol q}_a=" + _math(plant.arm_dq, printer)))
    lines.append(_equation(
        r"\boldsymbol q=\begin{bmatrix}\boldsymbol q_B\\\boldsymbol q_a\end{bmatrix},\qquad"
        r"\dot{\boldsymbol q}=\begin{bmatrix}\dot{\boldsymbol q}_B\\\dot{\boldsymbol q}_a\end{bmatrix}"
    ))
    lines.append(r"The NED world frame is used; gravity acts in the positive world $z$ direction.")
    lines.append(r"\subsection{Runtime Parameters}")
    lines.append(_parameter_table(_parameters(plant, actuation, constraint), printer))
    return "\n".join(lines)


def _kinematics(plant: SymbolicPlant, printer: SoftArmLatexPrinter) -> str:
    lines = [r"\section{Kinematics}"]
    if len(plant.base_q):
        lines.append(r"The floating-base attitude uses the ZYX convention $R_{WB}=R_z(\psi_B)R_y(\theta_B)R_x(\phi_B)$.")
        lines.append(_equation(r"H_{WB}=" + _math(plant.base_transform, printer), "eq:base-transform"))
    else:
        lines.append(_equation(r"H_{WB}=I_4"))
    mount = transform_rpy(
        tuple(sp.Float(str(value)) for value in plant.config.base.mount_xyz),
        tuple(sp.Float(str(value)) for value in plant.config.base.mount_rpy),
    )
    lines.append(_equation(r"H_{BA}=" + _math(mount, printer), "eq:mount-transform"))
    xi = sp.Symbol("xi", real=True, nonnegative=True)
    bx, by, length = sp.symbols("bx by L", real=True)
    if plant.config.family == "pcc":
        lines.extend([
            r"\subsection{PCC Section}",
            _equation(
                r"\Sfun(z)=\begin{cases}\dfrac{\sin\sqrt z}{\sqrt z},&z\ne0\\1,&z=0\end{cases},\qquad"
                r"\Cfun(z)=\begin{cases}\dfrac{1-\cos\sqrt z}{z},&z\ne0\\\dfrac12,&z=0\end{cases}"
            ),
            r"For one section, $b_x$ and $b_y$ are Cartesian bending-angle components and $L$ is its current length.",
            _equation(r"H_i(\xi)=" + _math(pcc_transform(bx, by, length, xi), printer), "eq:local-transform"),
        ])
    elif plant.config.family == "euler":
        ax, ay = sp.symbols("ax ay", real=True)
        psi_x = polynomial(plant.config.ritz_x or (), xi)
        psi_y = polynomial(plant.config.ritz_y or (), xi)
        lines.extend([
            r"\subsection{Euler--Bernoulli Ritz Section}",
            _equation(r"\psi_x(\xi)=" + _math(psi_x, printer) + r",\qquad\psi_y(\xi)=" + _math(psi_y, printer)),
            _equation(
                r"r_i(\xi)=\begin{bmatrix}a_x\psi_x(\xi)&a_y\psi_y(\xi)&L\xi\end{bmatrix}^T"
            ),
            _equation(
                r"\alpha_x(\xi)=\frac{a_x}{L}\psi_x'(\xi),\qquad"
                r"\alpha_y(\xi)=\frac{a_y}{L}\psi_y'(\xi)"
            ),
            _equation(
                r"H_i(\xi)=" + _math(euler_ritz_transform(
                    ax, ay, length, xi, psi_x, psi_y,
                    sp.diff(psi_x, xi), sp.diff(psi_y, xi),
                ), printer),
                "eq:local-transform",
            ),
            r"Products above first order in the arm Ritz coordinates are discarded; the floating-base attitude, when present, remains exact.",
        ])
    elif plant.config.family == "cosserat_pcs":
        lines.extend([
            r"\subsection{Cosserat Piecewise-Constant-Strain Section}",
            _equation(
                r"\Tfun(z)=\begin{cases}\dfrac{1-\Sfun(z)}{z},&z\ne0\\"
                r"\dfrac16,&z=0\end{cases}"
            ),
            r"The generalized coordinates are increments from the configured stress-free angular and linear strains.",
            _equation(
                r"\kappa_i=\kappa_{0,i}+\delta\kappa_i,\qquad"
                r"\nu_i=\nu_{0,i}+\delta\nu_i"
            ),
            _equation(
                r"\Omega_i(\xi)=L_i\xi\widehat{\kappa_i},\qquad"
                r"z_i(\xi)=(L_i\xi)^2\kappa_i^T\kappa_i"
            ),
            _equation(
                r"R_i=I_3+\Sfun(z_i)\Omega_i+\Cfun(z_i)\Omega_i^2,\qquad"
                r"r_i=\left[I_3+\Cfun(z_i)\Omega_i+\Tfun(z_i)\Omega_i^2\right]L_i\xi\nu_i"
            ),
            _equation(
                r"H_i(\xi)=\begin{bmatrix}R_i(\xi)&r_i(\xi)\\0&1\end{bmatrix}",
                "eq:local-transform",
            ),
        ])
    else:
        lines.append(r"This externally registered model family is documented from its public symbolic outputs.")
    lines.extend([
        _equation(r"H_{Wi}(\xi)=H_{WB}H_{BA}\left(\prod_{k=1}^{i-1}H_k(1)\right)H_i(\xi)"),
        _equation(
            r"J_{v,i}=\frac{\partial r_{Wi}}{\partial\boldsymbol q},\qquad"
            r"J_{\omega,i}^{(:,j)}=\operatorname{vex}\!\left(\operatorname{skew}\!\left("
            r"\frac{\partial R_{Wi}}{\partial q_j}R_{Wi}^{T}\right)\right)"
        ),
        _equation(r"J_e=\begin{bmatrix}J_{v,e}\\J_{\omega,e}\end{bmatrix}"),
    ])
    return "\n".join(lines)


def _energy_and_dynamics(plant: SymbolicPlant, printer: SoftArmLatexPrinter) -> str:
    if plant.config.inertia == "lumped":
        inertia_text = (
            r"For the lumped inertia option, each section contribution is evaluated at $\xi=1/2$."
        )
    else:
        inertia_text = (
            r"Each section uses distributed inertia over $\xi\in[0,1]$. "
            + (r"The configured integral is evaluated analytically."
               if plant.config.integration.method == "analytic"
               else rf"The configured {plant.config.integration.order}-point Gauss--Legendre rule is used.")
        )
        if plant.config.family == "euler":
            inertia_text += r" This is the inertia assembly for the configured Ritz coordinates."
    section_gravity = (
        r"\sum_{i=1}^{N}m_i z_i(1/2)"
        if plant.config.inertia == "lumped"
        else r"\sum_{i=1}^{N}m_i\int_0^1z_i(\xi)\,d\xi"
    )
    gravity_terms = section_gravity + r"+m_ez_e"
    if len(plant.base_q):
        gravity_terms = r"m_Bz_B+" + gravity_terms
    damping_diagonal = sp.Matrix([plant.damping[index, index] for index in range(len(plant.q))])
    lines = [
        r"\section{Energy and Dynamics}",
        inertia_text,
        _equation(
            r"M_i=\int_0^1\!\left(m_iJ_{v,i}^{T}J_{v,i}+"
            r"J_{\omega,i}^{T}R_{Wi}I_iR_{Wi}^{T}J_{\omega,i}\right)d\xi"
        ),
        _equation(
            r"M_e=m_eJ_{v,e}^{T}J_{v,e}+J_{\omega,e}^{T}R_eI_eR_e^TJ_{\omega,e}"
        ),
    ]
    if len(plant.base_q):
        lines.extend([
            _equation(
                r"M_B=J_B^T\operatorname{diag}\!\left(m_BI_3,"
                r"R_{WB}I_BR_{WB}^T\right)J_B"
            ),
            _equation(r"M=M_B+\sum_{i=1}^{N}M_i+M_e", "eq:mass-assembly"),
        ])
    else:
        lines.append(_equation(r"M=\sum_{i=1}^{N}M_i+M_e", "eq:mass-assembly"))
    lines.extend([
        _equation(r"T=\frac12\dot{\boldsymbol q}^{T}M(\boldsymbol q)\dot{\boldsymbol q}"),
        _equation(r"V=V_g+V_{\mathrm{elastic}},\qquad V_g=-g\!\left(" + gravity_terms + r"\right)"),
    ])
    if plant.config.family == "pcc":
        lines.append(_equation(
            r"V_{\mathrm{elastic}}=\frac12\sum_{i=1}^{N}\left("
            r"k_{b_x,i}b_{x,i}^{2}+k_{b_y,i}b_{y,i}^{2}+k_{L,i}(L_i-L_{0,i})^2\right)"
        ))
    elif plant.config.family == "euler":
        lines.append(_equation(
            r"V_{\mathrm{elastic}}=\frac12\sum_{i=1}^{N}\left["
            r"\frac{EI_{y,i}a_{x,i}^{2}}{L_i^3}\int_0^1(\psi_x'')^2d\xi+"
            r"\frac{EI_{x,i}a_{y,i}^{2}}{L_i^3}\int_0^1(\psi_y'')^2d\xi\right]"
        ))
    elif plant.config.family == "cosserat_pcs":
        lines.append(_equation(
            r"V_{\mathrm{elastic}}=\frac12\sum_{i=1}^{N}L_i\left("
            r"EI_{x,i}\delta\kappa_{x,i}^2+EI_{y,i}\delta\kappa_{y,i}^2+"
            r"GJ_i\delta\kappa_{z,i}^2+GA_{x,i}\delta\nu_{x,i}^2+"
            r"GA_{y,i}\delta\nu_{y,i}^2+EA_i\delta\nu_{z,i}^2\right)"
        ))
    lines.extend([
        _equation(r"D=\operatorname{diag}\!\left(" + _math(damping_diagonal, printer) + r"\right)"),
        _equation(
            r"c(\boldsymbol q,\dot{\boldsymbol q})="
            r"\frac{\partial(M\dot{\boldsymbol q})}{\partial\boldsymbol q}\dot{\boldsymbol q}-"
            r"\frac12\left(\frac{\partial(\dot{\boldsymbol q}^{T}M\dot{\boldsymbol q})}"
            r"{\partial\boldsymbol q}\right)^T"
        ),
        _equation(r"h=c+\nabla_{\boldsymbol q}V+D\dot{\boldsymbol q}"),
        r"The vehicle wrench $w_B$ is expressed in the vehicle body frame, whereas the end wrench $w_e$ is expressed in the NED world frame.",
        _equation(
            r"Q=S_a\tau_a+B_v(\boldsymbol q)w_B+J_e^Tw_e,\qquad"
            r"B_v=J_B^T\operatorname{diag}(R_{WB},R_{WB})"
        ),
        _equation(r"M\ddot{\boldsymbol q}+h=Q", "eq:dynamics"),
        _equation(
            r"\dot x=\begin{bmatrix}\dot{\boldsymbol q}\\M^{-1}(Q-h)\end{bmatrix},\qquad"
            r"x=\begin{bmatrix}\boldsymbol q\\\dot{\boldsymbol q}\end{bmatrix}"
        ),
    ])
    return "\n".join(lines)


def _actuation_section(actuation: ActuationModel, printer: SoftArmLatexPrinter) -> str:
    lines = [
        r"\section{Actuation}",
        rf"The configured actuator family is \texttt{{{_escape_text(actuation.family)}}} with {actuation.count} ordered channels.",
        _channel_table(actuation.channel_names, actuation.channel_kinds),
        _equation(r"y(\boldsymbol q_a)=" + _math(actuation.coordinates, printer)),
        _equation(
            r"J_a=\frac{\partial y}{\partial\boldsymbol q_a},\qquad"
            r"\gamma_a=\dot J_a\dot{\boldsymbol q}_a"
        ),
    ]
    if actuation.family == "tendon":
        lines.extend([
            r"The actuator coordinates are tendon lengths, and positive $T$ denotes tendon tension.",
            _equation(r"\tau_a=-J_a^TT"),
        ])
    if actuation.acceleration == "strict":
        lines.extend([
            r"Strict actuator-coordinate acceleration is enforced with the full-system Jacobian $J_{a,f}=[0\;J_a]$.",
            _equation(
                r"\begin{bmatrix}M&J_{a,f}^{T}\\J_{a,f}&0\end{bmatrix}"
                r"\begin{bmatrix}\ddot{\boldsymbol q}\\T\end{bmatrix}="
                r"\begin{bmatrix}Q-h\\\ddot y_{\mathrm{cmd}}-\gamma_a\end{bmatrix}"
            ),
        ])
    return "\n".join(lines)


def _constraint_section(constraint: ConstraintModel, printer: SoftArmLatexPrinter) -> str:
    lines = [
        r"\section{Active Acceleration Constraint}",
        rf"The configured constraint family is \texttt{{{_escape_text(constraint.family)}}} with {constraint.count} channels.",
        _channel_table(constraint.channel_names, constraint.channel_kinds),
        _equation(
            r"A=\frac{\partial\phi}{\partial\boldsymbol q},\qquad"
            r"\gamma=\dot A\dot{\boldsymbol q},\qquad Q_c=G(\boldsymbol q,\dot{\boldsymbol q})\lambda"
        ),
        _equation(
            r"\begin{bmatrix}M&-G\\A&0\end{bmatrix}"
            r"\begin{bmatrix}\ddot{\boldsymbol q}\\\lambda\end{bmatrix}="
            r"\begin{bmatrix}Q-h\\a_c-\gamma-2\zeta_c\omega_cA\dot{\boldsymbol q}-\omega_c^2\phi\end{bmatrix}"
        ),
    ]
    if constraint.family == "plane_point_contact":
        lines.extend([
            r"For the plane-point implementation, $n$ is the unit normal directed into the free half-space.",
            _equation(r"\phi=n^T(r_c-r_0)\ge0,\qquad A=n^TJ_c"),
            _equation(
                r"v_t=(I-nn^T)J_c\dot{\boldsymbol q},\qquad"
                r"f_t=-\mu\lambda\frac{v_t}{\sqrt{v_t^Tv_t+v_s^2}}"
            ),
            _equation(
                r"G=J_c^T\left(n-\mu\frac{v_t}{\sqrt{v_t^Tv_t+v_s^2}}\right),\qquad\lambda\ge0"
            ),
            r"The generated solver assumes that contact is active; a negative reaction marks the active mode as infeasible.",
        ])
    else:
        lines.append(_equation(r"\phi(\boldsymbol q)=" + _math(constraint.coordinates, printer)))
    return "\n".join(lines)


def _cse_data(
    expressions: list[sp.Expr],
    slug: str,
    backend: str,
    wolfram_kernel: str | None,
) -> tuple[list[tuple[sp.Symbol, sp.Expr]], list[sp.Expr]]:
    optimized = optimize(expressions, backend, wolfram_kernel)
    replacements, reduced = sp.cse(
        optimized, symbols=sp.numbered_symbols(f"cse_{slug}_"), order="canonical"
    )
    return list(replacements), list(reduced)


def _reconstruct_cse(
    replacements: list[tuple[sp.Symbol, sp.Expr]], reduced: list[sp.Expr]
) -> list[sp.Expr]:
    result = list(reduced)
    for symbol, expression in reversed(replacements):
        result = [item.xreplace({symbol: expression}) for item in result]
    return result


def _align(equations: list[str], chunk: int = 6) -> str:
    if not equations:
        return ""
    groups = []
    for start in range(0, len(equations), chunk):
        body = " \\\\\n".join(equations[start:start + chunk])
        groups.append("\\begin{align}\n" + body + "\n\\end{align}")
    return "\n".join(groups)


def _appendix_block(
    title: str,
    slug: str,
    matrix: sp.Matrix,
    symbol: str,
    printer: SoftArmLatexPrinter,
    backend: str,
    wolfram_kernel: str | None,
    upper: bool = False,
) -> str:
    indexed: list[tuple[int, int, sp.Expr]] = []
    for row in range(matrix.rows):
        for column in range(matrix.cols):
            if upper and column < row:
                continue
            indexed.append((row, column, matrix[row, column]))
    replacements, reduced = _cse_data(
        [item[2] for item in indexed], slug, backend, wolfram_kernel
    )
    definitions = [
        rf"{_math(temp, printer)}&={_math(expression, printer)}"
        for temp, expression in replacements
    ]
    values = []
    for (row, column, _), expression in zip(indexed, reduced, strict=True):
        subscript = str(row + 1) if matrix.cols == 1 else f"{row + 1},{column + 1}"
        values.append(rf"{symbol}_{{{subscript}}}&={_math(expression, printer)}")
    note = r"Only the upper triangle is listed; the lower triangle follows by symmetry." if upper else ""
    return "\n".join([
        rf"\subsection{{{_escape_text(title)}}}", note,
        r"\paragraph{Common subexpressions.}", _align(definitions) if definitions else r"No common subexpressions are required.",
        r"\paragraph{Elements.}", _align(values),
    ])


def _appendix(
    plant: SymbolicPlant,
    actuation: ActuationModel | None,
    constraint: ConstraintModel | None,
    printer: SoftArmLatexPrinter,
    backend: str,
    wolfram_kernel: str | None,
) -> str:
    blocks = [
        r"\appendix",
        r"\section{Exact Symbolic Appendix}",
        r"This appendix is generated from the exact symbolic outputs. Indices are one-based; the lower mass-matrix triangle follows by symmetry.",
        _equation(r"M_{ji}=M_{ij}\qquad(j>i)"),
        _appendix_block("Potential energy", "v", sp.Matrix([plant.potential]), "V", printer, backend, wolfram_kernel),
        _appendix_block("Damping matrix", "d", plant.damping, "D", printer, backend, wolfram_kernel),
        _appendix_block("Mass matrix", "m", plant.mass, "M", printer, backend, wolfram_kernel, upper=True),
        _appendix_block("Bias vector", "h", plant.bias, "h", printer, backend, wolfram_kernel),
        _appendix_block("End transform", "he", plant.end_transform, "H^e", printer, backend, wolfram_kernel),
        _appendix_block("End Jacobian", "je", plant.end_jacobian, "J^e", printer, backend, wolfram_kernel),
        _appendix_block("Vehicle wrench map", "bv", plant.vehicle_wrench_map, "B^v", printer, backend, wolfram_kernel),
    ]
    if actuation is not None:
        blocks.extend([
            _appendix_block("Actuator coordinates", "ay", actuation.coordinates, "y", printer, backend, wolfram_kernel),
            _appendix_block("Actuator Jacobian", "aja", actuation.jacobian, "J^a", printer, backend, wolfram_kernel),
            _appendix_block("Actuator velocity bias", "ag", actuation.velocity_bias, r"\gamma^a", printer, backend, wolfram_kernel),
        ])
    if constraint is not None:
        blocks.extend([
            _appendix_block("Constraint coordinates", "cp", constraint.coordinates, r"\phi", printer, backend, wolfram_kernel),
            _appendix_block("Constraint Jacobian", "ca", constraint.jacobian, "A", printer, backend, wolfram_kernel),
            _appendix_block("Constraint velocity bias", "cg", constraint.velocity_bias, r"\gamma", printer, backend, wolfram_kernel),
            _appendix_block("Constraint reaction map", "cr", constraint.reaction_map, "G", printer, backend, wolfram_kernel),
        ])
    return "\n".join(blocks)


def generate_latex_document(
    plant: SymbolicPlant,
    output: str | Path,
    actuation: ActuationModel | None = None,
    constraint: ConstraintModel | None = None,
    include_appendix: bool = False,
    backend: str = "sympy",
    wolfram_kernel: str | None = None,
) -> Path:
    """Generate a deterministic, standalone mathematical description of a model."""
    target = Path(output).resolve()
    target.mkdir(parents=True, exist_ok=True)
    printer = SoftArmLatexPrinter({"mat_delim": None, "mat_str": "bmatrix"})
    title = f"SoftArm {plant.config.family.upper()} Model ({plant.config.segments} Section"
    title += "s" if plant.config.segments != 1 else ""
    title += ")"
    sections = [
        r"\documentclass[11pt]{article}",
        r"\usepackage{amsmath,amssymb}",
        r"\usepackage[margin=24mm]{geometry}",
        r"\usepackage{longtable}",
        r"\allowdisplaybreaks",
        r"\newcommand{\Sfun}{\mathcal{S}}",
        r"\newcommand{\Cfun}{\mathcal{C}}",
        r"\newcommand{\Tfun}{\mathcal{T}}",
        rf"\title{{{_escape_text(title)}}}",
        r"\author{Generated by SoftArm Toolbox}",
        r"\date{}",
        r"\begin{document}",
        r"\maketitle",
        _model_summary(plant),
        _notation(plant, actuation, constraint, printer),
        _kinematics(plant, printer),
        _energy_and_dynamics(plant, printer),
    ]
    if actuation is not None:
        sections.append(_actuation_section(actuation, printer))
    if constraint is not None:
        sections.append(_constraint_section(constraint, printer))
    if include_appendix:
        sections.append(_appendix(
            plant, actuation, constraint, printer, backend, wolfram_kernel
        ))
    sections.append(r"\end{document}")
    path = target / "softarm_model.tex"
    path.write_text("\n\n".join(sections) + "\n", encoding="utf-8")
    return path
