from __future__ import annotations

from dataclasses import dataclass, replace
import json
from typing import Callable

import sympy as sp

from .config import ModelConfig, broadcast
from .geometry import angular_jacobian, euler_ritz_transform, pcc_transform, polynomial
from .integration import integrate_unit


@dataclass(frozen=True)
class RuntimeParameter:
    name: str
    symbol: sp.Symbol
    default: float


@dataclass
class SymbolicPlant:
    config: ModelConfig
    q: sp.Matrix
    dq: sp.Matrix
    parameters: tuple[RuntimeParameter, ...]
    mass: sp.Matrix
    potential: sp.Expr
    damping: sp.Matrix
    kinematics: sp.Matrix
    end_jacobian: sp.Matrix
    _bias: sp.Matrix | None = None

    @property
    def p(self) -> sp.Matrix:
        return sp.Matrix([item.symbol for item in self.parameters])

    @property
    def coordinate_names(self) -> list[str]:
        return [str(item) for item in self.q]

    @property
    def bias(self) -> sp.Matrix:
        if self._bias is None:
            kinetic_gradient = (self.dq.T * self.mass * self.dq)[0]
            coriolis = (self.mass * self.dq).jacobian(self.q) * self.dq - sp.Rational(1, 2) * sp.Matrix([kinetic_gradient]).jacobian(self.q).T
            conservative = sp.Matrix([self.potential]).jacobian(self.q).T
            self._bias = coriolis + conservative + self.damping * self.dq
        return self._bias


class _ParameterBuilder:
    def __init__(self, config: ModelConfig):
        self.config = config
        self.items: list[RuntimeParameter] = []

    def global_value(self, name: str, default: float) -> sp.Symbol:
        value = float(self.config.parameters.get(name, default))
        symbol = sp.Symbol(name, real=True)
        self.items.append(RuntimeParameter(name, symbol, value))
        return symbol

    def sections(self, name: str, default: float) -> list[sp.Symbol]:
        values = broadcast(self.config.parameters.get(name, default), self.config.segments, name)
        result = []
        for index, value in enumerate(values, start=1):
            item_name = f"s{index}_{name}"
            symbol = sp.Symbol(item_name, real=True)
            self.items.append(RuntimeParameter(item_name, symbol, value))
            result.append(symbol)
        return result


def _coordinates(config: ModelConfig) -> tuple[sp.Matrix, sp.Matrix]:
    names: list[str] = []
    velocity_names: list[str] = []
    labels = ("bx", "by", "l") if config.family == "pcc" else ("ax", "ay")
    for section in range(1, config.segments + 1):
        for label in labels:
            names.append(f"{label}{section}")
            velocity_names.append(f"d{label}{section}")
    return sp.Matrix(sp.symbols(" ".join(names), real=True)), sp.Matrix(sp.symbols(" ".join(velocity_names), real=True))


def _linearize_matrix(matrix: sp.Matrix, q: sp.Matrix) -> sp.Matrix:
    zero = {coordinate: 0 for coordinate in q}
    return matrix.applyfunc(
        lambda value: value.subs(zero) + sum(sp.diff(value, coordinate).subs(zero) * coordinate for coordinate in q)
    )


def _linear_angular_jacobian(rotation: sp.Matrix, q: sp.Matrix) -> sp.Matrix:
    columns = []
    for coordinate in q:
        rate = sp.diff(rotation, coordinate)
        columns.append(sp.Matrix([rate[2, 1] - rate[1, 2], rate[0, 2] - rate[2, 0], rate[1, 0] - rate[0, 1]]) / 2)
    return sp.Matrix.hstack(*columns)


def _point_terms(
    H: sp.Matrix, q: sp.Matrix, mass: sp.Expr, inertia: sp.Matrix, linear_rotation: bool = False
) -> tuple[sp.Matrix, sp.Expr]:
    position = H[:3, 3]
    rotation = H[:3, :3]
    jv = position.jacobian(q)
    if linear_rotation:
        jw = _linear_angular_jacobian(rotation, q)
        zero = {coordinate: 0 for coordinate in q}
        rotation_for_inertia = rotation.subs(zero)
    else:
        jw = angular_jacobian(rotation, q)
        rotation_for_inertia = rotation
    mass_term = mass * (jv.T * jv) + jw.T * (rotation_for_inertia * inertia * rotation_for_inertia.T) * jw
    return mass_term, position[2]


def _derive_common(
    config: ModelConfig,
    q: sp.Matrix,
    dq: sp.Matrix,
    parameters: _ParameterBuilder,
    local_transform: Callable[[int, sp.Expr], sp.Matrix],
    lengths: list[sp.Symbol],
    masses: list[sp.Symbol],
    inertias: list[sp.Matrix],
    damping: list[sp.Expr],
    elastic: sp.Expr,
    distributed: bool,
    linear_kinematics: bool = False,
) -> SymbolicPlant:
    gravity = parameters.global_value("gravity", 9.81)
    tip_mass = parameters.global_value("tip_mass", 0.1)
    tip_ixx = parameters.global_value("tip_Ixx", 0.01)
    tip_iyy = parameters.global_value("tip_Iyy", 0.01)
    tip_izz = parameters.global_value("tip_Izz", 0.01)

    nq = len(q)
    mass_matrix = sp.zeros(nq)
    potential = elastic
    base = sp.eye(4)
    transforms: list[sp.Matrix] = []
    xi = sp.Symbol("xi", real=True, nonnegative=True)

    for section in range(config.segments):
        local_end = local_transform(section, sp.S.One)
        end = base * local_end
        if linear_kinematics:
            end = _linearize_matrix(end, q)
        transforms.append(end)

        if distributed:
            material = base * local_transform(section, xi)
            if linear_kinematics:
                material = _linearize_matrix(material, q)
            point_mass, height = _point_terms(material, q, masses[section], inertias[section], linear_kinematics)
            mass_matrix += integrate_unit(point_mass, xi, config.integration, f"section {section + 1} inertia")
            potential -= masses[section] * gravity * integrate_unit(height, xi, config.integration, f"section {section + 1} gravity")
        else:
            midpoint = base * local_transform(section, sp.Rational(1, 2))
            if linear_kinematics:
                midpoint = _linearize_matrix(midpoint, q)
            point_mass, height = _point_terms(midpoint, q, masses[section], inertias[section], linear_kinematics)
            mass_matrix += point_mass
            potential -= masses[section] * gravity * height
        base = end

    end_position = base[:3, 3]
    end_rotation = base[:3, :3]
    jv_end = end_position.jacobian(q)
    jw_end = _linear_angular_jacobian(end_rotation, q) if linear_kinematics else angular_jacobian(end_rotation, q)
    tip_inertia = sp.diag(tip_ixx, tip_iyy, tip_izz)
    if linear_kinematics:
        zero = {coordinate: 0 for coordinate in q}
        tip_rotation = end_rotation.subs(zero)
    else:
        tip_rotation = end_rotation
    mass_matrix += tip_mass * (jv_end.T * jv_end) + jw_end.T * (tip_rotation * tip_inertia * tip_rotation.T) * jw_end
    potential -= tip_mass * gravity * end_position[2]

    kinematics = sp.Matrix.hstack(*transforms)
    end_jacobian = jv_end.col_join(jw_end)
    return SymbolicPlant(config, q, dq, tuple(parameters.items), mass_matrix, potential, sp.diag(*damping), kinematics, end_jacobian)


def _derive_pcc(config: ModelConfig) -> SymbolicPlant:
    q, dq = _coordinates(config)
    pb = _ParameterBuilder(config)
    lengths = pb.sections("rest_length", 0.5)
    masses = pb.sections("mass", 0.2)
    ixx = pb.sections("Ixx", 0.002)
    iyy = pb.sections("Iyy", 0.002)
    izz = pb.sections("Izz", 0.0001)
    kbx = pb.sections("k_bx", 1.2)
    kby = pb.sections("k_by", 1.2)
    kl = pb.sections("k_l", 50.0)
    dbx = pb.sections("d_bx", 0.05)
    dby = pb.sections("d_by", 0.05)
    dl = pb.sections("d_l", 0.1)
    inertias = [sp.diag(ixx[i], iyy[i], izz[i]) for i in range(config.segments)]

    def local(section: int, xi: sp.Expr) -> sp.Matrix:
        offset = 3 * section
        return pcc_transform(q[offset], q[offset + 1], q[offset + 2], xi)

    elastic = sp.S.Zero
    damping: list[sp.Expr] = []
    for section in range(config.segments):
        offset = 3 * section
        elastic += sp.Rational(1, 2) * (
            kbx[section] * q[offset] ** 2
            + kby[section] * q[offset + 1] ** 2
            + kl[section] * (q[offset + 2] - lengths[section]) ** 2
        )
        damping.extend([dbx[section], dby[section], dl[section]])
    return _derive_common(
        config, q, dq, pb, local, lengths, masses, inertias, damping, elastic, config.inertia == "distributed"
    )


def _derive_euler(config: ModelConfig) -> SymbolicPlant:
    q, dq = _coordinates(config)
    pb = _ParameterBuilder(config)
    lengths = pb.sections("length", 0.5)
    masses = pb.sections("mass", 0.2)
    ixx = pb.sections("Ixx", 0.002)
    iyy = pb.sections("Iyy", 0.002)
    izz = pb.sections("Izz", 0.0001)
    eix = pb.sections("EI_x", 1.2)
    eiy = pb.sections("EI_y", 1.2)
    dax = pb.sections("d_ax", 0.05)
    day = pb.sections("d_ay", 0.05)
    inertias = [sp.diag(ixx[i], iyy[i], izz[i]) for i in range(config.segments)]
    xi = sp.Symbol("xi", real=True, nonnegative=True)
    psi_x = polynomial(config.ritz_x or (), xi)
    psi_y = polynomial(config.ritz_y or (), xi)
    dpsi_x = sp.diff(psi_x, xi)
    dpsi_y = sp.diff(psi_y, xi)

    def local(section: int, local_xi: sp.Expr) -> sp.Matrix:
        offset = 2 * section
        return euler_ritz_transform(
            q[offset], q[offset + 1], lengths[section], local_xi,
            psi_x.subs(xi, local_xi), psi_y.subs(xi, local_xi),
            dpsi_x.subs(xi, local_xi), dpsi_y.subs(xi, local_xi),
        )

    curvature_x = sp.diff(psi_x, xi, 2)
    curvature_y = sp.diff(psi_y, xi, 2)
    integral_x = integrate_unit(curvature_x**2, xi, config.integration, "Euler x Ritz stiffness")
    integral_y = integrate_unit(curvature_y**2, xi, config.integration, "Euler y Ritz stiffness")
    elastic = sp.S.Zero
    damping: list[sp.Expr] = []
    for section in range(config.segments):
        offset = 2 * section
        elastic += sp.Rational(1, 2) * (
            eiy[section] * q[offset] ** 2 * integral_x / lengths[section] ** 3
            + eix[section] * q[offset + 1] ** 2 * integral_y / lengths[section] ** 3
        )
        damping.extend([dax[section], day[section]])
    return _derive_common(config, q, dq, pb, local, lengths, masses, inertias, damping, elastic, True, True)


_CACHE: dict[str, SymbolicPlant] = {}


_MODEL_BUILDERS: dict[str, Callable[[ModelConfig], SymbolicPlant]] = {
    "pcc": _derive_pcc,
    "euler": _derive_euler,
}


def register_model(name: str, builder: Callable[[ModelConfig], SymbolicPlant]) -> None:
    """Register a model builder that returns a SymPy-backed SymbolicPlant."""
    if not name or name in _MODEL_BUILDERS:
        raise ValueError(f"model family {name!r} is already registered or invalid")
    _MODEL_BUILDERS[name] = builder


def derive(config: ModelConfig) -> SymbolicPlant:
    key = json.dumps({
        "family": config.family,
        "segments": config.segments,
        "inertia": config.inertia,
        "integration": [config.integration.method, config.integration.order],
        "parameters": config.parameters,
        "ritz_x": config.ritz_x,
        "ritz_y": config.ritz_y,
    }, sort_keys=True)
    if key not in _CACHE:
        try:
            builder = _MODEL_BUILDERS[config.family]
        except KeyError as error:
            raise ValueError(f"unregistered model family: {config.family}") from error
        _CACHE[key] = builder(config)
    # Actuation is deliberately not part of the expensive physical-model cache,
    # but callers must retain the actuation attached to their own configuration.
    return replace(_CACHE[key], config=config)
