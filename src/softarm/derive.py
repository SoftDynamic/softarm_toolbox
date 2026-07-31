from __future__ import annotations

import json
from collections.abc import Callable
from dataclasses import dataclass, replace

import sympy as sp

from .config import ConfigError, ModelConfig, _validate_ritz, broadcast
from .dynamics import SympyBatchDifferentiator, assemble_bias
from .geometry import (
    angular_jacobian,
    pcs_transform,
    polynomial,
    ritz_transform,
    transform_rpy,
)
from .integration import integrate_unit, unit_gauss_rule


@dataclass(frozen=True)
class RuntimeParameter:
    name: str
    symbol: sp.Symbol
    default: float


ModelBuilder = Callable[[ModelConfig], "SymbolicPlant"]
ModelValidator = Callable[[ModelConfig], None]


@dataclass(frozen=True)
class RegisteredModel:
    builder: ModelBuilder
    validator: ModelValidator | None = None


@dataclass
class SymbolicPlant:
    config: ModelConfig
    q: sp.Matrix
    dq: sp.Matrix
    base_q: sp.Matrix
    base_dq: sp.Matrix
    arm_q: sp.Matrix
    arm_dq: sp.Matrix
    parameters: tuple[RuntimeParameter, ...]
    mass: sp.Matrix
    potential: sp.Expr
    damping: sp.Matrix
    kinematics: sp.Matrix
    end_jacobian: sp.Matrix
    base_transform: sp.Matrix
    end_transform: sp.Matrix
    base_jacobian: sp.Matrix
    vehicle_wrench_map: sp.Matrix
    arm_force_map: sp.Matrix
    _bias: sp.Matrix | None = None

    @property
    def p(self) -> sp.Matrix:
        return sp.Matrix([item.symbol for item in self.parameters])

    @property
    def coordinate_names(self) -> list[str]:
        return [str(item) for item in self.q]

    @property
    def base_coordinate_names(self) -> list[str]:
        return [str(item) for item in self.base_q]

    @property
    def arm_coordinate_names(self) -> list[str]:
        return [str(item) for item in self.arm_q]

    @property
    def bias(self) -> sp.Matrix:
        if self._bias is None:
            self._bias = assemble_bias(self, SympyBatchDifferentiator())
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


def _arm_coordinates(
    config: ModelConfig, labels: tuple[str, ...]
) -> tuple[sp.Matrix, sp.Matrix]:
    names: list[str] = []
    velocity_names: list[str] = []
    for section in range(1, config.segments + 1):
        for label in labels:
            names.append(f"{label}{section}")
            velocity_names.append(f"d{label}{section}")
    return sp.Matrix(sp.symbols(" ".join(names), real=True)), sp.Matrix(sp.symbols(" ".join(velocity_names), real=True))


def _base_coordinates(config: ModelConfig) -> tuple[sp.Matrix, sp.Matrix]:
    if config.base.mode == "fixed":
        return sp.zeros(0, 1), sp.zeros(0, 1)
    q = sp.Matrix(sp.symbols("base_x base_y base_z base_roll base_pitch base_yaw", real=True))
    dq = sp.Matrix(sp.symbols(
        "dbase_x dbase_y dbase_z dbase_roll dbase_pitch dbase_yaw", real=True
    ))
    return q, dq


def _linearize_matrix(matrix: sp.Matrix, q: sp.Matrix) -> sp.Matrix:
    zero = {coordinate: 0 for coordinate in q}
    return matrix.applyfunc(
        lambda value: value.subs(zero)
        + sum(
            sp.diff(value, coordinate).subs(zero) * coordinate
            for coordinate in q
        )
    )


def _mixed_linear_angular_jacobian(
    rotation: sp.Matrix,
    q: sp.Matrix,
    linear_coordinates: sp.Matrix,
) -> sp.Matrix:
    """Spatial angular Jacobian exact in the base attitude and first order in arm shape."""
    zero = {coordinate: 0 for coordinate in linear_coordinates}
    reference_rotation = rotation.subs(zero)
    columns = []
    for coordinate in q:
        rate = rotation.diff(coordinate).subs(zero) * reference_rotation.T
        skew = (rate - rate.T) / 2
        columns.append(sp.Matrix([skew[2, 1], skew[0, 2], skew[1, 0]]))
    return sp.Matrix.hstack(*columns)


def _point_terms(
    H: sp.Matrix,
    q: sp.Matrix,
    mass: sp.Expr,
    inertia: sp.Matrix,
    linear_rotation: bool = False,
    linear_coordinates: sp.Matrix | None = None,
) -> tuple[sp.Matrix, sp.Expr]:
    position = H[:3, 3]
    rotation = H[:3, :3]
    jv = position.jacobian(q)
    if linear_rotation:
        selected_coordinates = linear_coordinates if linear_coordinates is not None else q
        jw = _mixed_linear_angular_jacobian(
            rotation, q, selected_coordinates
        )
        zero = {coordinate: 0 for coordinate in selected_coordinates}
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

    base_q, base_dq = _base_coordinates(config)
    arm_q, arm_dq = q, dq
    q = base_q.col_join(arm_q)
    dq = base_dq.col_join(arm_dq)
    nq = len(q)
    mass_matrix = sp.zeros(nq)
    potential = elastic
    if config.base.mode == "floating_rpy":
        vehicle_mass = parameters.global_value("vehicle_mass", 1.5)
        vehicle_ixx = parameters.global_value("vehicle_Ixx", 0.03)
        vehicle_iyy = parameters.global_value("vehicle_Iyy", 0.03)
        vehicle_izz = parameters.global_value("vehicle_Izz", 0.05)
        base_transform = transform_rpy(base_q[:3, 0], tuple(base_q[3:6, 0]))
        vehicle_inertia = sp.diag(vehicle_ixx, vehicle_iyy, vehicle_izz)
        vehicle_term, vehicle_height = _point_terms(
            base_transform, q, vehicle_mass, vehicle_inertia,
        )
        mass_matrix += vehicle_term
        potential -= vehicle_mass * gravity * vehicle_height
    else:
        base_transform = sp.eye(4)

    mount = transform_rpy(
        tuple(sp.Float(str(value)) for value in config.base.mount_xyz),
        tuple(sp.Float(str(value)) for value in config.base.mount_rpy),
    )
    base = base_transform * mount
    transforms: list[sp.Matrix] = []
    xi = sp.Symbol("xi", real=True, nonnegative=True)

    for section in range(config.segments):
        local_end = local_transform(section, sp.S.One)
        end = base * local_end
        if linear_kinematics:
            end = _linearize_matrix(end, arm_q)
        transforms.append(end)

        if distributed:
            if config.integration.method == "gauss":
                for node, weight in unit_gauss_rule(config.integration):
                    material = base * local_transform(section, node)
                    if linear_kinematics:
                        material = _linearize_matrix(
                            material, arm_q
                        )
                    point_mass, height = _point_terms(
                        material, q, masses[section], inertias[section],
                        linear_kinematics,
                        arm_q if linear_kinematics else None,
                    )
                    mass_matrix += weight * point_mass
                    potential -= (
                        masses[section] * gravity * weight * height
                    )
            else:
                material = base * local_transform(section, xi)
                if linear_kinematics:
                    material = _linearize_matrix(
                        material, arm_q
                    )
                point_mass, height = _point_terms(
                    material, q, masses[section], inertias[section],
                    linear_kinematics,
                    arm_q if linear_kinematics else None,
                )
                mass_matrix += integrate_unit(
                    point_mass, xi, config.integration,
                    f"section {section + 1} inertia"
                )
                potential -= masses[section] * gravity * integrate_unit(
                    height, xi, config.integration,
                    f"section {section + 1} gravity"
                )
        else:
            midpoint = base * local_transform(section, sp.Rational(1, 2))
            if linear_kinematics:
                midpoint = _linearize_matrix(midpoint, arm_q)
            point_mass, height = _point_terms(
                midpoint, q, masses[section], inertias[section],
                linear_kinematics, arm_q if linear_kinematics else None,
            )
            mass_matrix += point_mass
            potential -= masses[section] * gravity * height
        base = end

    end_position = base[:3, 3]
    end_rotation = base[:3, :3]
    jv_end = end_position.jacobian(q)
    jw_end = (
        _mixed_linear_angular_jacobian(
            end_rotation, q, arm_q
        )
        if linear_kinematics else angular_jacobian(end_rotation, q)
    )
    tip_inertia = sp.diag(tip_ixx, tip_iyy, tip_izz)
    if linear_kinematics:
        zero = {coordinate: 0 for coordinate in arm_q}
        tip_rotation = end_rotation.subs(zero)
    else:
        tip_rotation = end_rotation
    mass_matrix += tip_mass * (jv_end.T * jv_end) + jw_end.T * (tip_rotation * tip_inertia * tip_rotation.T) * jw_end
    potential -= tip_mass * gravity * end_position[2]

    kinematics = sp.Matrix.hstack(*transforms)
    end_jacobian = jv_end.col_join(jw_end)
    base_position = base_transform[:3, 3]
    base_rotation = base_transform[:3, :3]
    if len(base_q):
        base_jv = base_position.jacobian(q)
        base_jw = angular_jacobian(base_rotation, q)
        base_jacobian = base_jv.col_join(base_jw)
        wrench_rotation = sp.diag(base_rotation, base_rotation)
        vehicle_wrench_map = base_jacobian.T * wrench_rotation
    else:
        base_jacobian = sp.zeros(6, nq)
        vehicle_wrench_map = sp.zeros(nq, 6)
    arm_force_map = sp.zeros(nq, len(arm_q))
    if len(arm_q):
        arm_force_map[len(base_q):, :] = sp.eye(len(arm_q))
    full_damping = [sp.S.Zero] * len(base_q) + damping
    return SymbolicPlant(
        config, q, dq, base_q, base_dq, arm_q, arm_dq,
        tuple(parameters.items), mass_matrix, potential, sp.diag(*full_damping),
        kinematics, end_jacobian, base_transform, base, base_jacobian,
        vehicle_wrench_map, arm_force_map,
    )


def _kirchhoff_pcs_strains(
    bx: sp.Expr,
    by: sp.Expr,
    current_length: sp.Expr,
    reference_length: sp.Expr,
) -> tuple[sp.Matrix, sp.Matrix]:
    return (
        sp.Matrix([-by / reference_length, bx / reference_length, 0]),
        sp.Matrix([0, 0, current_length / reference_length]),
    )


def _derive_extensible_kirchhoff_pcs(config: ModelConfig) -> SymbolicPlant:
    q, dq = _arm_coordinates(config, ("bx", "by", "l"))
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
        kappa, nu = _kirchhoff_pcs_strains(
            q[offset], q[offset + 1], q[offset + 2], lengths[section]
        )
        return pcs_transform(kappa, nu, lengths[section], xi)

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
        config, q, dq, pb, local, lengths, masses, inertias, damping,
        elastic, config.inertia == "distributed"
    )


def _derive_euler_bernoulli_ritz(config: ModelConfig) -> SymbolicPlant:
    q, dq = _arm_coordinates(config, ("ax", "ay"))
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
        return ritz_transform(
            q[offset], q[offset + 1], sp.S.Zero, lengths[section], local_xi,
            psi_x.subs(xi, local_xi), psi_y.subs(xi, local_xi),
            sp.S.Zero,
            dpsi_x.subs(xi, local_xi), dpsi_y.subs(xi, local_xi),
        )

    curvature_x = sp.diff(psi_x, xi, 2)
    curvature_y = sp.diff(psi_y, xi, 2)
    integral_x = integrate_unit(
        curvature_x**2, xi, config.integration,
        "Euler x Ritz stiffness"
    )
    integral_y = integrate_unit(
        curvature_y**2, xi, config.integration,
        "Euler y Ritz stiffness"
    )
    elastic = sp.S.Zero
    damping: list[sp.Expr] = []
    for section in range(config.segments):
        offset = 2 * section
        elastic += sp.Rational(1, 2) * (
            eiy[section] * q[offset] ** 2 * integral_x / lengths[section] ** 3
            + eix[section] * q[offset + 1] ** 2 * integral_y / lengths[section] ** 3
        )
        damping.extend([dax[section], day[section]])
    return _derive_common(
        config, q, dq, pb, local, lengths, masses, inertias, damping,
        elastic, True, True
    )


def _derive_euler_bernoulli_pcs(config: ModelConfig) -> SymbolicPlant:
    q, dq = _arm_coordinates(config, ("bx", "by"))
    pb = _ParameterBuilder(config)
    lengths = pb.sections("length", 0.5)
    masses = pb.sections("mass", 0.2)
    ixx = pb.sections("Ixx", 0.002)
    iyy = pb.sections("Iyy", 0.002)
    izz = pb.sections("Izz", 0.0001)
    eix = pb.sections("EI_x", 1.2)
    eiy = pb.sections("EI_y", 1.2)
    dbx = pb.sections("d_bx", 0.05)
    dby = pb.sections("d_by", 0.05)
    inertias = [sp.diag(ixx[i], iyy[i], izz[i]) for i in range(config.segments)]

    def local(section: int, xi: sp.Expr) -> sp.Matrix:
        offset = 2 * section
        kappa, nu = _kirchhoff_pcs_strains(
            q[offset], q[offset + 1], lengths[section], lengths[section]
        )
        return pcs_transform(kappa, nu, lengths[section], xi)

    elastic = sp.S.Zero
    damping: list[sp.Expr] = []
    for section in range(config.segments):
        offset = 2 * section
        elastic += sp.Rational(1, 2) * (
            eiy[section] * q[offset] ** 2 / lengths[section]
            + eix[section] * q[offset + 1] ** 2 / lengths[section]
        )
        damping.extend([dbx[section], dby[section]])
    return _derive_common(
        config, q, dq, pb, local, lengths, masses, inertias, damping,
        elastic, config.inertia == "distributed",
    )


def _derive_extensible_kirchhoff_ritz(config: ModelConfig) -> SymbolicPlant:
    q, dq = _arm_coordinates(config, ("ax", "ay", "az"))
    pb = _ParameterBuilder(config)
    lengths = pb.sections("rest_length", 0.5)
    masses = pb.sections("mass", 0.2)
    ixx = pb.sections("Ixx", 0.002)
    iyy = pb.sections("Iyy", 0.002)
    izz = pb.sections("Izz", 0.0001)
    eix = pb.sections("EI_x", 1.2)
    eiy = pb.sections("EI_y", 1.2)
    ea = pb.sections("EA", 50.0)
    dax = pb.sections("d_ax", 0.05)
    day = pb.sections("d_ay", 0.05)
    daz = pb.sections("d_az", 0.1)
    inertias = [sp.diag(ixx[i], iyy[i], izz[i]) for i in range(config.segments)]
    xi = sp.Symbol("xi", real=True, nonnegative=True)
    psi_x = polynomial(config.ritz_x or (), xi)
    psi_y = polynomial(config.ritz_y or (), xi)
    psi_z = polynomial(config.ritz_z or (), xi)
    dpsi_x = sp.diff(psi_x, xi)
    dpsi_y = sp.diff(psi_y, xi)

    def local(section: int, local_xi: sp.Expr) -> sp.Matrix:
        offset = 3 * section
        return ritz_transform(
            q[offset], q[offset + 1], q[offset + 2], lengths[section], local_xi,
            psi_x.subs(xi, local_xi), psi_y.subs(xi, local_xi),
            psi_z.subs(xi, local_xi),
            dpsi_x.subs(xi, local_xi), dpsi_y.subs(xi, local_xi),
        )

    curvature_x = sp.diff(psi_x, xi, 2)
    curvature_y = sp.diff(psi_y, xi, 2)
    axial_strain = sp.diff(psi_z, xi)
    integral_x = integrate_unit(
        curvature_x**2, xi, config.integration,
        "extensible Kirchhoff x Ritz stiffness",
    )
    integral_y = integrate_unit(
        curvature_y**2, xi, config.integration,
        "extensible Kirchhoff y Ritz stiffness",
    )
    integral_z = integrate_unit(
        axial_strain**2, xi, config.integration,
        "extensible Kirchhoff z Ritz stiffness",
    )
    elastic = sp.S.Zero
    damping: list[sp.Expr] = []
    for section in range(config.segments):
        offset = 3 * section
        elastic += sp.Rational(1, 2) * (
            eiy[section] * q[offset] ** 2 * integral_x / lengths[section] ** 3
            + eix[section] * q[offset + 1] ** 2 * integral_y / lengths[section] ** 3
            + ea[section] * q[offset + 2] ** 2 * integral_z / lengths[section]
        )
        damping.extend([dax[section], day[section], daz[section]])
    return _derive_common(
        config, q, dq, pb, local, lengths, masses, inertias, damping,
        elastic, True, True,
    )


def _derive_cosserat_pcs(config: ModelConfig) -> SymbolicPlant:
    q, dq = _arm_coordinates(config, ("kx", "ky", "kz", "vx", "vy", "vz"))
    pb = _ParameterBuilder(config)
    lengths = pb.sections("length", 0.5)
    masses = pb.sections("mass", 0.2)
    ixx = pb.sections("Ixx", 0.002)
    iyy = pb.sections("Iyy", 0.002)
    izz = pb.sections("Izz", 0.0001)
    kappa0_x = pb.sections("kappa0_x", 0.0)
    kappa0_y = pb.sections("kappa0_y", 0.0)
    kappa0_z = pb.sections("kappa0_z", 0.0)
    nu0_x = pb.sections("nu0_x", 0.0)
    nu0_y = pb.sections("nu0_y", 0.0)
    nu0_z = pb.sections("nu0_z", 1.0)
    eix = pb.sections("EI_x", 1.2)
    eiy = pb.sections("EI_y", 1.2)
    gj = pb.sections("GJ", 0.2)
    gax = pb.sections("GA_x", 20.0)
    gay = pb.sections("GA_y", 20.0)
    ea = pb.sections("EA", 50.0)
    dkx = pb.sections("d_kx", 0.05)
    dky = pb.sections("d_ky", 0.05)
    dkz = pb.sections("d_kz", 0.02)
    dvx = pb.sections("d_vx", 0.1)
    dvy = pb.sections("d_vy", 0.1)
    dvz = pb.sections("d_vz", 0.1)
    inertias = [sp.diag(ixx[i], iyy[i], izz[i]) for i in range(config.segments)]

    def local(section: int, xi: sp.Expr) -> sp.Matrix:
        offset = 6 * section
        kappa = sp.Matrix([
            kappa0_x[section] + q[offset],
            kappa0_y[section] + q[offset + 1],
            kappa0_z[section] + q[offset + 2],
        ])
        nu = sp.Matrix([
            nu0_x[section] + q[offset + 3],
            nu0_y[section] + q[offset + 4],
            nu0_z[section] + q[offset + 5],
        ])
        return pcs_transform(kappa, nu, lengths[section], xi)

    elastic = sp.S.Zero
    damping: list[sp.Expr] = []
    for section in range(config.segments):
        offset = 6 * section
        elastic += lengths[section] * sp.Rational(1, 2) * (
            eix[section] * q[offset] ** 2
            + eiy[section] * q[offset + 1] ** 2
            + gj[section] * q[offset + 2] ** 2
            + gax[section] * q[offset + 3] ** 2
            + gay[section] * q[offset + 4] ** 2
            + ea[section] * q[offset + 5] ** 2
        )
        damping.extend([
            lengths[section] * dkx[section],
            lengths[section] * dky[section],
            lengths[section] * dkz[section],
            lengths[section] * dvx[section],
            lengths[section] * dvy[section],
            lengths[section] * dvz[section],
        ])
    return _derive_common(
        config, q, dq, pb, local, lengths, masses, inertias, damping, elastic,
        config.inertia == "distributed",
    )


def _validate_no_ritz_options(config: ModelConfig) -> None:
    if any(item is not None for item in (config.ritz_x, config.ritz_y, config.ritz_z)):
        raise ConfigError("ritz coefficients are only applicable to Ritz parameterization")


def _validate_pcs_inertia(config: ModelConfig) -> None:
    _validate_no_ritz_options(config)
    if config.inertia not in {"distributed", "lumped"}:
        raise ConfigError("PCS inertia must be 'distributed' or 'lumped'")


def _validate_ritz_common(config: ModelConfig) -> None:
    if config.inertia != "distributed":
        raise ConfigError("Ritz parameterization requires distributed inertia")
    if config.ritz_x is None or config.ritz_y is None:
        raise ConfigError("Ritz parameterization requires ritz.x and ritz.y")
    _validate_ritz(config.ritz_x, "ritz.x")
    _validate_ritz(config.ritz_y, "ritz.y")


def _validate_euler_bernoulli_ritz(config: ModelConfig) -> None:
    _validate_ritz_common(config)
    if config.ritz_z is not None:
        raise ConfigError("euler_bernoulli + ritz does not accept ritz.z")


def _validate_extensible_kirchhoff_ritz(config: ModelConfig) -> None:
    _validate_ritz_common(config)
    if config.ritz_z is None:
        raise ConfigError("extensible_kirchhoff + ritz requires ritz.z")
    if abs(config.ritz_z[0]) > 1e-12:
        raise ConfigError("ritz.z must satisfy psi(0)=0")
    if abs(sum(config.ritz_z) - 1.0) > 1e-10:
        raise ConfigError("ritz.z must be normalized so psi(1)=1")


def _validate_cosserat_pcs(config: ModelConfig) -> None:
    _validate_pcs_inertia(config)
    if config.inertia == "distributed":
        if config.integration.method != "gauss":
            raise ConfigError(
                "distributed cosserat + pcs inertia requires integration.method='gauss'"
            )
        if config.integration.order is None or config.integration.order < 2:
            raise ConfigError(
                "distributed cosserat + pcs inertia requires Gauss order at least 2"
            )
    elif config.integration.method != "analytic":
        raise ConfigError("integration is not applicable to lumped cosserat + pcs inertia")


_CACHE: dict[str, SymbolicPlant] = {}


_MODELS: dict[tuple[str, str], RegisteredModel] = {
    ("euler_bernoulli", "ritz"): RegisteredModel(
        _derive_euler_bernoulli_ritz, _validate_euler_bernoulli_ritz
    ),
    ("euler_bernoulli", "pcs"): RegisteredModel(
        _derive_euler_bernoulli_pcs, _validate_pcs_inertia
    ),
    ("extensible_kirchhoff", "ritz"): RegisteredModel(
        _derive_extensible_kirchhoff_ritz, _validate_extensible_kirchhoff_ritz
    ),
    ("extensible_kirchhoff", "pcs"): RegisteredModel(
        _derive_extensible_kirchhoff_pcs, _validate_pcs_inertia
    ),
    ("cosserat", "pcs"): RegisteredModel(
        _derive_cosserat_pcs, _validate_cosserat_pcs
    ),
}


def register_model(
    rod: str,
    parameterization: str,
    builder: ModelBuilder,
    *,
    validator: ModelValidator | None = None,
) -> None:
    """Register one supported rod and spatial-parameterization combination."""
    key = (rod, parameterization)
    if not rod or not parameterization or key in _MODELS:
        raise ValueError(f"model combination {key!r} is already registered or invalid")
    _MODELS[key] = RegisteredModel(builder, validator)


def derive(config: ModelConfig) -> SymbolicPlant:
    combination = (config.rod, config.parameterization)
    try:
        registered = _MODELS[combination]
    except KeyError as error:
        raise ValueError(f"unregistered model combination: {combination!r}") from error
    if registered.validator is not None:
        registered.validator(config)
    key = json.dumps({
        "rod": config.rod,
        "parameterization": config.parameterization,
        "segments": config.segments,
        "inertia": config.inertia,
        "integration": [config.integration.method, config.integration.order],
        "parameters": config.parameters,
        "ritz_x": config.ritz_x,
        "ritz_y": config.ritz_y,
        "ritz_z": config.ritz_z,
        "base": [config.base.mode, config.base.mount_xyz, config.base.mount_rpy],
    }, sort_keys=True)
    if key not in _CACHE:
        _CACHE[key] = registered.builder(config)
    # Actuation is deliberately not part of the expensive physical-model cache,
    # but callers must retain the actuation attached to their own configuration.
    return replace(_CACHE[key], config=config, _bias=None)
