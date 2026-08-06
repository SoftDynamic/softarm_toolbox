from __future__ import annotations

import json
from collections.abc import Callable
from dataclasses import dataclass, replace

import sympy as sp

from .config import ConfigError, ModelConfig, _validate_ritz, broadcast
from .dynamics import SymbolicLagrangeAssembler
from .geometry import (
    pac_transform,
    pcs_transform,
    polynomial,
    ritz_transform,
)
from .integration import integrate_unit
from .modeling import (
    FunctionalSectionKinematics,
    ModelCoordinates,
    ModelDefinition,
    ModelDefinitionBuilder,
    SectionProperties,
)
from .models import PlantModel, RecursivePlant, RuntimeParameter, SymbolicPlant
from .recursive import RecursiveInverseDynamicsAssembler

ModelBuilder = Callable[[ModelConfig], "SymbolicPlant"]
ModelValidator = Callable[[ModelConfig], None]


@dataclass(frozen=True)
class _SymbolicDefinitionBuilder(ModelDefinitionBuilder):
    definition_factory: Callable[[ModelConfig], ModelDefinition]
    validator: ModelValidator | None = None

    def validate(self, config: ModelConfig) -> None:
        if self.validator is not None:
            self.validator(config)

    def define(self, config: ModelConfig) -> ModelDefinition:
        return self.definition_factory(config)


@dataclass(frozen=True)
class _LegacyPlantBuilder(ModelDefinitionBuilder):
    """Compatibility adapter for the established register_model API."""

    plant_factory: ModelBuilder
    validator: ModelValidator | None = None

    def validate(self, config: ModelConfig) -> None:
        if self.validator is not None:
            self.validator(config)

    def define(self, config: ModelConfig) -> ModelDefinition:
        raise NotImplementedError("legacy plant builders do not expose a model definition")

    def build(
        self,
        config: ModelConfig,
        assembler: SymbolicLagrangeAssembler,
    ) -> SymbolicPlant:
        self.validate(config)
        return self.plant_factory(config)


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


def _model_definition(
    config: ModelConfig,
    q: sp.Matrix,
    dq: sp.Matrix,
    parameters: _ParameterBuilder,
    local_transform: Callable[[int, sp.Expr], sp.Matrix],
    lengths: list[sp.Symbol],
    masses: list[sp.Symbol],
    inertias: list[sp.Matrix],
    damping: list[sp.Expr] | sp.Matrix,
    elastic: sp.Expr,
    distributed: bool,
    linear_kinematics: bool = False,
) -> ModelDefinition:
    return ModelDefinition(
        config=config,
        coordinates=ModelCoordinates(q, dq),
        parameters=tuple(parameters.items),
        kinematics=FunctionalSectionKinematics(local_transform),
        sections=SectionProperties(tuple(lengths), tuple(masses), tuple(inertias)),
        damping=tuple(damping) if isinstance(damping, list) else damping,
        elastic=elastic,
        distributed=distributed,
        linear_kinematics=linear_kinematics,
    )


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


def _pac_arm_coordinates(
    config: ModelConfig, extensible: bool
) -> tuple[sp.Matrix, sp.Matrix]:
    names: list[str] = []
    velocity_names: list[str] = []
    for section in range(1, config.segments + 1):
        section_names = [f"c0_{section}", f"c1_{section}", f"phi{section}"]
        if extensible:
            section_names.append(f"l{section}")
        names.extend(section_names)
        velocity_names.extend(f"d{name}" for name in section_names)
    return (
        sp.Matrix(sp.symbols(" ".join(names), real=True)),
        sp.Matrix(sp.symbols(" ".join(velocity_names), real=True)),
    )


def _euler_bernoulli_pcs_strains(
    bx: sp.Expr,
    by: sp.Expr,
    current_length: sp.Expr,
    reference_length: sp.Expr,
) -> tuple[sp.Matrix, sp.Matrix]:
    return (
        sp.Matrix([-by / reference_length, bx / reference_length, 0]),
        sp.Matrix([0, 0, current_length / reference_length]),
    )


def _define_extensible_euler_bernoulli_pcs(config: ModelConfig) -> ModelDefinition:
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
        kappa, nu = _euler_bernoulli_pcs_strains(
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
    return _model_definition(
        config, q, dq, pb, local, lengths, masses, inertias, damping,
        elastic, config.inertia == "distributed"
    )


def _define_euler_bernoulli_ritz(config: ModelConfig) -> ModelDefinition:
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
    return _model_definition(
        config, q, dq, pb, local, lengths, masses, inertias, damping,
        elastic, True, True
    )


def _define_euler_bernoulli_pcs(config: ModelConfig) -> ModelDefinition:
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
        kappa, nu = _euler_bernoulli_pcs_strains(
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
    return _model_definition(
        config, q, dq, pb, local, lengths, masses, inertias, damping,
        elastic, config.inertia == "distributed",
    )


def _define_extensible_euler_bernoulli_ritz(config: ModelConfig) -> ModelDefinition:
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
        "extensible Euler-Bernoulli x Ritz stiffness",
    )
    integral_y = integrate_unit(
        curvature_y**2, xi, config.integration,
        "extensible Euler-Bernoulli y Ritz stiffness",
    )
    integral_z = integrate_unit(
        axial_strain**2, xi, config.integration,
        "extensible Euler-Bernoulli z Ritz stiffness",
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
    return _model_definition(
        config, q, dq, pb, local, lengths, masses, inertias, damping,
        elastic, True, True,
    )


def _pac_hankel() -> sp.Matrix:
    return sp.Matrix([
        [sp.S.One, sp.Rational(1, 2)],
        [sp.Rational(1, 2), sp.Rational(1, 3)],
    ])


def _define_euler_bernoulli_pac(config: ModelConfig) -> ModelDefinition:
    q, dq = _pac_arm_coordinates(config, extensible=False)
    pb = _ParameterBuilder(config)
    lengths = pb.sections("length", 0.5)
    masses = pb.sections("mass", 0.2)
    ixx = pb.sections("Ixx", 0.002)
    iyy = pb.sections("Iyy", 0.002)
    izz = pb.sections("Izz", 0.0001)
    eix = pb.sections("EI_x", 1.2)
    eiy = pb.sections("EI_y", 1.2)
    gj = pb.sections("GJ", 0.2)
    dbx = pb.sections("d_bx", 0.05)
    dby = pb.sections("d_by", 0.05)
    dphi = pb.sections("d_phi", 0.02)
    inertias = [sp.diag(ixx[i], iyy[i], izz[i]) for i in range(config.segments)]

    def local(section: int, xi: sp.Expr) -> sp.Matrix:
        offset = 3 * section
        return pac_transform(
            q[offset], q[offset + 1], q[offset + 2], lengths[section], xi
        )

    hankel = _pac_hankel()
    elastic = sp.S.Zero
    damping = sp.zeros(len(q))
    for section in range(config.segments):
        offset = 3 * section
        c = q[offset:offset + 2, 0]
        phi = q[offset + 2]
        directional_stiffness = (
            eiy[section] * sp.cos(phi) ** 2
            + eix[section] * sp.sin(phi) ** 2
        ) / lengths[section]
        elastic += (
            sp.Rational(1, 2) * directional_stiffness * (c.T * hankel * c)[0]
            + sp.Rational(1, 2) * gj[section] * phi**2 / lengths[section]
        )
        directional_damping = (
            dbx[section] * sp.cos(phi) ** 2
            + dby[section] * sp.sin(phi) ** 2
        )
        damping[offset:offset + 2, offset:offset + 2] = directional_damping * hankel
        damping[offset + 2, offset + 2] = dphi[section]
    return _model_definition(
        config, q, dq, pb, local, lengths, masses, inertias, damping,
        elastic, config.inertia == "distributed",
    )


def _define_extensible_euler_bernoulli_pac(config: ModelConfig) -> ModelDefinition:
    q, dq = _pac_arm_coordinates(config, extensible=True)
    pb = _ParameterBuilder(config)
    lengths = pb.sections("rest_length", 0.5)
    masses = pb.sections("mass", 0.2)
    ixx = pb.sections("Ixx", 0.002)
    iyy = pb.sections("Iyy", 0.002)
    izz = pb.sections("Izz", 0.0001)
    kbx = pb.sections("k_bx", 1.2)
    kby = pb.sections("k_by", 1.2)
    kphi = pb.sections("k_phi", 0.2)
    kl = pb.sections("k_l", 50.0)
    dbx = pb.sections("d_bx", 0.05)
    dby = pb.sections("d_by", 0.05)
    dphi = pb.sections("d_phi", 0.02)
    dl = pb.sections("d_l", 0.1)
    inertias = [sp.diag(ixx[i], iyy[i], izz[i]) for i in range(config.segments)]

    def local(section: int, xi: sp.Expr) -> sp.Matrix:
        offset = 4 * section
        return pac_transform(
            q[offset], q[offset + 1], q[offset + 2], q[offset + 3], xi
        )

    hankel = _pac_hankel()
    elastic = sp.S.Zero
    damping = sp.zeros(len(q))
    for section in range(config.segments):
        offset = 4 * section
        c = q[offset:offset + 2, 0]
        phi = q[offset + 2]
        current_length = q[offset + 3]
        directional_stiffness = (
            kbx[section] * sp.cos(phi) ** 2
            + kby[section] * sp.sin(phi) ** 2
        )
        elastic += (
            sp.Rational(1, 2) * directional_stiffness * (c.T * hankel * c)[0]
            + sp.Rational(1, 2) * kphi[section] * phi**2
            + sp.Rational(1, 2) * kl[section]
            * (current_length - lengths[section]) ** 2
        )
        directional_damping = (
            dbx[section] * sp.cos(phi) ** 2
            + dby[section] * sp.sin(phi) ** 2
        )
        damping[offset:offset + 2, offset:offset + 2] = directional_damping * hankel
        damping[offset + 2, offset + 2] = dphi[section]
        damping[offset + 3, offset + 3] = dl[section]
    return _model_definition(
        config, q, dq, pb, local, lengths, masses, inertias, damping,
        elastic, config.inertia == "distributed",
    )


def _define_cosserat_pcs(config: ModelConfig) -> ModelDefinition:
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
    return _model_definition(
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


def _validate_extensible_euler_bernoulli_ritz(config: ModelConfig) -> None:
    _validate_ritz_common(config)
    if config.ritz_z is None:
        raise ConfigError("extensible_euler_bernoulli + ritz requires ritz.z")
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


def _validate_pac(config: ModelConfig) -> None:
    _validate_no_ritz_options(config)
    if config.inertia == "distributed":
        if config.integration.method != "gauss":
            raise ConfigError("distributed PAC inertia requires integration.method='gauss'")
        if config.integration.order is None or config.integration.order < 4:
            raise ConfigError("distributed PAC inertia requires Gauss order at least 4")
    elif config.inertia == "lumped":
        if config.integration.method != "analytic":
            raise ConfigError("integration is not applicable to lumped PAC inertia")
    else:
        raise ConfigError("PAC inertia must be 'distributed' or 'lumped'")


_CACHE: dict[str, PlantModel] = {}


_MODELS: dict[tuple[str, str], ModelDefinitionBuilder] = {
    ("euler_bernoulli", "ritz"): _SymbolicDefinitionBuilder(
        _define_euler_bernoulli_ritz, _validate_euler_bernoulli_ritz
    ),
    ("euler_bernoulli", "pcs"): _SymbolicDefinitionBuilder(
        _define_euler_bernoulli_pcs, _validate_pcs_inertia
    ),
    ("extensible_euler_bernoulli", "ritz"): _SymbolicDefinitionBuilder(
        _define_extensible_euler_bernoulli_ritz,
        _validate_extensible_euler_bernoulli_ritz,
    ),
    ("extensible_euler_bernoulli", "pcs"): _SymbolicDefinitionBuilder(
        _define_extensible_euler_bernoulli_pcs, _validate_pcs_inertia
    ),
    ("euler_bernoulli", "pac"): _SymbolicDefinitionBuilder(
        _define_euler_bernoulli_pac, _validate_pac
    ),
    ("extensible_euler_bernoulli", "pac"): _SymbolicDefinitionBuilder(
        _define_extensible_euler_bernoulli_pac, _validate_pac
    ),
    ("cosserat", "pcs"): _SymbolicDefinitionBuilder(
        _define_cosserat_pcs, _validate_cosserat_pcs
    ),
}


def register_model(
    rod: str,
    parameterization: str,
    builder: ModelBuilder | ModelDefinitionBuilder,
    *,
    validator: ModelValidator | None = None,
) -> None:
    """Register one supported rod and spatial-parameterization combination."""
    key = (rod, parameterization)
    if not rod or not parameterization or key in _MODELS:
        raise ValueError(f"model combination {key!r} is already registered or invalid")
    if isinstance(builder, ModelDefinitionBuilder):
        if validator is not None:
            raise ValueError(
                "a ModelDefinitionBuilder owns its validation; do not pass validator"
            )
        _MODELS[key] = builder
    else:
        _MODELS[key] = _LegacyPlantBuilder(builder, validator)


def derive(config: ModelConfig) -> PlantModel:
    combination = (config.rod, config.parameterization)
    try:
        builder = _MODELS[combination]
    except KeyError as error:
        raise ValueError(f"unregistered model combination: {combination!r}") from error
    key = json.dumps({
        "rod": config.rod,
        "parameterization": config.parameterization,
        "segments": config.segments,
        "dynamics": config.dynamics.formulation,
        "inertia": config.inertia,
        "integration": [config.integration.method, config.integration.order],
        "parameters": config.parameters,
        "ritz_x": config.ritz_x,
        "ritz_y": config.ritz_y,
        "ritz_z": config.ritz_z,
        "base": [config.base.mode, config.base.mount_xyz, config.base.mount_rpy],
    }, sort_keys=True)
    if key not in _CACHE:
        if config.dynamics.formulation == "symbolic_lagrange":
            assembler = SymbolicLagrangeAssembler()
        else:
            if isinstance(builder, _LegacyPlantBuilder):
                raise ValueError(
                    "legacy function model builders support only symbolic_lagrange"
                )
            assembler = RecursiveInverseDynamicsAssembler()
        plant = builder.build(config, assembler)
        _CACHE[key] = plant
    # Actuation is deliberately not part of the expensive physical-model cache,
    # but callers must retain the actuation attached to their own configuration.
    cached = _CACHE[key]
    if isinstance(cached, SymbolicPlant):
        return replace(cached, config=config, _bias=None)
    if isinstance(cached, RecursivePlant):
        return replace(
            cached,
            config=config,
            definition=replace(cached.definition, config=config),
        )
    raise TypeError("registered model builder returned an unsupported plant type")
