from __future__ import annotations

import math
import re
import tomllib
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any


class ConfigError(ValueError):
    pass


@dataclass(frozen=True)
class IntegrationConfig:
    method: str = "analytic"
    order: int | None = None


@dataclass(frozen=True)
class BaseConfig:
    mode: str = "fixed"
    mount_xyz: tuple[float, float, float] = (0.0, 0.0, 0.0)
    mount_rpy: tuple[float, float, float] = (0.0, 0.0, 0.0)


@dataclass(frozen=True)
class TendonSpanConfig:
    section: int
    radius: float
    angle: float


@dataclass(frozen=True)
class TendonChannelConfig:
    name: str
    kind: str
    spans: tuple[TendonSpanConfig, ...]


@dataclass(frozen=True)
class ActuationConfig:
    family: str
    acceleration: str
    channels: tuple[TendonChannelConfig, ...] = ()
    data: dict[str, Any] = field(default_factory=dict)


@dataclass(frozen=True)
class ConstraintConfig:
    family: str
    data: dict[str, Any] = field(default_factory=dict)


@dataclass(frozen=True)
class ModelConfig:
    rod: str
    parameterization: str
    segments: int
    inertia: str = "distributed"
    integration: IntegrationConfig = field(default_factory=IntegrationConfig)
    parameters: dict[str, Any] = field(default_factory=dict)
    ritz_x: tuple[float, ...] | None = None
    ritz_y: tuple[float, ...] | None = None
    ritz_z: tuple[float, ...] | None = None
    base: BaseConfig = field(default_factory=BaseConfig)
    actuation: ActuationConfig | None = None
    constraint: ConstraintConfig | None = None
    source: Path | None = None


def _coefficients(value: Any, name: str) -> tuple[float, ...]:
    if not isinstance(value, list) or not value:
        raise ConfigError(f"{name} must be a non-empty coefficient array")
    return tuple(float(item) for item in value)


def _validate_ritz(coefficients: tuple[float, ...], name: str) -> None:
    # psi(xi) = sum c[k] xi^k. A clamped base requires c0=c1=0.
    if len(coefficients) < 3:
        raise ConfigError(f"{name} must contain at least quadratic coefficients")
    if abs(coefficients[0]) > 1e-12 or abs(coefficients[1]) > 1e-12:
        raise ConfigError(f"{name} must satisfy psi(0)=psi'(0)=0")
    if abs(sum(coefficients) - 1.0) > 1e-10:
        raise ConfigError(f"{name} must be normalized so psi(1)=1")


def _load_actuation(raw: Any, segments: int) -> ActuationConfig | None:
    if raw is None:
        return None
    if not isinstance(raw, dict):
        raise ConfigError("actuation must be a TOML table")
    family = str(raw.get("family", "")).lower()
    if not family:
        raise ConfigError("actuation.family is required")
    acceleration = str(raw.get("acceleration", "")).lower()
    if acceleration not in {"strict", "none"}:
        raise ConfigError("actuation.acceleration must be explicitly 'strict' or 'none'")
    if family != "tendon":
        return ActuationConfig(family, acceleration, data=dict(raw))

    channels_raw = raw.get("channels")
    if not isinstance(channels_raw, list) or not channels_raw:
        raise ConfigError("tendon actuation requires at least one channel")
    channels: list[TendonChannelConfig] = []
    names: set[str] = set()
    for channel_index, channel_raw in enumerate(channels_raw, start=1):
        if not isinstance(channel_raw, dict):
            raise ConfigError(f"actuation channel {channel_index} must be a table")
        name = str(channel_raw.get("name", ""))
        if not re.fullmatch(r"[A-Za-z][A-Za-z0-9_]*", name):
            raise ConfigError(f"invalid actuation channel name: {name!r}")
        if name in names:
            raise ConfigError(f"duplicate actuation channel name: {name}")
        names.add(name)
        kind = str(channel_raw.get("kind", "")).lower()
        if kind not in {"unilateral", "signed"}:
            raise ConfigError(f"actuation channel {name!r} kind must be 'unilateral' or 'signed'")
        spans_raw = channel_raw.get("spans")
        if not isinstance(spans_raw, list) or not spans_raw:
            raise ConfigError(f"actuation channel {name!r} requires at least one span")
        spans: list[TendonSpanConfig] = []
        seen_sections: set[int] = set()
        for span_raw in spans_raw:
            if not isinstance(span_raw, dict):
                raise ConfigError(f"actuation channel {name!r} span must be a table")
            section = int(span_raw.get("section", 0))
            if section < 1 or section > segments:
                raise ConfigError(f"actuation channel {name!r} has invalid section {section}")
            if section in seen_sections:
                raise ConfigError(f"actuation channel {name!r} repeats section {section}")
            seen_sections.add(section)
            radius = float(span_raw.get("radius", 0.0))
            angle = float(span_raw.get("angle", float("nan")))
            if not math.isfinite(radius) or radius <= 0:
                raise ConfigError(f"actuation channel {name!r} radius must be positive and finite")
            if not math.isfinite(angle):
                raise ConfigError(f"actuation channel {name!r} angle must be finite radians")
            spans.append(TendonSpanConfig(section, radius, angle))
        channels.append(TendonChannelConfig(name, kind, tuple(spans)))
    return ActuationConfig(family, acceleration, tuple(channels), dict(raw))


def _vector3(value: Any, name: str, default: tuple[float, float, float]) -> tuple[float, float, float]:
    selected = default if value is None else value
    if not isinstance(selected, (list, tuple)) or len(selected) != 3:
        raise ConfigError(f"{name} must be an array of three finite numbers")
    result = tuple(float(item) for item in selected)
    if not all(math.isfinite(item) for item in result):
        raise ConfigError(f"{name} must be an array of three finite numbers")
    return result  # type: ignore[return-value]


def _load_base(raw: Any) -> BaseConfig:
    if raw is None:
        return BaseConfig()
    if not isinstance(raw, dict):
        raise ConfigError("base must be a TOML table")
    mode = str(raw.get("mode", "fixed")).lower()
    if mode not in {"fixed", "floating_rpy"}:
        raise ConfigError("base.mode must be 'fixed' or 'floating_rpy'")
    return BaseConfig(
        mode,
        _vector3(raw.get("mount_xyz"), "base.mount_xyz", (0.0, 0.0, 0.0)),
        _vector3(raw.get("mount_rpy"), "base.mount_rpy", (0.0, 0.0, 0.0)),
    )


def _load_constraint(raw: Any) -> ConstraintConfig | None:
    if raw is None:
        return None
    if not isinstance(raw, dict):
        raise ConfigError("constraint must be a TOML table")
    family = str(raw.get("family", "")).lower()
    if not family:
        raise ConfigError("constraint.family is required")
    if family == "plane_point_contact":
        _vector3(raw.get("tool_offset"), "constraint.tool_offset", (0.0, 0.0, 0.0))
        _vector3(raw.get("plane_point"), "constraint.plane_point", (0.0, 0.0, 0.0))
        normal = _vector3(raw.get("plane_normal"), "constraint.plane_normal", (0.0, 0.0, -1.0))
        if math.sqrt(sum(item * item for item in normal)) <= 1e-12:
            raise ConfigError("constraint.plane_normal must be nonzero")
        for key, default, lower, strict in (
            ("friction", 0.0, 0.0, False),
            ("friction_velocity", 0.01, 0.0, True),
            ("stabilization_frequency", 20.0, 0.0, True),
            ("stabilization_ratio", 1.0, 0.0, False),
        ):
            value = float(raw.get(key, default))
            if not math.isfinite(value) or value < lower or (strict and value <= lower):
                qualifier = "positive" if strict else "nonnegative"
                raise ConfigError(f"constraint.{key} must be finite and {qualifier}")
    return ConstraintConfig(family, dict(raw))


def load_config(path: str | Path) -> ModelConfig:
    source = Path(path).resolve()
    with source.open("rb") as stream:
        raw = tomllib.load(stream)

    model = raw.get("model", {})
    rod = str(model.get("rod", "")).lower()
    if not rod:
        raise ConfigError("model.rod is required")
    parameterization = str(model.get("parameterization", "")).lower()
    if not parameterization:
        raise ConfigError("model.parameterization is required")
    segments = int(model.get("segments", 0))
    if segments < 1:
        raise ConfigError("model.segments must be positive")
    inertia = str(model.get("inertia", "distributed")).lower()

    integration_raw = raw.get("integration", {})
    method = str(integration_raw.get("method", "analytic")).lower()
    if method not in {"analytic", "gauss"}:
        raise ConfigError("integration.method must be 'analytic' or 'gauss'")
    order = integration_raw.get("order")
    if method == "gauss":
        if order is None or int(order) < 1:
            raise ConfigError("Gauss integration requires a positive order")
        order = int(order)
    else:
        order = None
    ritz_x = ritz_y = ritz_z = None
    if "ritz" in raw:
        ritz = raw["ritz"]
        if not isinstance(ritz, dict):
            raise ConfigError("ritz must be a TOML table")
        if "x" in ritz:
            ritz_x = _coefficients(ritz["x"], "ritz.x")
        if "y" in ritz:
            ritz_y = _coefficients(ritz["y"], "ritz.y")
        if "z" in ritz:
            ritz_z = _coefficients(ritz["z"], "ritz.z")

    parameters = dict(raw.get("parameters", {}))
    base = _load_base(raw.get("base"))
    actuation = _load_actuation(raw.get("actuation"), segments)
    constraint = _load_constraint(raw.get("constraint"))
    return ModelConfig(
        rod=rod,
        parameterization=parameterization,
        segments=segments,
        inertia=inertia,
        integration=IntegrationConfig(method, order),
        parameters=parameters,
        ritz_x=ritz_x,
        ritz_y=ritz_y,
        ritz_z=ritz_z,
        base=base,
        actuation=actuation,
        constraint=constraint,
        source=source,
    )


def broadcast(value: Any, count: int, name: str) -> list[float]:
    if isinstance(value, (int, float)):
        return [float(value)] * count
    if isinstance(value, list) and len(value) == count:
        return [float(item) for item in value]
    raise ConfigError(f"parameter {name!r} must be a scalar or an array of length {count}")
