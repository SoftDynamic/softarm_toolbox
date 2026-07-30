from __future__ import annotations

from dataclasses import dataclass, field
import math
from pathlib import Path
import re
import tomllib
from typing import Any


class ConfigError(ValueError):
    pass


@dataclass(frozen=True)
class IntegrationConfig:
    method: str = "analytic"
    order: int | None = None


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
class ModelConfig:
    family: str
    segments: int
    inertia: str = "distributed"
    integration: IntegrationConfig = field(default_factory=IntegrationConfig)
    parameters: dict[str, Any] = field(default_factory=dict)
    ritz_x: tuple[float, ...] | None = None
    ritz_y: tuple[float, ...] | None = None
    actuation: ActuationConfig | None = None
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


def load_config(path: str | Path) -> ModelConfig:
    source = Path(path).resolve()
    with source.open("rb") as stream:
        raw = tomllib.load(stream)

    model = raw.get("model", {})
    family = str(model.get("family", "")).lower()
    if not family:
        raise ConfigError("model.family is required")
    segments = int(model.get("segments", 0))
    if segments < 1:
        raise ConfigError("model.segments must be positive")
    inertia = str(model.get("inertia", "distributed")).lower()
    if family == "pcc" and inertia not in {"distributed", "lumped"}:
        raise ConfigError("PCC inertia must be 'distributed' or 'lumped'")
    if family == "euler":
        inertia = "distributed"

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

    ritz_x = ritz_y = None
    if family == "euler":
        ritz = raw.get("ritz", {})
        ritz_x = _coefficients(ritz.get("x"), "ritz.x")
        ritz_y = _coefficients(ritz.get("y"), "ritz.y")
        _validate_ritz(ritz_x, "ritz.x")
        _validate_ritz(ritz_y, "ritz.y")

    parameters = dict(raw.get("parameters", {}))
    actuation = _load_actuation(raw.get("actuation"), segments)
    return ModelConfig(
        family=family,
        segments=segments,
        inertia=inertia,
        integration=IntegrationConfig(method, order),
        parameters=parameters,
        ritz_x=ritz_x,
        ritz_y=ritz_y,
        actuation=actuation,
        source=source,
    )


def broadcast(value: Any, count: int, name: str) -> list[float]:
    if isinstance(value, (int, float)):
        return [float(value)] * count
    if isinstance(value, list) and len(value) == count:
        return [float(item) for item in value]
    raise ConfigError(f"parameter {name!r} must be a scalar or an array of length {count}")
