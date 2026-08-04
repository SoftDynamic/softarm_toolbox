from __future__ import annotations

import re
import tomllib
from dataclasses import dataclass
from pathlib import Path
from typing import TYPE_CHECKING, Any

import sympy as sp

from .codegen.optimization import FunctionOptimizer, RenderedFunction
from .config import ModelConfig

if TYPE_CHECKING:
    from .derive import SymbolicPlant
    from .pipeline import BuildOptions, DerivedSystem


@dataclass(frozen=True)
class StageTiming:
    name: str
    seconds: float
    failed: bool


@dataclass(frozen=True)
class FunctionMetrics:
    name: str
    shape: tuple[int, ...]
    raw_operations: int
    temporary_count: int
    generated_operations: int
    source_bytes: int


@dataclass(frozen=True)
class ArtifactMetrics:
    files: tuple[Path, ...]
    counts: dict[str, int]
    sizes: dict[str, int]


def _format_value(value: Any) -> str:
    if isinstance(value, float):
        return f"{value:.12g}"
    if isinstance(value, (tuple, list)):
        return "[" + ", ".join(_format_value(item) for item in value) + "]"
    if value is None:
        return "none"
    return str(value)


def _format_shape(shape: tuple[int, ...]) -> str:
    return " x ".join(str(item) for item in shape)


def _format_source_size(size: int) -> str:
    if size < 1024:
        return f"{size} bytes"
    return f"{size / 1024:.1f} KiB"


def _setting(
    label: str,
    value: Any,
    *,
    indent: int = 0,
    width: int = 24,
) -> str:
    return f"{' ' * indent}{label:<{width}}: {_format_value(value)}"


def _with_default(value: Any, is_default: bool) -> str:
    suffix = " (default)" if is_default else ""
    return _format_value(value) + suffix


def _section(lines: list[str], title: str) -> None:
    if lines:
        lines.append("")
    lines.extend([title, "-" * len(title)])


def _operation_count(expressions: tuple[sp.Expr, ...]) -> int:
    return sum(int(sp.count_ops(expression, visual=False)) for expression in expressions)


class BuildDiagnostics:
    def __init__(
        self,
        config: ModelConfig,
        output: str | Path,
        options: BuildOptions,
    ):
        self.config = config
        self.output = Path(output).resolve()
        self.options = options
        self.timings: list[StageTiming] = []
        self.system: DerivedSystem | None = None
        self.plant: SymbolicPlant | None = None
        self.optimizer: FunctionOptimizer | None = None
        self.wolfram_executable: Path | None = None
        self.wolfram_version: str | None = None
        self.function_metrics: tuple[FunctionMetrics, ...] = ()
        self.artifact_metrics: ArtifactMetrics | None = None
        self._raw = self._load_raw_config()

    def _load_raw_config(self) -> dict[str, Any] | None:
        source = self.config.source
        if source is None or not source.is_file():
            return None
        try:
            with source.open("rb") as stream:
                value = tomllib.load(stream)
        except (OSError, tomllib.TOMLDecodeError):
            return None
        return value

    def record_timing(self, name: str, seconds: float, failed: bool) -> None:
        self.timings.append(StageTiming(name, seconds, failed))

    def capture_system(self, system: DerivedSystem) -> None:
        self.system = system
        self.plant = system.plant

    def capture_plant(self, plant: SymbolicPlant) -> None:
        self.plant = plant

    def capture_optimizer(self, optimizer: FunctionOptimizer) -> None:
        self.optimizer = optimizer

    def capture_wolfram(self, kernel: Any) -> None:
        executable = getattr(kernel, "executable", None)
        if executable is not None:
            self.wolfram_executable = Path(executable).resolve()
        version = getattr(kernel, "version", None)
        if version is not None:
            self.wolfram_version = str(version)

    def analyze(self) -> None:
        if self.optimizer is not None:
            self.function_metrics = tuple(
                self._function_metrics(item)
                for item in self.optimizer.rendered_functions
            )
        if self.output.is_dir():
            files = tuple(sorted(
                (path for path in self.output.iterdir() if path.is_file()),
                key=lambda path: path.name,
            ))
            counts: dict[str, int] = {}
            sizes: dict[str, int] = {}
            for path in files:
                suffix = path.suffix.lower()
                counts[suffix] = counts.get(suffix, 0) + 1
                sizes[suffix] = sizes.get(suffix, 0) + path.stat().st_size
            self.artifact_metrics = ArtifactMetrics(files, counts, sizes)

    def _function_metrics(self, item: RenderedFunction) -> FunctionMetrics:
        generated = tuple(
            expression for _, expression in item.optimized.replacements
        ) + item.optimized.expressions
        return FunctionMetrics(
            name=item.name,
            shape=item.optimized.shape,
            raw_operations=_operation_count(item.expressions),
            temporary_count=len(item.optimized.replacements),
            generated_operations=_operation_count(generated),
            source_bytes=item.path.stat().st_size,
        )

    def render(self, *, success: bool) -> str:
        lines = ["SoftArm Build Diagnostics", "=" * 25]
        self._configuration(lines)
        self._symbolic_strategy(lines)
        self._build_stages(lines)
        if self.system is not None:
            self._resolved_model(lines)
        if self.function_metrics:
            self._exported_functions(lines)
        if success and self.artifact_metrics is not None:
            self._generated_artifacts(lines)
        self._build_result(lines, success)
        return "\n".join(lines) + "\n"

    def _configuration(self, lines: list[str]) -> None:
        _section(lines, "Configuration")
        if self.config.source is not None:
            lines.append(_setting("Source file", self.config.source.resolve()))
        lines.extend([
            _setting("Output directory", self.output),
            _setting("Target", "MATLAB"),
            _setting("TeX appendix", "enabled" if self.options.tex_appendix else "disabled"),
            "",
            "Model",
            _setting("Rod", self.config.rod, indent=2),
            _setting("Parameterization", self.config.parameterization, indent=2),
            _setting("Sections", self.config.segments, indent=2),
        ])
        raw_model = None if self._raw is None else self._raw.get("model", {})
        inertia_default = raw_model is not None and "inertia" not in raw_model
        lines.append(_setting(
            "Inertia model",
            _with_default(self.config.inertia, inertia_default),
            indent=2,
        ))

        raw_base = None if self._raw is None else self._raw.get("base", {})
        orientation = _format_value(self.config.base.mount_rpy) + " rad"
        if raw_base is not None and "mount_rpy" not in raw_base:
            orientation += " (default)"
        lines.extend([
            "",
            "Base",
            _setting(
                "Mode",
                _with_default(
                    self.config.base.mode,
                    raw_base is not None and "mode" not in raw_base,
                ),
                indent=2,
            ),
            _setting(
                "Mount position",
                _with_default(
                    self.config.base.mount_xyz,
                    raw_base is not None and "mount_xyz" not in raw_base,
                ),
                indent=2,
            ),
            _setting(
                "Mount orientation RPY",
                orientation,
                indent=2,
            ),
        ])

        raw_integration = None if self._raw is None else self._raw.get("integration", {})
        lines.extend([
            "",
            "Integration",
            _setting(
                "Method",
                _with_default(
                    self.config.integration.method,
                    raw_integration is not None and "method" not in raw_integration,
                ),
                indent=2,
            ),
            _setting(
                "Quadrature order",
                self.config.integration.order
                if self.config.integration.order is not None else "not applicable",
                indent=2,
            ),
            "",
            "Ritz basis",
            _setting("X coefficients", self.config.ritz_x or "not configured", indent=2),
            _setting("Y coefficients", self.config.ritz_y or "not configured", indent=2),
            _setting("Z coefficients", self.config.ritz_z or "not configured", indent=2),
            "",
            "Model parameter inputs",
        ])
        if self.config.parameters:
            for name, value in self.config.parameters.items():
                lines.append(_setting(name, value, indent=2))
        else:
            lines.append("  none")

        lines.extend(["", "Actuation"])
        actuation = self.config.actuation
        if actuation is None:
            lines.append(_setting("Family", "none", indent=2))
        else:
            lines.extend([
                _setting("Family", actuation.family, indent=2),
                _setting("Acceleration mode", actuation.acceleration, indent=2),
                _setting("Channels", len(actuation.channels), indent=2),
            ])
            for name, value in actuation.data.items():
                if name not in {"family", "acceleration", "channels"}:
                    lines.append(_setting(name, value, indent=2))
            for channel in actuation.channels:
                lines.extend([
                    "",
                    f"  Channel {channel.name}",
                    _setting("Kind", channel.kind, indent=4),
                ])
                for span in channel.spans:
                    detail = (
                        f"radius {_format_value(span.radius)}, "
                        f"angle {_format_value(span.angle)} rad"
                    )
                    lines.append(_setting(
                        f"Section {span.section} routing", detail, indent=4
                    ))

        lines.extend(["", "Constraint"])
        constraint = self.config.constraint
        if constraint is None:
            lines.append(_setting("Family", "none", indent=2))
        else:
            lines.append(_setting("Family", constraint.family, indent=2))
            for name, value in constraint.data.items():
                if name != "family":
                    lines.append(_setting(name, value, indent=2))

    def _symbolic_strategy(self, lines: list[str]) -> None:
        _section(lines, "Symbolic Strategy")
        bias = "Wolfram" if self.options.backend == "wolfram" else f"SymPy {sp.__version__}"
        normalization = (
            "Wolfram FactorTerms (enabled)"
            if self.options.wolfram_factor_terms else "disabled"
        )
        cse = "Wolfram experimental CSE" if self.options.wolfram_cse else "SymPy CSE"
        strategy_width = 32
        lines.extend([
            _setting("Model construction", f"SymPy {sp.__version__}", width=strategy_width),
            _setting("Bias differentiation", bias, width=strategy_width),
            _setting("Expression normalization", normalization, width=strategy_width),
            _setting("Common-subexpression elimination", cse, width=strategy_width),
            _setting(
                "Fallback policy",
                "disabled; selected operations stop on failure",
                width=strategy_width,
            ),
        ])
        if self.options.backend == "wolfram":
            lines.extend(["", "Wolfram Kernel"])
            if self.wolfram_executable is not None:
                lines.append(_setting("Executable", self.wolfram_executable, indent=2))
            if self.wolfram_version is not None:
                lines.append(_setting("Version", self.wolfram_version, indent=2))
            lines.append(_setting("Operation timeout", f"{self.options.wolfram_timeout:g} s", indent=2))

    def _build_stages(self, lines: list[str]) -> None:
        _section(lines, "Build Stages")
        width = max((len(item.name) for item in self.timings), default=0)
        for item in self.timings:
            status = "FAIL" if item.failed else "OK"
            suffix = " (failed)" if item.failed else ""
            lines.append(
                f"  {status:<4}  {item.name:<{width}}  {item.seconds:>9.3f} s{suffix}"
            )
        if any(item.name == "Diagnostics analysis" for item in self.timings):
            lines.extend([
                "        Diagnostics analysis is excluded from total build time."
            ])

    def _resolved_model(self, lines: list[str]) -> None:
        assert self.system is not None
        plant = self.plant or self.system.plant
        actuation = self.system.actuation
        constraint = self.system.constraint
        nq = len(plant.q)
        _section(lines, "Resolved Model")
        lines.extend([
            "Coordinates",
            _setting("Base coordinates", len(plant.base_q), indent=2),
            _setting("Arm coordinates", len(plant.arm_q), indent=2),
            _setting("Generalized coordinates", nq, indent=2),
            _setting("State dimension", 2 * nq, indent=2),
            "",
            "Channels",
            _setting("Actuator channels", 0 if actuation is None else actuation.count, indent=2),
            _setting("Constraint channels", 0 if constraint is None else constraint.count, indent=2),
        ])
        plant_parameters = plant.parameters
        actuation_parameters = () if actuation is None else actuation.parameters
        constraint_parameters = () if constraint is None else constraint.parameters
        lines.extend([
            "",
            "Runtime parameters",
            _setting("Plant", len(plant_parameters), indent=2),
            _setting("Actuation", len(actuation_parameters), indent=2),
            _setting("Constraint", len(constraint_parameters), indent=2),
            _setting(
                "Total",
                len(plant_parameters) + len(actuation_parameters) + len(constraint_parameters),
                indent=2,
            ),
            "",
            "Key symbolic outputs",
            _setting("Mass matrix M", _format_shape(plant.mass.shape), indent=2),
            _setting("Bias vector h", _format_shape((nq, 1)), indent=2),
            _setting("Damping matrix D", _format_shape(plant.damping.shape), indent=2),
            _setting(
                "Section kinematics H",
                _format_shape((4, 4, plant.config.segments)),
                indent=2,
            ),
            _setting("End Jacobian Je", _format_shape(plant.end_jacobian.shape), indent=2),
            _setting("Base Jacobian Jb", _format_shape(plant.base_jacobian.shape), indent=2),
            _setting(
                "Vehicle wrench map Bv",
                _format_shape(plant.vehicle_wrench_map.shape),
                indent=2,
            ),
            _setting("Arm force map Ba", _format_shape(plant.arm_force_map.shape), indent=2),
        ])
        if actuation is not None:
            lines.extend([
                _setting("Actuator coordinates", _format_shape(actuation.coordinates.shape), indent=2),
                _setting("Actuator Jacobian Ja", _format_shape(actuation.jacobian.shape), indent=2),
                _setting("Actuator velocity bias", _format_shape(actuation.velocity_bias.shape), indent=2),
            ])
        if constraint is not None:
            lines.extend([
                _setting("Constraint coordinates", _format_shape(constraint.coordinates.shape), indent=2),
                _setting("Constraint Jacobian A", _format_shape(constraint.jacobian.shape), indent=2),
                _setting("Constraint velocity bias", _format_shape(constraint.velocity_bias.shape), indent=2),
                _setting("Constraint reaction map", _format_shape(constraint.reaction_map.shape), indent=2),
            ])
        lines.extend([
            "",
            "Resolved runtime parameters include configured values and model defaults.",
        ])
        self._parameter_group(lines, "Plant parameters", plant_parameters)
        if actuation_parameters:
            self._parameter_group(lines, "Actuation parameters", actuation_parameters)
        if constraint_parameters:
            self._parameter_group(lines, "Constraint parameters", constraint_parameters)

    def _parameter_group(self, lines: list[str], title: str, parameters: tuple[Any, ...]) -> None:
        lines.extend(["", title])
        width = max((len(item.name) for item in parameters), default=0)
        groups: list[list[Any]] = []
        if title == "Plant parameters":
            global_parameters = [
                item for item in parameters
                if re.fullmatch(r"s\d+_.+", item.name) is None
            ]
            if global_parameters:
                groups.append(global_parameters)
            section_numbers = sorted({
                int(match.group(1))
                for item in parameters
                if (match := re.fullmatch(r"s(\d+)_.+", item.name)) is not None
            })
            groups.extend([
                [
                    item for item in parameters
                    if re.fullmatch(rf"s{section}_.+", item.name) is not None
                ]
                for section in section_numbers
            ])
        else:
            groups.append(list(parameters))
        for group_index, group in enumerate(groups):
            if group_index:
                lines.append("")
            for item in group:
                lines.append(f"  {item.name:<{width}} = {_format_value(item.default)}")

    def _exported_functions(self, lines: list[str]) -> None:
        _section(lines, "Exported Symbolic Functions")
        lines.extend([
            "Operation counts are symbolic estimates, not measured hardware FLOPs.",
            '"Generated ops" includes CSE temporary definitions and final output expressions.',
            "",
        ])
        name_width = max(32, max(len(item.name) for item in self.function_metrics))
        shape_width = max(10, max(len(_format_shape(item.shape)) for item in self.function_metrics))
        lines.append(
            f"{'Function':<{name_width}}  {'Output':<{shape_width}}  "
            f"{'Raw ops':>12}  {'CSE temps':>10}  {'Generated ops':>14}  "
            f"{'Reduction':>10}  {'MATLAB size':>11}"
        )
        for item in self.function_metrics:
            reduction = (
                "-" if item.raw_operations == 0
                else f"{100 * (item.raw_operations - item.generated_operations) / item.raw_operations:.1f}%"
            )
            lines.append(
                f"{item.name:<{name_width}}  {_format_shape(item.shape):<{shape_width}}  "
                f"{item.raw_operations:>12,}  {item.temporary_count:>10,}  "
                f"{item.generated_operations:>14,}  {reduction:>10}  "
                f"{_format_source_size(item.source_bytes):>11}"
            )

    def _generated_artifacts(self, lines: list[str]) -> None:
        assert self.artifact_metrics is not None
        metrics = self.artifact_metrics
        _section(lines, "Generated Artifacts")
        lines.extend(["Bundle directory", f"  {self.output}", "", "Contents"])
        labels = ((".m", "MATLAB files"), (".json", "JSON files"), (".tex", "TeX files"))
        for suffix, label in labels:
            count = metrics.counts.get(suffix, 0)
            size = metrics.sizes.get(suffix, 0)
            noun = "file" if count == 1 else "files"
            lines.append(_setting(label, f"{count} {noun}, {size:,} bytes", indent=2))
        total_size = sum(path.stat().st_size for path in metrics.files)
        total_noun = "file" if len(metrics.files) == 1 else "files"
        lines.append(_setting(
            "Total",
            f"{len(metrics.files)} {total_noun}, {total_size:,} bytes",
            indent=2,
        ))
        names = {path.name for path in metrics.files}
        lines.extend([
            "",
            "Optional outputs",
            _setting(
                "Centerline functions",
                "generated" if "softarm_centerline.m" in names else "not generated",
                indent=2,
            ),
            _setting(
                "Actuator functions",
                "generated" if any(name.startswith("softarm_actuator_") for name in names)
                else "not applicable",
                indent=2,
            ),
            _setting(
                "Constraint functions",
                "generated" if any(name.startswith("softarm_constraint_") for name in names)
                else "not applicable",
                indent=2,
            ),
            _setting(
                "Exact TeX appendix",
                "generated" if self.options.tex_appendix else "not requested",
                indent=2,
            ),
        ])

    def _build_result(self, lines: list[str], success: bool) -> None:
        _section(lines, "Build Result")
        lines.append("SUCCESS" if success else "FAILED")
        if success:
            lines.extend(["", "Generated MATLAB bundle:", f"  {self.output}"])
