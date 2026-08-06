from __future__ import annotations

import time
from collections.abc import Iterator
from contextlib import contextmanager
from dataclasses import dataclass, replace
from pathlib import Path
from typing import Literal

from .actuation import ActuationModel, derive_actuation
from .backends.wolfram import WolframError, WolframKernel
from .codegen.matlab import generate_matlab_bundle
from .codegen.optimization import FunctionOptimizer, SympyCse
from .config import ModelConfig
from .constraints import ConstraintModel, derive_constraint
from .derive import derive
from .diagnostics import BuildDiagnostics
from .dynamics import SympyBatchDifferentiator, assemble_bias
from .models import PlantModel, RecursivePlant, SymbolicPlant


@dataclass(frozen=True)
class BuildOptions:
    backend: Literal["sympy", "wolfram"] = "sympy"
    wolfram_kernel: str | None = None
    wolfram_timeout: float = 600.0
    wolfram_cse: bool = False
    wolfram_factor_terms: bool = False
    tex_appendix: bool = False
    verbose: bool = False

    def __post_init__(self) -> None:
        if self.backend not in {"sympy", "wolfram"}:
            raise ValueError(f"unknown symbolic backend: {self.backend}")
        if self.wolfram_timeout <= 0:
            raise ValueError("Wolfram timeout must be positive")
        wolfram_only = (
            self.wolfram_kernel is not None
            or self.wolfram_cse
            or self.wolfram_factor_terms
            or self.wolfram_timeout != 600.0
        )
        if self.backend != "wolfram" and wolfram_only:
            raise ValueError(
                "Wolfram options require backend='wolfram'"
            )


@dataclass(frozen=True)
class DerivedSystem:
    plant: PlantModel
    actuation: ActuationModel | None
    constraint: ConstraintModel | None


DEFAULT_BUILD_OPTIONS = BuildOptions()


class BuildError(RuntimeError):
    def __init__(self, stage: str, strategy: str, cause: Exception):
        self.stage = stage
        self.strategy = strategy
        self.cause = cause
        detail = f"{stage} [{strategy}] failed: {cause}"
        lowered = str(cause).lower()
        if "factor_terms" in lowered or "factorterms" in lowered:
            detail += (
                "\nNo alternative strategy was attempted. Retry without "
                "--wolfram-factor-terms or increase --wolfram-timeout."
            )
        elif "cse" in lowered or "optimized-expression" in lowered:
            detail += (
                "\nNo alternative strategy was attempted. Retry without "
                "--wolfram-cse to use the default SymPy CSE."
            )
        elif strategy == "Wolfram":
            detail += "\nNo SymPy fallback was attempted."
        super().__init__(detail)


@contextmanager
def _timed_stage(
    diagnostics: BuildDiagnostics | None,
    stage: str,
) -> Iterator[None]:
    if diagnostics is None:
        yield
        return
    start = time.perf_counter()
    try:
        yield
    except BaseException:
        elapsed = time.perf_counter() - start
        diagnostics.record_timing(stage, elapsed, True)
        raise
    elapsed = time.perf_counter() - start
    diagnostics.record_timing(stage, elapsed, False)


def derive_system(config: ModelConfig) -> DerivedSystem:
    plant = derive(config)
    if isinstance(plant, RecursivePlant) and config.constraint is not None:
        raise ValueError("recursive dynamics do not yet support constraints")
    return DerivedSystem(
        plant=plant,
        actuation=derive_actuation(plant),
        constraint=derive_constraint(plant),
    )


def symbolic_plan(options: BuildOptions) -> str:
    bias = "Wolfram" if options.backend == "wolfram" else "SymPy"
    factor_terms = "on" if options.wolfram_factor_terms else "off"
    cse = "Wolfram (experimental)" if options.wolfram_cse else "SymPy"
    return (
        "Symbolic plan: derive=SymPy, "
        f"bias={bias}, FactorTerms={factor_terms}, CSE={cse}"
    )


def _materialize_and_generate(
    system: DerivedSystem,
    output: str | Path,
    options: BuildOptions,
    kernel: WolframKernel | None,
    diagnostics: BuildDiagnostics | None,
) -> Path:
    if isinstance(system.plant, RecursivePlant):
        optimizer = FunctionOptimizer(
            eliminator=SympyCse(),
            collect_diagnostics=diagnostics is not None,
        )
        if diagnostics is not None:
            diagnostics.capture_plant(system.plant)
            diagnostics.capture_optimizer(optimizer)
        try:
            with _timed_stage(diagnostics, "MATLAB recursive-kernel generation/CSE"):
                return generate_matlab_bundle(
                    system.plant,
                    output,
                    actuation=system.actuation,
                    constraint=None,
                    tex_appendix=options.tex_appendix,
                    optimizer=optimizer,
                )
        except Exception as error:
            raise BuildError("code generation", "recursive + SymPy CSE", error) from error
    if not isinstance(system.plant, SymbolicPlant):
        raise TypeError("unsupported plant representation")
    differentiator = kernel or SympyBatchDifferentiator()
    strategy = "Wolfram" if kernel is not None else "SymPy"
    try:
        with _timed_stage(diagnostics, f"Bias differentiation [{strategy}]"):
            bias = assemble_bias(system.plant, differentiator)
    except Exception as error:
        raise BuildError("bias differentiation", strategy, error) from error

    plant = replace(system.plant, _bias=bias)
    if diagnostics is not None:
        diagnostics.capture_plant(plant)
    normalizer = kernel if options.wolfram_factor_terms else None
    eliminator = kernel if options.wolfram_cse else SympyCse()
    optimizer = FunctionOptimizer(
        normalizer=normalizer,
        eliminator=eliminator,
        collect_diagnostics=diagnostics is not None,
    )
    if diagnostics is not None:
        diagnostics.capture_optimizer(optimizer)
    try:
        with _timed_stage(diagnostics, "MATLAB generation/CSE"):
            return generate_matlab_bundle(
                plant,
                output,
                actuation=system.actuation,
                constraint=system.constraint,
                tex_appendix=options.tex_appendix,
                optimizer=optimizer,
            )
    except Exception as error:
        steps = []
        if options.wolfram_factor_terms:
            steps.append("Wolfram FactorTerms")
        steps.append(
            "Wolfram experimental CSE"
            if options.wolfram_cse
            else "SymPy CSE"
        )
        raise BuildError("code generation", " + ".join(steps), error) from error


def build_bundle(
    config: ModelConfig,
    output: str | Path,
    options: BuildOptions = DEFAULT_BUILD_OPTIONS,
) -> Path:
    diagnostics = (
        BuildDiagnostics(config, output, options) if options.verbose else None
    )
    try:
        with _timed_stage(diagnostics, "Total build"):
            result = _build_bundle(config, output, options, diagnostics)
    except BaseException:
        if diagnostics is not None:
            print(diagnostics.render(success=False), end="")
        raise
    if diagnostics is not None:
        try:
            with _timed_stage(diagnostics, "Diagnostics analysis"):
                diagnostics.analyze()
        except BaseException:
            print(diagnostics.render(success=False), end="")
            raise
        print(diagnostics.render(success=True), end="")
    return result


def _build_bundle(
    config: ModelConfig,
    output: str | Path,
    options: BuildOptions,
    diagnostics: BuildDiagnostics | None,
) -> Path:
    if config.dynamics.formulation == "recursive":
        if options.backend != "sympy":
            raise ValueError("recursive dynamics require backend='sympy'")
        if options.tex_appendix:
            raise ValueError("recursive dynamics do not support --tex-appendix")
    try:
        with _timed_stage(diagnostics, "Model derivation [SymPy]"):
            system = derive_system(config)
    except Exception as error:
        raise BuildError("model derivation", "SymPy", error) from error
    if diagnostics is not None:
        diagnostics.capture_system(system)

    if options.backend == "sympy":
        return _materialize_and_generate(
            system, output, options, None, diagnostics
        )

    try:
        with _timed_stage(diagnostics, "Wolfram Kernel startup"):
            kernel_context = WolframKernel(
                options.wolfram_kernel,
                timeout=options.wolfram_timeout,
            )
        if diagnostics is not None:
            diagnostics.capture_wolfram(kernel_context)
        with kernel_context as kernel:
            if diagnostics is not None:
                diagnostics.capture_wolfram(kernel)
            return _materialize_and_generate(
                system, output, options, kernel, diagnostics
            )
    except BuildError:
        raise
    except WolframError as error:
        raise BuildError("Wolfram startup", "Wolfram", error) from error
