from __future__ import annotations

from dataclasses import dataclass, replace
from pathlib import Path
from typing import Literal

from .actuation import ActuationModel, derive_actuation
from .backends.wolfram import WolframError, WolframKernel
from .codegen.matlab import generate_matlab_bundle
from .codegen.optimization import FunctionOptimizer, SympyCse
from .config import ModelConfig
from .constraints import ConstraintModel, derive_constraint
from .derive import SymbolicPlant, derive
from .dynamics import SympyBatchDifferentiator, assemble_bias


@dataclass(frozen=True)
class BuildOptions:
    backend: Literal["sympy", "wolfram"] = "sympy"
    wolfram_kernel: str | None = None
    wolfram_timeout: float = 600.0
    wolfram_cse: bool = False
    wolfram_factor_terms: bool = False
    tex_appendix: bool = False

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
    plant: SymbolicPlant
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


def derive_system(config: ModelConfig) -> DerivedSystem:
    plant = derive(config)
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
) -> Path:
    differentiator = kernel or SympyBatchDifferentiator()
    try:
        bias = assemble_bias(system.plant, differentiator)
    except Exception as error:
        strategy = "Wolfram" if kernel is not None else "SymPy"
        raise BuildError("bias differentiation", strategy, error) from error

    plant = replace(system.plant, _bias=bias)
    normalizer = kernel if options.wolfram_factor_terms else None
    eliminator = kernel if options.wolfram_cse else SympyCse()
    optimizer = FunctionOptimizer(
        normalizer=normalizer,
        eliminator=eliminator,
    )
    try:
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
    try:
        system = derive_system(config)
    except Exception as error:
        raise BuildError("model derivation", "SymPy", error) from error

    if options.backend == "sympy":
        return _materialize_and_generate(system, output, options, None)

    try:
        with WolframKernel(
            options.wolfram_kernel,
            timeout=options.wolfram_timeout,
        ) as kernel:
            return _materialize_and_generate(system, output, options, kernel)
    except BuildError:
        raise
    except WolframError as error:
        raise BuildError("Wolfram startup", "Wolfram", error) from error
