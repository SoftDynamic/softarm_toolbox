from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass
from typing import Protocol

import sympy as sp


class ExpressionNormalizer(Protocol):
    def factor_terms(self, expressions: Sequence[sp.Expr]) -> list[sp.Expr]: ...


class CommonSubexpressionEliminator(Protocol):
    def cse(
        self,
        expressions: Sequence[sp.Expr],
        prefix: str = "t",
        order: str = "none",
    ) -> tuple[list[tuple[sp.Symbol, sp.Expr]], list[sp.Expr]]: ...


class SympyCse:
    def cse(
        self,
        expressions: Sequence[sp.Expr],
        prefix: str = "t",
        order: str = "none",
    ) -> tuple[list[tuple[sp.Symbol, sp.Expr]], list[sp.Expr]]:
        replacements, reduced = sp.cse(
            expressions,
            symbols=sp.numbered_symbols(prefix),
            order=order,
        )
        return list(replacements), list(reduced)


@dataclass(frozen=True)
class OptimizedFunction:
    replacements: tuple[tuple[sp.Symbol, sp.Expr], ...]
    expressions: tuple[sp.Expr, ...]
    shape: tuple[int, ...]


class FunctionOptimizer:
    """Function-boundary optimization with an in-memory build cache."""

    def __init__(
        self,
        *,
        normalizer: ExpressionNormalizer | None = None,
        eliminator: CommonSubexpressionEliminator | None = None,
    ):
        self._normalizer = normalizer
        self._eliminator = eliminator or SympyCse()
        self._cache: dict[
            tuple[tuple[sp.Expr, ...], tuple[int, ...]], OptimizedFunction
        ] = {}

    def optimize(
        self,
        expressions: Sequence[sp.Expr],
        shape: tuple[int, ...],
    ) -> OptimizedFunction:
        canonical = tuple(expressions)
        key = canonical, shape
        cached = self._cache.get(key)
        if cached is not None:
            return cached
        normalized = (
            list(canonical)
            if self._normalizer is None
            else self._normalizer.factor_terms(canonical)
        )
        if len(normalized) != len(canonical):
            raise ValueError(
                "expression normalizer changed the function output dimension"
            )
        replacements, reduced = self._eliminator.cse(
            normalized, prefix="t", order="none"
        )
        if len(reduced) != len(canonical):
            raise ValueError("CSE changed the function output dimension")
        result = OptimizedFunction(
            tuple(replacements), tuple(reduced), shape
        )
        self._cache[key] = result
        return result
