from __future__ import annotations

from collections.abc import Sequence
from typing import TYPE_CHECKING, Protocol

import sympy as sp

if TYPE_CHECKING:
    from .derive import SymbolicPlant


class BatchDifferentiator(Protocol):
    """Narrow symbolic capability required by dynamics assembly."""

    def differentiate(
        self,
        expressions: Sequence[sp.Expr],
        variables: Sequence[sp.Symbol],
    ) -> list[sp.Expr]: ...


class SympyBatchDifferentiator:
    """Reference implementation using the canonical SymPy expressions."""

    def differentiate(
        self,
        expressions: Sequence[sp.Expr],
        variables: Sequence[sp.Symbol],
    ) -> list[sp.Expr]:
        return [
            sp.diff(expression, variable)
            for expression in expressions
            for variable in variables
        ]


def assemble_bias(
    plant: SymbolicPlant,
    differentiator: BatchDifferentiator,
) -> sp.Matrix:
    """Assemble Coriolis, conservative, and damping terms from one formula."""
    nq = len(plant.q)
    upper = [
        (row, column)
        for row in range(nq)
        for column in range(row, nq)
    ]
    derivatives = differentiator.differentiate(
        [plant.mass[row, column] for row, column in upper] + [plant.potential],
        list(plant.q),
    )
    expected = (len(upper) + 1) * nq
    if len(derivatives) != expected:
        raise ValueError(
            "batch differentiator returned "
            f"{len(derivatives)} values; expected {expected}"
        )

    mass_derivatives: dict[tuple[int, int, int], sp.Expr] = {}
    for pair_index, (row, column) in enumerate(upper):
        for coordinate in range(nq):
            value = derivatives[pair_index * nq + coordinate]
            mass_derivatives[coordinate, row, column] = value
            mass_derivatives[coordinate, column, row] = value

    coriolis = sp.Matrix([
        sp.Add(*(
            mass_derivatives[k, i, j] * plant.dq[j] * plant.dq[k]
            - sp.Rational(1, 2)
            * mass_derivatives[i, j, k]
            * plant.dq[j]
            * plant.dq[k]
            for j in range(nq)
            for k in range(nq)
        ))
        for i in range(nq)
    ])
    conservative = sp.Matrix(derivatives[-nq:])
    return coriolis + conservative + plant.damping * plant.dq
