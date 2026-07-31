from __future__ import annotations

import sympy as sp
from sympy.integrals.quadrature import gauss_legendre

from .backends.session import SymbolicSession
from .config import IntegrationConfig


class IntegrationError(RuntimeError):
    pass


def unit_gauss_rule(config: IntegrationConfig) -> list[tuple[sp.Expr, sp.Expr]]:
    nodes, weights = gauss_legendre(config.order or 1, 30)
    return [
        ((node + 1) / 2, sp.N(weight / 2, 17))
        for node, weight in zip(nodes, weights, strict=True)
    ]


def integrate_unit(
    expr: sp.Expr | sp.Matrix,
    xi: sp.Symbol,
    config: IntegrationConfig,
    label: str,
    symbolic: SymbolicSession | None = None,
):
    executor = symbolic or SymbolicSession()
    if config.method == "gauss":
        result = expr * 0
        for node, weight in unit_gauss_rule(config):
            result += weight * executor.substitute(expr, {xi: node})
        return result

    result = executor.integrate(expr, xi, sp.S.Zero, sp.S.One)
    values = list(result) if isinstance(result, sp.MatrixBase) else [result]
    if any(value.has(sp.Integral) for value in values):
        raise IntegrationError(
            f"analytic integration failed for {label}; explicitly select integration.method='gauss' and an order"
        )
    return result
