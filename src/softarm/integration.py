from __future__ import annotations

import sympy as sp
from sympy.integrals.quadrature import gauss_legendre

from .config import IntegrationConfig


class IntegrationError(RuntimeError):
    pass


def integrate_unit(expr: sp.Expr | sp.Matrix, xi: sp.Symbol, config: IntegrationConfig, label: str):
    if config.method == "gauss":
        nodes, weights = gauss_legendre(config.order or 1, 30)
        result = expr * 0
        for node, weight in zip(nodes, weights, strict=True):
            mapped = (node + 1) / 2
            result += sp.N(weight / 2, 17) * expr.subs(xi, mapped)
        return result

    result = expr.applyfunc(lambda item: sp.integrate(item, (xi, 0, 1))) if isinstance(expr, sp.MatrixBase) else sp.integrate(expr, (xi, 0, 1))
    values = list(result) if isinstance(result, sp.MatrixBase) else [result]
    if any(value.has(sp.Integral) for value in values):
        raise IntegrationError(
            f"analytic integration failed for {label}; explicitly select integration.method='gauss' and an order"
        )
    return result

