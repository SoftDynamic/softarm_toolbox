from __future__ import annotations

from collections.abc import Iterable
from typing import Any

import sympy as sp

from ..special import (
    AffineCosMoment,
    AffineSinMoment,
    CoscSqrt,
    CoscSqrtD,
    CoscSqrtDD,
    Sinc3Sqrt,
    Sinc3SqrtD,
    Sinc3SqrtDD,
    SincSqrt,
    SincSqrtD,
    SincSqrtDD,
)

_HEADS = {
    sp.Add: "add",
    sp.Mul: "mul",
    sp.Pow: "pow",
    sp.sin: "sin",
    sp.cos: "cos",
    sp.Abs: "abs",
}

_CONSTRUCTORS = {
    "add": sp.Add,
    "mul": sp.Mul,
    "pow": sp.Pow,
    "sin": sp.sin,
    "cos": sp.cos,
    "abs": sp.Abs,
    "SincSqrt": SincSqrt,
    "SincSqrtD": SincSqrtD,
    "SincSqrtDD": SincSqrtDD,
    "CoscSqrt": CoscSqrt,
    "CoscSqrtD": CoscSqrtD,
    "CoscSqrtDD": CoscSqrtDD,
    "Sinc3Sqrt": Sinc3Sqrt,
    "Sinc3SqrtD": Sinc3SqrtD,
    "Sinc3SqrtDD": Sinc3SqrtDD,
    "AffineCosMoment": AffineCosMoment,
    "AffineSinMoment": AffineSinMoment,
}


def _decode_real(value: str) -> sp.Float:
    # Wolfram InputForm appends arbitrary-precision marks such as
    # 1.23`28.7 and uses *^ for scientific notation.
    return sp.Float(value.split("`", 1)[0].replace("*^", "e"))


def _symbol_assumptions(symbol: sp.Symbol) -> dict[str, bool]:
    """Return only assumptions that affect the supported CAS operations."""
    return {
        name: True
        for name in ("real", "positive", "nonnegative")
        if symbol.assumptions0.get(name) is True
    }


def encode_dag(expressions: Iterable[sp.Expr]) -> dict[str, Any]:
    """Encode expressions as a post-order DAG instead of repeated trees."""
    nodes: list[dict[str, Any]] = []
    indices: dict[sp.Basic, int] = {}

    def visit(expr: sp.Basic) -> int:
        cached = indices.get(expr)
        if cached is not None:
            return cached
        args = [visit(arg) for arg in expr.args]
        if isinstance(expr, sp.Symbol):
            node: dict[str, Any] = {
                "head": "symbol",
                "name": str(expr),
                "assumptions": _symbol_assumptions(expr),
            }
        elif isinstance(expr, sp.Integer):
            node = {"head": "integer", "value": str(expr)}
        elif isinstance(expr, sp.Rational):
            node = {
                "head": "rational",
                "numerator": str(expr.p),
                "denominator": str(expr.q),
            }
        elif isinstance(expr, sp.Float):
            node = {"head": "real", "value": str(expr)}
        else:
            node = {
                "head": _HEADS.get(expr.func, expr.func.__name__),
                "args": args,
            }
        index = len(nodes)
        nodes.append(node)
        indices[expr] = index
        return index

    roots = [visit(expr) for expr in expressions]
    return {"nodes": nodes, "roots": roots}


def decode_dag(
    graph: dict[str, Any],
    symbols: dict[str, sp.Symbol] | None = None,
) -> list[sp.Expr]:
    """Decode a post-order DAG, reusing canonical input symbols by name."""
    symbol_table = {} if symbols is None else dict(symbols)
    values: list[sp.Expr] = []
    for node in graph["nodes"]:
        head = node["head"]
        if head == "symbol":
            name = node["name"]
            value = symbol_table.get(name)
            if value is None:
                assumptions = {
                    key: bool(item)
                    for key, item in node.get("assumptions", {}).items()
                    if item
                }
                value = sp.Symbol(name, **assumptions)
                symbol_table[name] = value
        elif head == "integer":
            value = sp.Integer(node["value"])
        elif head == "rational":
            value = sp.Rational(node["numerator"], node["denominator"])
        elif head == "real":
            value = _decode_real(node["value"])
        else:
            args = [values[index] for index in node.get("args", [])]
            value = (_CONSTRUCTORS.get(head) or sp.Function(head))(*args)
        values.append(value)
    return [values[index] for index in graph["roots"]]
