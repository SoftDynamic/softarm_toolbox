from __future__ import annotations

from typing import Any, Iterable
import sympy as sp

from .special import (
    CoscSqrt, CoscSqrtD, CoscSqrtDD,
    Sinc3Sqrt, Sinc3SqrtD, Sinc3SqrtDD,
    SincSqrt, SincSqrtD, SincSqrtDD,
)


_HEADS = {
    sp.Add: "add",
    sp.Mul: "mul",
    sp.Pow: "pow",
    sp.sin: "sin",
    sp.cos: "cos",
    sp.Abs: "abs",
}


def _decode_real(value: str) -> sp.Float:
    # Wolfram InputForm appends arbitrary-precision marks such as
    # 1.23`28.7 and uses *^ for scientific notation.
    return sp.Float(value.split("`", 1)[0].replace("*^", "e"))


def encode(expr: sp.Expr) -> dict[str, Any]:
    if isinstance(expr, sp.Symbol):
        return {"head": "symbol", "name": str(expr)}
    if isinstance(expr, sp.Integer):
        return {"head": "integer", "value": str(expr)}
    if isinstance(expr, sp.Rational):
        return {"head": "rational", "numerator": str(expr.p), "denominator": str(expr.q)}
    if isinstance(expr, sp.Float):
        return {"head": "real", "value": str(expr)}
    head = _HEADS.get(expr.func)
    if head is None:
        head = expr.func.__name__
    return {"head": head, "args": [encode(arg) for arg in expr.args]}


def decode(node: dict[str, Any]) -> sp.Expr:
    if not isinstance(node, dict):
        raise TypeError(f"AST node must be an object, got {node!r}")
    head = node["head"]
    if head == "symbol":
        return sp.Symbol(node["name"], real=True)
    if head == "integer":
        return sp.Integer(node["value"])
    if head == "rational":
        return sp.Rational(node["numerator"], node["denominator"])
    if head == "real":
        return _decode_real(node["value"])
    try:
        args = [decode(arg) for arg in node.get("args", [])]
    except TypeError as error:
        raise TypeError(f"invalid child of {head}: {error}") from error
    constructors = {
        "add": sp.Add, "mul": sp.Mul, "pow": sp.Pow, "sin": sp.sin, "cos": sp.cos, "abs": sp.Abs,
        "SincSqrt": SincSqrt, "SincSqrtD": SincSqrtD, "SincSqrtDD": SincSqrtDD,
        "CoscSqrt": CoscSqrt, "CoscSqrtD": CoscSqrtD, "CoscSqrtDD": CoscSqrtDD,
        "Sinc3Sqrt": Sinc3Sqrt, "Sinc3SqrtD": Sinc3SqrtD, "Sinc3SqrtDD": Sinc3SqrtDD,
    }
    constructor = constructors.get(head) or sp.Function(head)
    return constructor(*args)


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
    constructors = {
        "add": sp.Add, "mul": sp.Mul, "pow": sp.Pow,
        "sin": sp.sin, "cos": sp.cos, "abs": sp.Abs,
        "SincSqrt": SincSqrt, "SincSqrtD": SincSqrtD,
        "SincSqrtDD": SincSqrtDD, "CoscSqrt": CoscSqrt,
        "CoscSqrtD": CoscSqrtD, "CoscSqrtDD": CoscSqrtDD,
        "Sinc3Sqrt": Sinc3Sqrt, "Sinc3SqrtD": Sinc3SqrtD,
        "Sinc3SqrtDD": Sinc3SqrtDD,
    }
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
            value = (constructors.get(head) or sp.Function(head))(*args)
        values.append(value)
    return [values[index] for index in graph["roots"]]
