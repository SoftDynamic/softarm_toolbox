from __future__ import annotations

from typing import Any
import sympy as sp

from .special import CoscSqrt, CoscSqrtD, CoscSqrtDD, SincSqrt, SincSqrtD, SincSqrtDD


_HEADS = {
    sp.Add: "add",
    sp.Mul: "mul",
    sp.Pow: "pow",
    sp.sin: "sin",
    sp.cos: "cos",
    sp.Abs: "abs",
}


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
        return sp.Float(node["value"])
    try:
        args = [decode(arg) for arg in node.get("args", [])]
    except TypeError as error:
        raise TypeError(f"invalid child of {head}: {error}") from error
    constructors = {
        "add": sp.Add, "mul": sp.Mul, "pow": sp.Pow, "sin": sp.sin, "cos": sp.cos, "abs": sp.Abs,
        "SincSqrt": SincSqrt, "SincSqrtD": SincSqrtD, "SincSqrtDD": SincSqrtDD,
        "CoscSqrt": CoscSqrt, "CoscSqrtD": CoscSqrtD, "CoscSqrtDD": CoscSqrtDD,
    }
    constructor = constructors.get(head) or sp.Function(head)
    return constructor(*args)
