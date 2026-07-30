from __future__ import annotations

import sympy as sp


def optimize(expressions: list[sp.Expr]) -> list[sp.Expr]:
    # Derivation already preserves useful matrix structure. Global factoring is
    # disproportionately expensive for multi-section Coriolis expressions; CSE
    # is performed once by the target emitter instead.
    return expressions
