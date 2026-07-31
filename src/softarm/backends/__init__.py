from __future__ import annotations

from .sympy_backend import optimize as optimize_sympy
from .wolfram_backend import optimize as optimize_wolfram
from .session import (
    SymbolicBackendError,
    SymbolicSession,
    WolframSession,
    create_session,
)

__all__ = [
    "SymbolicBackendError",
    "SymbolicSession",
    "WolframSession",
    "create_session",
    "optimize",
]


def optimize(expressions, backend: str, wolfram_kernel: str | None = None):
    if backend == "sympy":
        return optimize_sympy(expressions)
    if backend == "wolfram":
        return optimize_wolfram(expressions, wolfram_kernel)
    raise ValueError(f"unknown CAS backend: {backend}")
