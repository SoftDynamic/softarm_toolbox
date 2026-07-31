from __future__ import annotations

from contextlib import AbstractContextManager
import hashlib
import json
from pathlib import Path
import threading
import time
from typing import Iterable, Sequence

import sympy as sp

from ..ast import decode_dag, encode_dag
from ..special import SPECIAL_DERIVATIVE_HEADS
from .wolfram_backend import _kernel_path


class SymbolicBackendError(RuntimeError):
    pass


def _flatten(matrix: sp.MatrixBase) -> list[sp.Expr]:
    return [
        matrix[row, column]
        for row in range(matrix.rows)
        for column in range(matrix.cols)
    ]


def _symbols(expressions: Iterable[sp.Expr]) -> dict[str, sp.Symbol]:
    result: dict[str, sp.Symbol] = {}
    for expression in expressions:
        for symbol in expression.free_symbols:
            existing = result.get(str(symbol))
            if existing is not None and existing != symbol:
                raise SymbolicBackendError(
                    f"inconsistent assumptions for symbol {symbol}"
                )
            result[str(symbol)] = symbol
    return result


class SymbolicSession(AbstractContextManager["SymbolicSession"]):
    backend = "sympy"
    version = sp.__version__

    def close(self) -> None:
        pass

    def __exit__(self, exc_type, exc, traceback) -> None:
        self.close()

    @property
    def cache_identity(self) -> tuple[str, str]:
        return self.backend, self.version

    def optimize(self, expressions: Sequence[sp.Expr]) -> list[sp.Expr]:
        return list(expressions)

    def optimize_matrix(self, matrix: sp.MatrixBase) -> sp.Matrix:
        return sp.Matrix(matrix.rows, matrix.cols, self.optimize(_flatten(matrix)))

    def diff(self, expression: sp.Expr, variable: sp.Symbol) -> sp.Expr:
        return self.differentiate([expression], [variable])[0]

    def diff_matrix(
        self, matrix: sp.MatrixBase, variable: sp.Symbol
    ) -> sp.Matrix:
        return sp.Matrix(
            matrix.rows,
            matrix.cols,
            self.differentiate(_flatten(matrix), [variable]),
        )

    def differentiate(
        self, expressions: Sequence[sp.Expr], variables: Sequence[sp.Symbol]
    ) -> list[sp.Expr]:
        return [
            sp.diff(expression, variable)
            for expression in expressions
            for variable in variables
        ]

    def jacobian(
        self, matrix: sp.MatrixBase, variables: Sequence[sp.Symbol]
    ) -> sp.Matrix:
        if not variables:
            return sp.zeros(matrix.rows * matrix.cols, 0)
        return sp.Matrix(_flatten(matrix)).jacobian(list(variables))

    def integrate(
        self,
        expression: sp.Expr | sp.MatrixBase,
        variable: sp.Symbol,
        lower: sp.Expr,
        upper: sp.Expr,
    ) -> sp.Expr | sp.Matrix:
        if isinstance(expression, sp.MatrixBase):
            values = [
                sp.integrate(item, (variable, lower, upper))
                for item in _flatten(expression)
            ]
            return sp.Matrix(expression.rows, expression.cols, values)
        return sp.integrate(expression, (variable, lower, upper))

    def substitute(
        self,
        expression: sp.Expr | sp.MatrixBase,
        substitutions: dict[sp.Symbol, sp.Expr],
    ) -> sp.Expr | sp.Matrix:
        return expression.subs(substitutions)

    def cse(
        self,
        expressions: Sequence[sp.Expr],
        prefix: str = "t",
        order: str = "none",
    ) -> tuple[list[tuple[sp.Symbol, sp.Expr]], list[sp.Expr]]:
        replacements, reduced = sp.cse(
            expressions, symbols=sp.numbered_symbols(prefix), order=order
        )
        return list(replacements), list(reduced)

    def rank(self, matrix: sp.MatrixBase) -> int:
        return int(matrix.rank())


class WolframSession(SymbolicSession):
    backend = "wolfram"
    _SCHEMA = 2

    def __init__(self, kernel: str | None = None, timeout: float = 600.0):
        executable = _kernel_path(kernel)
        if executable.name.lower() in {"math", "math.exe"}:
            session_executable = executable.with_name("WolframKernel.exe")
            if session_executable.is_file():
                executable = session_executable
        bridge = Path(__file__).with_name("wolfram_bridge.wls")
        self.timeout = timeout
        self._lock = threading.Lock()
        self._cache: dict[str, dict] = {}
        try:
            from wolframclient.evaluation import WolframLanguageSession
            from wolframclient.language import wl, wlexpr
        except ImportError as error:
            raise SymbolicBackendError(
                "persistent Wolfram sessions require the 'wolfram' optional "
                "dependency (wolframclient)"
            ) from error
        self._wl = wl
        last_error: Exception | None = None
        for attempt in range(2):
            self._session = WolframLanguageSession(str(executable))
            try:
                self._session.start()
                source = str(bridge).replace("\\", "/").replace('"', '\\"')
                self._session.evaluate(
                    wlexpr(f'$SoftArmLibraryMode=True; Get["{source}"]')
                )
                self.version = str(self._session.evaluate(wlexpr("$Version")))
                break
            except Exception as error:
                last_error = error
                self.close()
                if attempt == 0:
                    time.sleep(0.5)
        else:
            raise SymbolicBackendError(
                f"Wolfram session failed to start: {last_error}"
            ) from last_error

    def _request(
        self,
        operation: str,
        expressions: Sequence[sp.Expr],
        **parameters,
    ) -> list[sp.Expr] | dict:
        canonical = list(expressions)
        calculus_parameters = {
            "variable", "variables", "values", "lower", "upper"
        }
        normalized_parameters = {
            key: (
                [sp.sympify(item) for item in value]
                if key in calculus_parameters
                and isinstance(value, (list, tuple))
                else sp.sympify(value)
                if key in calculus_parameters
                else value
            )
            for key, value in parameters.items()
        }
        symbol_table = _symbols(canonical)
        for value in normalized_parameters.values():
            if isinstance(value, sp.Basic):
                symbol_table.update(_symbols([value]))
            elif isinstance(value, (list, tuple)):
                symbol_table.update(_symbols([
                    item for item in value if isinstance(item, sp.Basic)
                ]))
        payload = {
            "schema": self._SCHEMA,
            "operation": operation,
            "expressions": encode_dag(canonical),
            "parameters": {
                key: (
                    encode_dag([value])
                    if isinstance(value, sp.Basic)
                    else encode_dag(value)
                    if isinstance(value, (list, tuple))
                    and all(isinstance(item, sp.Basic) for item in value)
                    else value
                )
                for key, value in normalized_parameters.items()
            },
            "derivatives": [
                {"head": head, "derivative": derivative}
                for head, derivative in SPECIAL_DERIVATIVE_HEADS
            ],
            "timeout": self.timeout,
        }
        serialized = json.dumps(payload, separators=(",", ":"), sort_keys=True)
        key = hashlib.sha256(serialized.encode("utf-8")).hexdigest()
        cached = self._cache.get(key)
        if cached is None:
            with self._lock:
                try:
                    raw = self._session.evaluate(
                        self._wl.Global.softarmExecuteJSON(serialized)
                    )
                    result = json.loads(raw)
                except Exception as error:
                    raise SymbolicBackendError(
                        f"Wolfram {operation} transport failed: {error}"
                    ) from error
            if not result.get("ok"):
                raise SymbolicBackendError(
                    f"Wolfram {operation} failed: {result.get('error', 'unknown error')}"
                )
            self._cache[key] = result
        else:
            result = cached
        if operation == "cse" or "expressions" not in result:
            return result
        decoded = decode_dag(result["expressions"], symbol_table)
        allowed = set(symbol_table)
        for expression in decoded:
            unexpected = {str(item) for item in expression.free_symbols} - allowed
            if unexpected:
                raise SymbolicBackendError(
                    f"Wolfram {operation} introduced symbols: {sorted(unexpected)}"
                )
            if expression.has(sp.Derivative, sp.Integral):
                raise SymbolicBackendError(
                    f"Wolfram {operation} returned an unevaluated calculus operator"
                )
        return decoded

    def optimize(self, expressions: Sequence[sp.Expr]) -> list[sp.Expr]:
        # FactorTerms can be dramatically slower than calculus on already
        # expanded multi-section dynamics. DAG size is cheap to compute and
        # avoids sending pathological normalization jobs to the kernel.
        if len(encode_dag(expressions)["nodes"]) > 4000:
            return list(expressions)
        return list(self._request("optimize", expressions))

    def differentiate(
        self, expressions: Sequence[sp.Expr], variables: Sequence[sp.Symbol]
    ) -> list[sp.Expr]:
        return list(self._request(
            "differentiate", expressions, variables=list(variables)
        ))

    def jacobian(
        self, matrix: sp.MatrixBase, variables: Sequence[sp.Symbol]
    ) -> sp.Matrix:
        if not variables:
            return sp.zeros(matrix.rows * matrix.cols, 0)
        values = self.differentiate(_flatten(matrix), variables)
        return sp.Matrix(matrix.rows * matrix.cols, len(variables), values)

    def integrate(
        self,
        expression: sp.Expr | sp.MatrixBase,
        variable: sp.Symbol,
        lower: sp.Expr,
        upper: sp.Expr,
    ) -> sp.Expr | sp.Matrix:
        matrix = expression if isinstance(expression, sp.MatrixBase) else sp.Matrix([expression])
        values = list(self._request(
            "integrate",
            _flatten(matrix),
            variable=variable,
            lower=lower,
            upper=upper,
        ))
        result = sp.Matrix(matrix.rows, matrix.cols, values)
        return result if isinstance(expression, sp.MatrixBase) else result[0]

    def substitute(
        self,
        expression: sp.Expr | sp.MatrixBase,
        substitutions: dict[sp.Symbol, sp.Expr],
    ) -> sp.Expr | sp.Matrix:
        matrix = expression if isinstance(expression, sp.MatrixBase) else sp.Matrix([expression])
        variables = list(substitutions)
        values = list(substitutions.values())
        result = list(self._request(
            "substitute",
            _flatten(matrix),
            variables=variables,
            values=values,
        ))
        decoded = sp.Matrix(matrix.rows, matrix.cols, result)
        return decoded if isinstance(expression, sp.MatrixBase) else decoded[0]

    def cse(
        self,
        expressions: Sequence[sp.Expr],
        prefix: str = "t",
        order: str = "none",
    ) -> tuple[list[tuple[sp.Symbol, sp.Expr]], list[sp.Expr]]:
        # Experimental`OptimizeExpression is intentionally best-effort. The
        # stable, correctness-preserving fallback remains SymPy CSE.
        try:
            result = self._request("cse", expressions, prefix=prefix)
            if not isinstance(result, dict) or "replacements" not in result:
                raise SymbolicBackendError("unsupported optimized-expression shape")
            symbols = _symbols(expressions)
            replacements_raw = decode_dag(result["replacements"], symbols)
            reduced = decode_dag(result["expressions"], symbols)
            names = result["replacement_names"]
            replacements = [
                (sp.Symbol(name), expression)
                for name, expression in zip(names, replacements_raw, strict=True)
            ]
            original_symbols = set().union(
                *(expression.free_symbols for expression in expressions)
            )
            available = set(original_symbols)
            allowed_functions = {
                node.func.__name__
                for expression in expressions
                for node in sp.preorder_traversal(expression)
                if isinstance(node, sp.Function)
            }
            for symbol, expression in replacements:
                if not expression.free_symbols <= available:
                    raise SymbolicBackendError("CSE replacement dependency is invalid")
                introduced_functions = {
                    node.func.__name__
                    for node in sp.preorder_traversal(expression)
                    if isinstance(node, sp.Function)
                } - allowed_functions
                if introduced_functions:
                    raise SymbolicBackendError(
                        f"CSE introduced functions: {sorted(introduced_functions)}"
                    )
                available.add(symbol)
            if any(not expression.free_symbols <= available for expression in reduced):
                raise SymbolicBackendError("CSE reduced expression dependency is invalid")
            return replacements, reduced
        except (SymbolicBackendError, KeyError, TypeError, ValueError):
            return super().cse(expressions, prefix, order)

    def rank(self, matrix: sp.MatrixBase) -> int:
        result = self._request(
            "rank", _flatten(matrix), rows=matrix.rows, cols=matrix.cols
        )
        if not isinstance(result, dict) or "rank" not in result:
            raise SymbolicBackendError("Wolfram rank returned no result")
        return int(result["rank"])

    def close(self) -> None:
        session = getattr(self, "_session", None)
        if session is None:
            return
        try:
            session.terminate()
        except Exception:
            pass
        self._session = None

    def __del__(self):
        self.close()


def create_session(
    backend: str = "sympy",
    wolfram_kernel: str | None = None,
) -> SymbolicSession:
    if backend == "sympy":
        return SymbolicSession()
    if backend == "wolfram":
        return WolframSession(wolfram_kernel)
    raise ValueError(f"unknown CAS backend: {backend}")
