from __future__ import annotations

import hashlib
import json
import os
import shutil
import threading
import tomllib
from collections.abc import Iterable, Sequence
from contextlib import AbstractContextManager
from pathlib import Path

import sympy as sp

from ..special import SPECIAL_DERIVATIVE_HEADS, SPECIAL_FUNCTION_HEADS
from .protocol import decode_dag, encode_dag


class WolframError(RuntimeError):
    pass


class WolframUnavailableError(WolframError):
    pass


class WolframTransportError(WolframError):
    pass


class WolframOperationError(WolframError):
    def __init__(self, operation: str, detail: str):
        self.operation = operation
        super().__init__(f"Wolfram {operation} failed: {detail}")


class WolframProtocolError(WolframError):
    pass


def _kernel_path(explicit: str | None) -> Path:
    candidate = explicit or os.environ.get("SOFTARM_WOLFRAM_KERNEL")
    local = Path.cwd() / ".softarm.local.toml"
    if not candidate and local.is_file():
        with local.open("rb") as stream:
            candidate = tomllib.load(stream).get("tools", {}).get("wolfram_math")
    if not candidate:
        candidate = shutil.which("math") or shutil.which("wolframscript")
    path = Path(candidate) if candidate else Path("__missing_wolfram_kernel__")
    if not path.is_file():
        raise WolframUnavailableError(
            "Wolfram backend requested but no kernel was found; pass "
            "--wolfram-kernel or set SOFTARM_WOLFRAM_KERNEL"
        )
    return path


def _symbols(expressions: Iterable[sp.Expr]) -> dict[str, sp.Symbol]:
    result: dict[str, sp.Symbol] = {}
    for expression in expressions:
        for symbol in expression.free_symbols:
            existing = result.get(str(symbol))
            if existing is not None and existing != symbol:
                raise WolframProtocolError(
                    f"inconsistent assumptions for symbol {symbol}"
                )
            result[str(symbol)] = symbol
    return result


class WolframKernel(AbstractContextManager["WolframKernel"]):
    """One explicit persistent Kernel implementing generic symbolic tasks."""

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
            raise WolframUnavailableError(
                "persistent Wolfram sessions require the 'wolfram' optional "
                "dependency (wolframclient)"
            ) from error
        self._wl = wl
        self._session = WolframLanguageSession(
            str(executable), STARTUP_TIMEOUT=self.timeout
        )
        try:
            self._session.start()
            source = str(bridge).replace("\\", "/").replace('"', '\\"')
            self._session.evaluate(
                wlexpr(f'$SoftArmLibraryMode=True; Get["{source}"]')
            )
            self.version = str(self._session.evaluate(wlexpr("$Version")))
        except Exception as error:
            self.close()
            raise WolframUnavailableError(
                f"Wolfram Kernel failed to start: {error}"
            ) from error

    def __exit__(self, exc_type, exc, traceback) -> None:
        self.close(graceful=exc_type is None)

    def _request(
        self,
        operation: str,
        expressions: Sequence[sp.Expr],
        **parameters,
    ) -> list[sp.Expr] | dict:
        canonical = list(expressions)
        calculus_parameters = {"variables"}
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
                    raise WolframTransportError(
                        f"Wolfram {operation} transport failed: {error}"
                    ) from error
            if not result.get("ok"):
                raise WolframOperationError(
                    operation, result.get("error", "unknown error")
                )
            self._cache[key] = result
        else:
            result = cached
        if operation == "cse" or "expressions" not in result:
            return result
        decoded = decode_dag(result["expressions"], symbol_table)
        allowed = set(symbol_table)
        allowed_functions = {
            node.func.__name__
            for expression in canonical
            for node in sp.preorder_traversal(expression)
            if isinstance(node, sp.Function)
        } | {derivative for _, derivative in SPECIAL_DERIVATIVE_HEADS} | set(
            SPECIAL_FUNCTION_HEADS
        )
        for expression in decoded:
            unexpected = {str(item) for item in expression.free_symbols} - allowed
            if unexpected:
                raise WolframProtocolError(
                    f"Wolfram {operation} introduced symbols: {sorted(unexpected)}"
                )
            introduced_functions = {
                node.func.__name__
                for node in sp.preorder_traversal(expression)
                if isinstance(node, sp.Function)
            } - allowed_functions
            if introduced_functions:
                raise WolframProtocolError(
                    f"Wolfram {operation} introduced functions: "
                    f"{sorted(introduced_functions)}"
                )
            if expression.has(sp.Derivative, sp.Integral):
                raise WolframProtocolError(
                    f"Wolfram {operation} returned an unevaluated calculus operator"
                )
        return decoded

    def factor_terms(self, expressions: Sequence[sp.Expr]) -> list[sp.Expr]:
        """Run the explicitly requested FactorTerms pass without heuristics."""
        result = list(self._request("factor_terms", expressions))
        if len(result) != len(expressions):
            raise WolframProtocolError(
                "FactorTerms returned an unexpected number of expressions"
            )
        return result

    def differentiate(
        self, expressions: Sequence[sp.Expr], variables: Sequence[sp.Symbol]
    ) -> list[sp.Expr]:
        result = list(self._request(
            "differentiate", expressions, variables=list(variables)
        ))
        expected = len(expressions) * len(variables)
        if len(result) != expected:
            raise WolframProtocolError(
                "differentiate returned "
                f"{len(result)} expressions; expected {expected}"
            )
        return result

    def cse(
        self,
        expressions: Sequence[sp.Expr],
        prefix: str = "t",
        order: str = "none",
    ) -> tuple[list[tuple[sp.Symbol, sp.Expr]], list[sp.Expr]]:
        del order
        result = self._request("cse", expressions, prefix=prefix)
        if not isinstance(result, dict) or "replacements" not in result:
            raise WolframProtocolError("unsupported optimized-expression shape")
        symbols = _symbols(expressions)
        try:
            replacements_raw = decode_dag(result["replacements"], symbols)
            reduced = decode_dag(result["expressions"], symbols)
            names = result["replacement_names"]
            replacements = [
                (sp.Symbol(name), expression)
                for name, expression in zip(names, replacements_raw, strict=True)
            ]
        except (KeyError, TypeError, ValueError) as error:
            raise WolframProtocolError(
                f"invalid optimized-expression structure: {error}"
            ) from error
        if len(reduced) != len(expressions):
            raise WolframProtocolError(
                "CSE returned an unexpected number of expressions"
            )
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
                raise WolframProtocolError("CSE replacement dependency is invalid")
            introduced_functions = {
                node.func.__name__
                for node in sp.preorder_traversal(expression)
                if isinstance(node, sp.Function)
            } - allowed_functions
            if introduced_functions:
                raise WolframProtocolError(
                    f"CSE introduced functions: {sorted(introduced_functions)}"
                )
            available.add(symbol)
        if any(not expression.free_symbols <= available for expression in reduced):
            raise WolframProtocolError("CSE reduced expression dependency is invalid")
        reconstructed = list(reduced)
        for symbol, expression in reversed(replacements):
            reconstructed = [item.xreplace({symbol: expression}) for item in reconstructed]
        for original, rebuilt in zip(expressions, reconstructed, strict=True):
            if original != rebuilt and sp.simplify(original - rebuilt) != 0:
                raise WolframProtocolError("CSE reconstruction is not equivalent")
        return replacements, reduced

    def close(self, *, graceful: bool = False) -> None:
        session = getattr(self, "_session", None)
        if session is None:
            return
        try:
            if graceful:
                session.stop()
            else:
                session.terminate()
        except Exception:
            pass
        self._session = None
