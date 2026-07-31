from pathlib import Path

import numpy as np
import pytest
import sympy as sp

from softarm.ast import decode_dag, encode_dag
from softarm.backends.session import SymbolicSession, WolframSession
from softarm.config import IntegrationConfig, ModelConfig
from softarm.derive import derive
from softarm.special import LAMBDA_MODULES, SincSqrt, SincSqrtD


ROOT = Path(__file__).parents[1]


def test_dag_round_trip_preserves_sharing_and_symbol_assumptions():
    x = sp.Symbol("x", real=True, nonnegative=True)
    common = (x + 1) ** 2
    graph = encode_dag([common + 2, common + 3, common])
    assert graph["roots"][2] in graph["nodes"][graph["roots"][0]]["args"]
    decoded = decode_dag(graph, {"x": x})
    assert decoded == [common + 2, common + 3, common]
    assert next(iter(decoded[0].free_symbols)) is x
    wolfram_real = {
        "nodes": [{
            "head": "real",
            "value": "0.1250000000000000000000`18.5",
        }],
        "roots": [0],
    }
    assert float(decode_dag(wolfram_real)[0]) == 0.125


def test_sympy_session_generic_calculus_and_cse():
    x = sp.Symbol("x", real=True)
    session = SymbolicSession()
    assert session.diff(x**3, x) == 3 * x**2
    assert session.integrate(x**2, x, 0, 1) == sp.Rational(1, 3)
    replacements, reduced = session.cse(
        [x**2 + sp.sin(x**2), x**2 + sp.cos(x**2)]
    )
    assert replacements
    rebuilt = list(reduced)
    for symbol, expression in reversed(replacements):
        rebuilt = [item.xreplace({symbol: expression}) for item in rebuilt]
    assert rebuilt == [x**2 + sp.sin(x**2), x**2 + sp.cos(x**2)]


def test_persistent_wolfram_session_calculus_assumptions_and_cse():
    if not (ROOT / ".softarm.local.toml").is_file():
        pytest.skip("no local Wolfram tool configuration")
    x = sp.Symbol("x", real=True, nonnegative=True)
    with WolframSession() as session:
        assert session.diff(SincSqrt(x**2), x) == 2 * x * SincSqrtD(x**2)
        assert session.integrate(x**2, x, 0, 1) == sp.Rational(1, 3)
        before = len(session._cache)
        session.optimize([x * (x + 1)])
        after = len(session._cache)
        session.optimize([x * (x + 1)])
        assert len(session._cache) == after == before + 1
        replacements, reduced = session.cse(
            [x**2 + sp.sin(x**2), x**2 + sp.cos(x**2)]
        )
        assert replacements == [(sp.Symbol("t0"), x**2)]
        assert reduced == [
            sp.Symbol("t0") + sp.sin(sp.Symbol("t0")),
            sp.Symbol("t0") + sp.cos(sp.Symbol("t0")),
        ]


def test_wolfram_full_pcc_results_match_sympy():
    if not (ROOT / ".softarm.local.toml").is_file():
        pytest.skip("no local Wolfram tool configuration")
    config = ModelConfig(
        family="pcc",
        segments=1,
        inertia="lumped",
        integration=IntegrationConfig("analytic"),
    )
    reference = derive(config)
    with WolframSession() as session:
        candidate = derive(config, symbolic=session)
        reference_values = (
            list(reference.q) + list(reference.dq)
            + [item.symbol for item in reference.parameters]
        )
        candidate_values = (
            list(candidate.q) + list(candidate.dq)
            + [item.symbol for item in candidate.parameters]
        )
        numeric = [0.03, -0.02, 0.51, 0.1, -0.05, 0.02]
        numeric += [item.default for item in reference.parameters]
        for expected, actual in (
            (reference.mass, candidate.mass),
            (reference.bias, candidate.bias),
            (reference.end_jacobian, candidate.end_jacobian),
        ):
            expected_function = sp.lambdify(
                reference_values, expected, [LAMBDA_MODULES, "numpy"]
            )
            actual_function = sp.lambdify(
                candidate_values, actual, [LAMBDA_MODULES, "numpy"]
            )
            np.testing.assert_allclose(
                expected_function(*numeric),
                actual_function(*numeric),
                rtol=1e-10,
                atol=1e-12,
            )
