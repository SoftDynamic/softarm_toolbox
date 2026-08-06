import pytest
import sympy as sp

from softarm.backends.protocol import decode_dag, encode_dag
from softarm.backends.wolfram import WolframKernel, WolframProtocolError
from softarm.codegen.optimization import FunctionOptimizer, SympyCse
from softarm.config import DynamicsConfig, IntegrationConfig, ModelConfig
from softarm.derive import derive
from softarm.dynamics import SympyBatchDifferentiator, assemble_bias
from softarm.pipeline import BuildOptions, symbolic_plan


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


def test_wolfram_capabilities_reject_output_dimension_changes(monkeypatch):
    kernel = object.__new__(WolframKernel)
    monkeypatch.setattr(kernel, "_request", lambda *args, **kwargs: [])
    x = sp.Symbol("x")
    with pytest.raises(WolframProtocolError, match="FactorTerms returned"):
        kernel.factor_terms([x])
    with pytest.raises(WolframProtocolError, match="expected 1"):
        kernel.differentiate([x], [x])


def test_sympy_differentiator_and_function_optimizer():
    x = sp.Symbol("x", real=True)
    differentiator = SympyBatchDifferentiator()
    assert differentiator.differentiate([x**3], [x]) == [3 * x**2]
    optimized = FunctionOptimizer().optimize(
        [x**2 + sp.sin(x**2), x**2 + sp.cos(x**2)], (2,)
    )
    assert optimized.replacements
    rebuilt = list(optimized.expressions)
    for symbol, expression in reversed(optimized.replacements):
        rebuilt = [item.xreplace({symbol: expression}) for item in rebuilt]
    assert rebuilt == [x**2 + sp.sin(x**2), x**2 + sp.cos(x**2)]


def test_bias_uses_the_single_batch_derivative_formula():
    plant = derive(ModelConfig(
        rod="euler_bernoulli",
        parameterization="ritz",
        segments=1,
        dynamics=DynamicsConfig("symbolic_lagrange"),
        integration=IntegrationConfig("analytic"),
        ritz_x=(0.0, 0.0, 1.5, -0.5),
        ritz_y=(0.0, 0.0, 1.5, -0.5),
    ))
    assert assemble_bias(plant, SympyBatchDifferentiator()) == plant.bias


def test_explicit_optimizer_order_and_no_fallback():
    events: list[str] = []

    class RecordingNormalizer:
        def factor_terms(self, expressions):
            events.append("factor_terms")
            return list(expressions)

    class FailingCse:
        def cse(self, expressions, prefix="t", order="none"):
            events.append("wolfram_cse")
            raise RuntimeError("selected CSE failed")

    optimizer = FunctionOptimizer(
        normalizer=RecordingNormalizer(), eliminator=FailingCse()
    )
    with pytest.raises(RuntimeError, match="selected CSE failed"):
        optimizer.optimize([sp.Symbol("x") + 1], (1,))
    assert events == ["factor_terms", "wolfram_cse"]


def test_default_optimizer_does_not_normalize():
    x = sp.Symbol("x")
    optimized = FunctionOptimizer(eliminator=SympyCse()).optimize([x + 1], (1,))
    assert optimized.expressions == (x + 1,)


def test_function_optimizer_reuses_build_lifetime_result():
    calls = 0

    class CountingCse:
        def cse(self, expressions, prefix="t", order="none"):
            nonlocal calls
            calls += 1
            return [], list(expressions)

    x = sp.Symbol("x")
    optimizer = FunctionOptimizer(eliminator=CountingCse())
    first = optimizer.optimize([x + 1], (1,))
    second = optimizer.optimize([x + 1], (1,))
    assert first is second
    assert calls == 1


def test_build_options_reject_implicit_wolfram_strategy():
    with pytest.raises(ValueError, match="require backend='wolfram'"):
        BuildOptions(wolfram_cse=True)
    with pytest.raises(ValueError, match="require backend='wolfram'"):
        BuildOptions(wolfram_factor_terms=True)
    assert symbolic_plan(BuildOptions(backend="wolfram")) == (
        "Symbolic plan: derive=SymPy, bias=Wolfram, "
        "FactorTerms=off, CSE=SymPy"
    )
