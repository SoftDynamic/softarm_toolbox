from __future__ import annotations

from pathlib import Path

import pytest
import sympy as sp

from softarm import cli as cli_module
from softarm import pipeline as pipeline_module
from softarm.backends.wolfram import WolframUnavailableError
from softarm.cli import main as cli_main
from softarm.codegen.optimization import SympyCse
from softarm.config import DynamicsConfig, IntegrationConfig, ModelConfig, load_config
from softarm.dynamics import SympyBatchDifferentiator
from softarm.pipeline import BuildError, BuildOptions, build_bundle


def _config() -> ModelConfig:
    return ModelConfig(
        rod="euler_bernoulli",
        parameterization="ritz",
        segments=1,
        dynamics=DynamicsConfig("symbolic_lagrange"),
        integration=IntegrationConfig("analytic"),
        ritz_x=(0.0, 0.0, 1.5, -0.5),
        ritz_y=(0.0, 0.0, 1.5, -0.5),
    )


class StepClock:
    def __init__(self):
        self.value = -1.0

    def __call__(self):
        self.value += 1.0
        return self.value


def _timing_line(output: str, stage: str) -> str:
    return next(line for line in output.splitlines() if stage in line)


def _assert_sections_in_order(output: str, sections: list[str]) -> None:
    positions = [output.index(f"\n{section}\n") for section in sections]
    assert positions == sorted(positions)


class RecordingKernel:
    events: list[str] = []
    fail_cse = False

    def __init__(self, kernel=None, timeout=600.0):
        self.events.append(f"start:{kernel}:{timeout}")

    def __enter__(self):
        self.events.append("enter")
        return self

    def __exit__(self, exc_type, exc, traceback):
        self.events.append("close")

    def differentiate(self, expressions, variables):
        self.events.append("differentiate")
        return SympyBatchDifferentiator().differentiate(expressions, variables)

    def factor_terms(self, expressions):
        self.events.append("factor_terms")
        return list(expressions)

    def cse(self, expressions, prefix="t", order="none"):
        self.events.append("wolfram_cse")
        if self.fail_cse:
            raise RuntimeError("experimental CSE failed")
        return SympyCse().cse(expressions, prefix, order)


def _fake_generate(
    plant,
    output,
    actuation=None,
    constraint=None,
    tex_appendix=False,
    optimizer=None,
):
    del plant, actuation, constraint, tex_appendix
    x = sp.Symbol("x")
    optimizer.optimize([x * (x + 1)], (1,))
    return Path(output)


def test_explicit_wolfram_steps_share_one_kernel(monkeypatch, tmp_path):
    RecordingKernel.events = []
    RecordingKernel.fail_cse = False
    monkeypatch.setattr(pipeline_module, "WolframKernel", RecordingKernel)
    monkeypatch.setattr(
        pipeline_module, "generate_matlab_bundle", _fake_generate
    )
    result = build_bundle(
        _config(),
        tmp_path,
        BuildOptions(
            backend="wolfram",
            wolfram_timeout=12.0,
            wolfram_factor_terms=True,
            wolfram_cse=True,
        ),
    )
    assert result == tmp_path
    assert RecordingKernel.events == [
        "start:None:12.0",
        "enter",
        "differentiate",
        "factor_terms",
        "wolfram_cse",
        "close",
    ]


def test_selected_wolfram_cse_failure_does_not_fallback(monkeypatch, tmp_path):
    RecordingKernel.events = []
    RecordingKernel.fail_cse = True
    monkeypatch.setattr(pipeline_module, "WolframKernel", RecordingKernel)
    monkeypatch.setattr(
        pipeline_module, "generate_matlab_bundle", _fake_generate
    )
    with pytest.raises(BuildError, match="No alternative strategy was attempted"):
        build_bundle(
            _config(),
            tmp_path,
            BuildOptions(backend="wolfram", wolfram_cse=True),
        )
    assert RecordingKernel.events.count("wolfram_cse") == 1
    assert RecordingKernel.events[-1] == "close"


def test_cli_rejects_wolfram_flags_with_sympy(capsys):
    assert cli_main([
        "build",
        "unused.toml",
        "--out",
        "unused",
        "--wolfram-cse",
    ]) == 2
    assert "require --backend wolfram" in capsys.readouterr().err


def test_cli_passes_verbose_build_option(monkeypatch, tmp_path):
    captured = None

    def fake_build(config, output, options):
        nonlocal captured
        del config
        captured = options
        return Path(output)

    monkeypatch.setattr(cli_module, "load_config", lambda path: _config())
    monkeypatch.setattr(cli_module, "build_bundle", fake_build)
    assert cli_main([
        "build",
        "unused.toml",
        "--out",
        str(tmp_path),
        "--verbose",
    ]) == 0
    assert captured is not None
    assert captured.verbose is True


def test_cli_success_is_silent_without_verbose(monkeypatch, tmp_path, capsys):
    monkeypatch.setattr(cli_module, "load_config", lambda path: _config())
    monkeypatch.setattr(
        cli_module,
        "build_bundle",
        lambda config, output, options: Path(output),
    )
    assert cli_main([
        "build",
        "unused.toml",
        "--out",
        str(tmp_path),
    ]) == 0
    captured = capsys.readouterr()
    assert captured.out == ""
    assert captured.err == ""


def test_default_build_does_not_print_timings(monkeypatch, tmp_path, capsys):
    monkeypatch.setattr(
        pipeline_module, "generate_matlab_bundle", _fake_generate
    )
    build_bundle(_config(), tmp_path)
    assert "Timing:" not in capsys.readouterr().out


def test_verbose_sympy_build_prints_stage_timings(
    monkeypatch, tmp_path, capsys
):
    monkeypatch.setattr(pipeline_module.time, "perf_counter", StepClock())
    monkeypatch.setattr(
        pipeline_module, "generate_matlab_bundle", _fake_generate
    )
    build_bundle(_config(), tmp_path, BuildOptions(verbose=True))
    output = capsys.readouterr().out
    _assert_sections_in_order(output, [
        "Configuration",
        "Symbolic Strategy",
        "Build Stages",
        "Resolved Model",
        "Generated Artifacts",
        "Build Result",
    ])
    assert "SymPy " + sp.__version__ in output
    assert "1.000 s" in _timing_line(output, "Model derivation [SymPy]")
    assert "1.000 s" in _timing_line(output, "Bias differentiation [SymPy]")
    assert "1.000 s" in _timing_line(output, "MATLAB generation/CSE")
    assert "7.000 s" in _timing_line(output, "Total build")
    assert "1.000 s" in _timing_line(output, "Diagnostics analysis")
    assert "SUCCESS" in output


def test_verbose_wolfram_build_prints_kernel_timing(
    monkeypatch, tmp_path, capsys
):
    RecordingKernel.events = []
    RecordingKernel.fail_cse = False
    monkeypatch.setattr(pipeline_module.time, "perf_counter", StepClock())
    monkeypatch.setattr(pipeline_module, "WolframKernel", RecordingKernel)
    monkeypatch.setattr(
        pipeline_module, "generate_matlab_bundle", _fake_generate
    )
    build_bundle(
        _config(), tmp_path,
        BuildOptions(backend="wolfram", verbose=True),
    )
    output = capsys.readouterr().out
    assert "Wolfram Kernel" in output
    assert "1.000 s" in _timing_line(output, "Model derivation [SymPy]")
    assert "1.000 s" in _timing_line(output, "Wolfram Kernel startup")
    assert "1.000 s" in _timing_line(output, "Bias differentiation [Wolfram]")
    assert "1.000 s" in _timing_line(output, "MATLAB generation/CSE")
    assert "9.000 s" in _timing_line(output, "Total build")
    assert "1.000 s" in _timing_line(output, "Diagnostics analysis")


def test_verbose_model_derivation_failure_prints_timings(
    monkeypatch, tmp_path, capsys
):
    monkeypatch.setattr(pipeline_module.time, "perf_counter", StepClock())
    monkeypatch.setattr(
        pipeline_module,
        "derive_system",
        lambda config: (_ for _ in ()).throw(ValueError("bad model")),
    )
    with pytest.raises(BuildError, match="bad model"):
        build_bundle(_config(), tmp_path, BuildOptions(verbose=True))
    output = capsys.readouterr().out
    model_line = _timing_line(output, "Model derivation [SymPy]")
    total_line = _timing_line(output, "Total build")
    assert "FAIL" in model_line and "1.000 s (failed)" in model_line
    assert "FAIL" in total_line and "3.000 s (failed)" in total_line
    assert "Resolved Model" not in output
    assert output.rstrip().endswith("FAILED")


def test_verbose_kernel_startup_failure_prints_timings(
    monkeypatch, tmp_path, capsys
):
    def fail_kernel(*args, **kwargs):
        raise WolframUnavailableError("no kernel")

    monkeypatch.setattr(pipeline_module.time, "perf_counter", StepClock())
    monkeypatch.setattr(pipeline_module, "WolframKernel", fail_kernel)
    with pytest.raises(BuildError, match="no kernel"):
        build_bundle(
            _config(), tmp_path,
            BuildOptions(backend="wolfram", verbose=True),
        )
    output = capsys.readouterr().out
    kernel_line = _timing_line(output, "Wolfram Kernel startup")
    total_line = _timing_line(output, "Total build")
    assert "FAIL" in kernel_line and "1.000 s (failed)" in kernel_line
    assert "FAIL" in total_line and "5.000 s (failed)" in total_line
    assert "Resolved Model" in output


def test_verbose_bias_failure_prints_timings(monkeypatch, tmp_path, capsys):
    monkeypatch.setattr(pipeline_module.time, "perf_counter", StepClock())
    monkeypatch.setattr(
        pipeline_module,
        "assemble_bias",
        lambda plant, differentiator: (_ for _ in ()).throw(
            ValueError("bad bias")
        ),
    )
    with pytest.raises(BuildError, match="bad bias"):
        build_bundle(_config(), tmp_path, BuildOptions(verbose=True))
    output = capsys.readouterr().out
    bias_line = _timing_line(output, "Bias differentiation [SymPy]")
    total_line = _timing_line(output, "Total build")
    assert "FAIL" in bias_line and "1.000 s (failed)" in bias_line
    assert "FAIL" in total_line and "5.000 s (failed)" in total_line


def test_verbose_codegen_failure_prints_timings(monkeypatch, tmp_path, capsys):
    def fail_generate(*args, **kwargs):
        raise ValueError("bad codegen")

    monkeypatch.setattr(pipeline_module.time, "perf_counter", StepClock())
    monkeypatch.setattr(
        pipeline_module, "generate_matlab_bundle", fail_generate
    )
    with pytest.raises(BuildError, match="bad codegen"):
        build_bundle(_config(), tmp_path, BuildOptions(verbose=True))
    output = capsys.readouterr().out
    codegen_line = _timing_line(output, "MATLAB generation/CSE")
    total_line = _timing_line(output, "Total build")
    assert "FAIL" in codegen_line and "1.000 s (failed)" in codegen_line
    assert "FAIL" in total_line and "7.000 s (failed)" in total_line


def test_verbose_report_resolves_defaults_and_function_metrics(
    tmp_path, capsys
):
    config_path = tmp_path / "model.toml"
    config_path.write_text(
        "# raw comments must not be echoed\n"
        '[model]\nrod="euler_bernoulli"\nparameterization="ritz"\nsegments=1\n'
        '[dynamics]\nformulation="symbolic_lagrange"\n'
        '[ritz]\nx=[0,0,1.5,-0.5]\ny=[0,0,1.5,-0.5]\n',
        encoding="utf-8",
    )
    output_path = tmp_path / "bundle"
    build_bundle(
        load_config(config_path),
        output_path,
        BuildOptions(verbose=True),
    )
    output = capsys.readouterr().out
    _assert_sections_in_order(output, [
        "Configuration",
        "Symbolic Strategy",
        "Build Stages",
        "Resolved Model",
        "Exported Symbolic Functions",
        "Generated Artifacts",
        "Build Result",
    ])
    assert str(config_path.resolve()) in output
    assert str(output_path.resolve()) in output
    assert "raw comments must not be echoed" not in output
    assert "distributed (default)" in output
    assert "tip_mass" in output
    assert "softarm_mass" in output
    assert "softarm_bias" in output
    assert "Raw ops" in output
    assert "CSE temps" in output
    assert "Generated ops" in output
    assert "MATLAB files" in output
    assert output.rstrip().endswith(str(output_path.resolve()))


def test_verbose_wolfram_metadata_uses_kernel_values(
    monkeypatch, tmp_path, capsys
):
    executable = tmp_path / "actual-kernel.exe"

    class MetadataKernel(RecordingKernel):
        def __init__(self, kernel=None, timeout=600.0):
            super().__init__(kernel, timeout)
            self.executable = executable
            self.version = "Actual Test Kernel 1.2.3"

    monkeypatch.setattr(pipeline_module, "WolframKernel", MetadataKernel)
    monkeypatch.setattr(
        pipeline_module, "generate_matlab_bundle", _fake_generate
    )
    build_bundle(
        _config(), tmp_path,
        BuildOptions(backend="wolfram", verbose=True),
    )
    output = capsys.readouterr().out
    assert str(executable.resolve()) in output
    assert "Actual Test Kernel 1.2.3" in output


def test_verbose_omits_unavailable_wolfram_metadata(
    monkeypatch, tmp_path, capsys
):
    monkeypatch.setattr(pipeline_module, "WolframKernel", RecordingKernel)
    monkeypatch.setattr(
        pipeline_module, "generate_matlab_bundle", _fake_generate
    )
    build_bundle(
        _config(), tmp_path,
        BuildOptions(backend="wolfram", verbose=True),
    )
    output = capsys.readouterr().out
    kernel_section = output.split("Wolfram Kernel", 1)[1].split("Build Stages", 1)[0]
    assert "Executable" not in kernel_section
    assert "Version" not in kernel_section


def test_verbose_strategy_reports_selected_wolfram_switches(
    monkeypatch, tmp_path, capsys
):
    RecordingKernel.events = []
    RecordingKernel.fail_cse = False
    monkeypatch.setattr(pipeline_module, "WolframKernel", RecordingKernel)
    monkeypatch.setattr(
        pipeline_module, "generate_matlab_bundle", _fake_generate
    )
    build_bundle(
        _config(),
        tmp_path,
        BuildOptions(
            backend="wolfram",
            wolfram_factor_terms=True,
            wolfram_cse=True,
            tex_appendix=True,
            verbose=True,
        ),
    )
    output = capsys.readouterr().out
    assert "Wolfram FactorTerms (enabled)" in output
    assert "Wolfram experimental CSE" in output
    assert "TeX appendix" in output and "enabled" in output


def test_verbose_diagnostics_do_not_change_generated_bundle(
    tmp_path, capsys
):
    quiet_output = tmp_path / "quiet"
    verbose_output = tmp_path / "verbose"
    build_bundle(_config(), quiet_output)
    build_bundle(
        _config(), verbose_output, BuildOptions(verbose=True)
    )
    capsys.readouterr()
    quiet_files = {
        path.name: path.read_bytes()
        for path in quiet_output.iterdir() if path.is_file()
    }
    verbose_files = {
        path.name: path.read_bytes()
        for path in verbose_output.iterdir() if path.is_file()
    }
    assert verbose_files == quiet_files


def test_factor_terms_error_guides_without_retry():
    error = BuildError(
        "code generation",
        "Wolfram FactorTerms",
        RuntimeError("Wolfram factor_terms failed: timed out"),
    )
    assert "No alternative strategy was attempted" in str(error)
    assert "--wolfram-factor-terms" in str(error)
