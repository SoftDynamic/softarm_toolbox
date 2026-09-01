from pathlib import Path
from types import SimpleNamespace

from softarm import cli as cli_module
from softarm.cli import _matlab_path
from softarm.cli import main as cli_main


def test_matlab_path_reads_local_tool_config(tmp_path, monkeypatch):
    executable = tmp_path / "matlab.exe"
    executable.touch()
    (tmp_path / ".softarm.local.toml").write_text(
        f'[tools]\nmatlab = "{executable.as_posix()}"\n',
        encoding="utf-8",
    )
    monkeypatch.setattr(
        cli_module,
        "local_tool_config_path",
        lambda: tmp_path / ".softarm.local.toml",
    )

    assert _matlab_path() == executable


def test_matlab_path_requires_local_tool_config(tmp_path, monkeypatch, capsys):
    monkeypatch.setattr(
        cli_module,
        "local_tool_config_path",
        lambda: tmp_path / ".softarm.local.toml",
    )

    assert cli_main(["matlab", "-batch", "disp(1)"]) == 2
    assert "tools.matlab" in capsys.readouterr().err


def test_matlab_path_requires_existing_configured_executable(tmp_path, monkeypatch, capsys):
    (tmp_path / ".softarm.local.toml").write_text(
        '[tools]\nmatlab = "missing-matlab.exe"\n',
        encoding="utf-8",
    )
    monkeypatch.setattr(
        cli_module,
        "local_tool_config_path",
        lambda: tmp_path / ".softarm.local.toml",
    )

    assert cli_main(["matlab", "-batch", "disp(1)"]) == 2
    assert "does not exist" in capsys.readouterr().err


def test_matlab_command_forwards_arguments_and_exit_code(monkeypatch):
    executable = Path("configured-matlab.exe")
    calls = []
    monkeypatch.setattr(cli_module, "_matlab_path", lambda: executable)
    monkeypatch.setattr(
        cli_module.subprocess,
        "run",
        lambda command, check: calls.append((command, check)) or SimpleNamespace(returncode=7),
    )

    assert cli_main(["matlab", "-batch", "disp(1)"]) == 7
    assert calls == [([str(executable), "-batch", "disp(1)"], False)]


def test_casadi_target_passes_affine_terms_to_build_options(
    monkeypatch, tmp_path
):
    captured = None

    def fake_build(config, output, options):
        nonlocal captured
        del config, output
        captured = options
        return tmp_path

    monkeypatch.setattr(cli_module, "load_config", lambda path: object())
    monkeypatch.setattr(cli_module, "build_bundle", fake_build)
    assert cli_main([
        "build",
        "model.toml",
        "--target",
        "casadi",
        "--casadi-affine-terms",
        "64",
        "--out",
        str(tmp_path),
    ]) == 0
    assert captured is not None
    assert captured.target == "casadi"
    assert captured.casadi_affine_terms == 64


def test_casadi_target_rejects_matlab_only_optimization_flags(capsys):
    assert cli_main([
        "build",
        "model.toml",
        "--target",
        "casadi",
        "--backend",
        "wolfram",
        "--wolfram-cse",
        "--out",
        "unused",
    ]) == 2
    assert "CasADi target does not support" in capsys.readouterr().err
