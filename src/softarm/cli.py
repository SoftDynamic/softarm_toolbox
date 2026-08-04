from __future__ import annotations

import argparse
import json
import subprocess
import sys
import tomllib
from pathlib import Path

from .actuation import derive_actuation
from .config import ConfigError, load_config
from .constraints import derive_constraint
from .derive import derive
from .pipeline import BuildOptions, build_bundle, symbolic_plan


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog="softarm")
    commands = parser.add_subparsers(dest="command", required=True)
    build = commands.add_parser("build", help="derive a model and generate numerical code")
    build.add_argument("config")
    build.add_argument("--backend", choices=("sympy", "wolfram"), default="sympy")
    build.add_argument("--target", choices=("matlab",), default="matlab")
    build.add_argument("--out", required=True)
    build.add_argument("--wolfram-kernel")
    build.add_argument(
        "--wolfram-timeout",
        type=float,
        help="seconds allowed for each explicitly selected Wolfram operation",
    )
    build.add_argument(
        "--wolfram-cse",
        action="store_true",
        help="use experimental Wolfram CSE; failures stop the build",
    )
    build.add_argument(
        "--wolfram-factor-terms",
        action="store_true",
        help="run Wolfram FactorTerms before CSE; may increase build time",
    )
    build.add_argument(
        "--tex-appendix",
        action="store_true",
        help="append CAS-optimized exact symbolic expressions to softarm_model.tex",
    )
    validate = commands.add_parser("validate", help="validate a model configuration")
    validate.add_argument("config")
    inspect = commands.add_parser("inspect", help="print a generated bundle manifest")
    inspect.add_argument("bundle")
    commands.add_parser("matlab", help="run the configured MATLAB executable")
    return parser


def _toolbox_root() -> Path:
    return Path(__file__).resolve().parents[2]


def _matlab_path() -> Path:
    local = _toolbox_root() / ".softarm.local.toml"
    if not local.is_file():
        raise FileNotFoundError(
            f"missing {local}; copy the toolbox .softarm.local.toml.example "
            "and set tools.matlab"
        )
    with local.open("rb") as stream:
        tools = tomllib.load(stream).get("tools")
    candidate = tools.get("matlab") if isinstance(tools, dict) else None
    if not isinstance(candidate, str) or not candidate.strip():
        raise ValueError(f"{local} has no tools.matlab entry")
    executable = Path(candidate)
    if not executable.is_file():
        raise FileNotFoundError(f"configured MATLAB executable does not exist: {executable}")
    return executable


def main(argv: list[str] | None = None) -> int:
    command_args = sys.argv[1:] if argv is None else argv
    try:
        matlab_help = command_args[:1] == ["matlab"] and command_args[1:] in (
            ["-h"],
            ["--help"],
        )
        if command_args[:1] == ["matlab"] and not matlab_help:
            return subprocess.run([str(_matlab_path()), *command_args[1:]], check=False).returncode
        args = _parser().parse_args(command_args)
        if args.command == "validate":
            config = load_config(args.config)
            plant = derive(config)
            actuation = derive_actuation(plant)
            constraint = derive_constraint(plant)
            suffix = "" if actuation is None else f" and {actuation.count} actuator channel(s)"
            if constraint is not None:
                suffix += f" and {constraint.count} constraint channel(s)"
            print(
                f"valid {config.rod} + {config.parameterization} configuration "
                f"with {config.segments} segment(s){suffix}"
            )
            return 0
        if args.command == "inspect":
            manifest = Path(args.bundle, "manifest.json")
            print(json.dumps(json.loads(manifest.read_text(encoding="utf-8")), indent=2))
            return 0
        wolfram_option_selected = (
            args.wolfram_kernel is not None
            or args.wolfram_timeout is not None
            or args.wolfram_cse
            or args.wolfram_factor_terms
        )
        if args.backend != "wolfram" and wolfram_option_selected:
            raise ValueError(
                "--wolfram-kernel, --wolfram-timeout, --wolfram-cse, and "
                "--wolfram-factor-terms require --backend wolfram"
            )
        options = BuildOptions(
            backend=args.backend,
            wolfram_kernel=args.wolfram_kernel,
            wolfram_timeout=(
                600.0 if args.wolfram_timeout is None else args.wolfram_timeout
            ),
            wolfram_cse=args.wolfram_cse,
            wolfram_factor_terms=args.wolfram_factor_terms,
            tex_appendix=args.tex_appendix,
        )
        config = load_config(args.config)
        print(symbolic_plan(options))
        output = build_bundle(config, args.out, options)
        print(output)
        return 0
    except (ConfigError, RuntimeError, OSError, ValueError) as error:
        print(f"softarm: {error}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
