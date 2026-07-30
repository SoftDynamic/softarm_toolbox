from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys

from .codegen import generate_matlab_bundle
from .actuation import derive_actuation
from .config import ConfigError, load_config
from .derive import derive


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog="softarm")
    commands = parser.add_subparsers(dest="command", required=True)
    build = commands.add_parser("build", help="derive a model and generate numerical code")
    build.add_argument("config")
    build.add_argument("--backend", choices=("sympy", "wolfram"), default="sympy")
    build.add_argument("--target", choices=("matlab",), default="matlab")
    build.add_argument("--out", required=True)
    build.add_argument("--wolfram-kernel")
    validate = commands.add_parser("validate", help="validate a model configuration")
    validate.add_argument("config")
    inspect = commands.add_parser("inspect", help="print a generated bundle's minimal manifest")
    inspect.add_argument("bundle")
    return parser


def main(argv: list[str] | None = None) -> int:
    args = _parser().parse_args(argv)
    try:
        if args.command == "validate":
            config = load_config(args.config)
            plant = derive(config)
            actuation = derive_actuation(plant)
            suffix = "" if actuation is None else f" and {actuation.count} actuator channel(s)"
            print(f"valid {config.family} configuration with {config.segments} segment(s){suffix}")
            return 0
        if args.command == "inspect":
            manifest = Path(args.bundle, "manifest.json")
            print(json.dumps(json.loads(manifest.read_text(encoding="utf-8")), indent=2))
            return 0
        config = load_config(args.config)
        plant = derive(config)
        actuation = derive_actuation(plant)
        output = generate_matlab_bundle(
            plant, args.out, args.backend, args.wolfram_kernel, actuation
        )
        print(output)
        return 0
    except (ConfigError, RuntimeError, OSError, ValueError) as error:
        print(f"softarm: {error}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
