#!/usr/bin/env python3
"""Generate candidate formulas for every composition space in config.yaml.

Reads composition-space definitions from config.yaml and generates
charge-balanced candidate formulas. No chemistry is hardcoded — users
define their own spaces via the ``composition_spaces`` section.

For a single system with full control over novelty filtering and diversity
options, prefer the dedicated script in that system's example folder.

Usage:
    python examples/generate_candidates.py                  # 20 per space
    python examples/generate_candidates.py -n 100           # 100 per space
    python examples/generate_candidates.py --filter-known   # exclude known
    python examples/generate_candidates.py -n 100 -o out.txt  # out.<space>.txt

Requires the package to be installed:
    pip install -e .
"""

import argparse
from pathlib import Path

from candidate_cli import build_arg_parser, run
from config import get as cfg_get
from composition import CompositionSpace


def load_space(name: str) -> CompositionSpace:
    """Load a composition space by name from config.yaml."""
    spaces = cfg_get("composition_spaces", {}) or {}
    if name not in spaces:
        raise KeyError(f"Composition space {name!r} not found in config.yaml. "
                       f"Available: {sorted(spaces)}")
    return CompositionSpace(name, spaces[name])


def _per_space_output(output: str, space_name: str) -> str:
    """Insert the space name into an output path: out.txt -> out.<space>.txt."""
    p = Path(output)
    return str(p.with_name(f"{p.stem}.{space_name}{p.suffix or '.txt'}"))


def main() -> None:
    parser = build_arg_parser(
        "Generate candidate formulas for all composition spaces in config.yaml."
    )
    args = parser.parse_args()

    spaces_cfg = cfg_get("composition_spaces", {}) or {}
    if not spaces_cfg:
        print("No composition_spaces defined in config.yaml")
        return

    for space_name in sorted(spaces_cfg):
        space = load_space(space_name)
        space_args = argparse.Namespace(**vars(args))
        if args.output:
            space_args.output = _per_space_output(args.output, space_name)
        run(space, space_args, title=f"Composition Space: {space_name}")
        print()


if __name__ == "__main__":
    main()
