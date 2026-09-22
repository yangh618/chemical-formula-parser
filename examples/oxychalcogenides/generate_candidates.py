#!/usr/bin/env python3
"""Generate candidate formulas for the oxychalcogenide composition space.

Oxychalcogenides contain both oxygen and a chalcogen (S, Se, Te) as anions,
combined with rare-earth and/or transition-metal cations. These materials are
of interest for photocatalysis, thermoelectrics, and energy-storage applications.

This script reads the ``oxychalcogenides`` composition space from config.yaml,
enumerates charge-balanced candidate formulas, screens each one with SMACT
(charge neutrality under an oxidation-state set plus the Pauling
electronegativity test), optionally discards those already reported in
COD/MP/OQMD, and returns a maximally spread-out subset so the output covers many
distinct chemical species rather than clustering on a few.

Usage:
    python generate_candidates.py -n 1000 --filter-known -o novel_1000.txt
    python generate_candidates.py -n 100 --filter-known --db cod mp
    python generate_candidates.py -n 100 --no-diverse      # enumeration order
    python generate_candidates.py -n 100 --no-smact        # skip SMACT screen

Requires the package to be installed:
    pip install -e .
"""

from candidate_cli import build_arg_parser, run
from config import get as cfg_get
from composition import CompositionSpace

SPACE_NAME = "oxychalcogenides"
TITLE = "Oxychalcogenide Composition Space"

EPILOG = """\
examples:
  # 1000 novel, maximally diverse oxychalcogenides saved to a file
  python generate_candidates.py -n 1000 --filter-known -o novel_1000.txt

  # filter against COD and MP only
  python generate_candidates.py -n 500 --filter-known --db cod mp

  # enumerate in generation order (clustered, not diverse)
  python generate_candidates.py -n 500 --filter-known --no-diverse

  # emphasise distinct chemical species over distinct stoichiometry
  python generate_candidates.py -n 500 --w-element-set 3
"""


def load_space() -> CompositionSpace:
    """Load the oxychalcogenide composition space from config.yaml."""
    spaces = cfg_get("composition_spaces", {}) or {}
    if SPACE_NAME not in spaces:
        raise KeyError(
            f"Composition space {SPACE_NAME!r} not found in config.yaml. "
            f"Available: {sorted(spaces)}"
        )
    return CompositionSpace(SPACE_NAME, spaces[SPACE_NAME])


def main() -> None:
    parser = build_arg_parser(
        "Generate oxychalcogenide candidate formulas.", epilog=EPILOG
    )
    args = parser.parse_args()
    run(load_space(), args, TITLE)


if __name__ == "__main__":
    main()
