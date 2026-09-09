#!/usr/bin/env bash
# Simple end-to-end novel-formula search:
#   1) generate SMACT candidates
#   2) screen them for novelty against a known-formulas list
#   3) optionally check the novel results against the live databases
#
# Usage:
#   ./run_formula_search.sh                 # use config.yaml defaults
#   ./run_formula_search.sh "Al Ti Nb W" 3  # elements, num-elements
#   CHECK_LIVE=1 ./run_formula_search.sh    # also run the live DB lookup

set -euo pipefail

# --- resolve project root (so this works from any directory) --------------
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

PYTHON="${PYTHON:-python}"

# --- inputs / outputs ------------------------------------------------------
ELEMENTS="${1:-}"            # optional: space-separated element universe
NUM_ELEMENTS="${2:-}"        # optional: combination sizes, e.g. 3 or "2 3 4"
MAX_STOICH="${MAX_STOICH:-8}"

KNOWN="${KNOWN:-data/formulas/known_formulas.csv}"     # your reference list
CAND="${CAND:-/tmp/candidates.csv}"                    # step 1 output
NOVEL="${NOVEL:-data/formulas/novel_formulas.csv}"     # step 2 output
CHECK_LIVE="${CHECK_LIVE:-0}"                          # set to 1 for step 3

# --- 1) generate candidates ------------------------------------------------
gen_args=()
[ -n "$ELEMENTS" ]    && gen_args+=(--elements $ELEMENTS)
[ -n "$NUM_ELEMENTS" ] && gen_args+=(--num-elements $NUM_ELEMENTS)
gen_args+=(--max-stoich "$MAX_STOICH" --out "$CAND")

echo ">> [1/2] Generating SMACT candidates -> $CAND"
"$PYTHON" scripts/generate_formulas.py "${gen_args[@]}"

# --- 2) screen for novelty -------------------------------------------------
echo ">> [2/2] Screening candidates against $KNOWN -> $NOVEL"
mkdir -p "$(dirname "$NOVEL")"
"$PYTHON" scripts/find_novel_formulas.py \
    --candidates "$CAND" \
    --known "$KNOWN" \
    --out "$NOVEL" \
    --novel-only

echo
echo "Done. Novel candidates written to: $NOVEL"
echo "Total novel rows: $(($(wc -l < "$NOVEL") - 1))"

# --- 3) optional live database lookup -------------------------------------
if [ "$CHECK_LIVE" = "1" ]; then
    echo
    echo ">> [3] Checking novel formulas against live databases (may be slow/flaky)"
    "$PYTHON" scripts/search_known_formulas.py "$NOVEL" --output "${NOVEL%.csv}_live.csv"
fi
