# FAQ

## What does "novel" mean here?

Only *database-unseen relative to your reference list*: the reduced formula does
not appear in `known_formulas.csv` (or, with `search_known_formulas.py`, in the
reachable live databases). It is not a claim of a new material or stable phase.

## Why is my known compound flagged novel?

Because it is not in the reference list you supplied. `Al2O3` is a real compound
but will be marked novel if `known_formulas.csv` does not contain it. Use the live
lookup to cross-check real databases.

## I cannot download the full databases. Do I need them?

No. Generation and novelty screening are fully local. Only the optional live
lookup touches the network, and it queries per-formula (no bulk download).

## Why are some generated formulas strange (e.g. parentheses)?

Canonical reduced formulas come from pymatgen; exotic stoichiometries sometimes
serialise with grouping. Treat candidates as chemically plausible hypotheses, not
verified compounds.

## It is slow.

Reduce the element universe, combination sizes and `max_stoich`:
`python scripts/generate_formulas.py --elements Al Ti Nb --num-elements 2 --max-stoich 6`

## The live lookup times out / errors.

COD/OQMD/JARVIS throttle sustained queries. Retry later or query fewer
formulas; a failing database is reported per-row as `error`, not a crash.

## Where are the 1000 novel formulas?

`data/formulas/novel_1000.csv` (committed). Recreate with the tutorial steps or
`./run_formula_search.sh`.
