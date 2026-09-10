# Manual

Reference for the CrystalFNO novel-composition discovery toolkit. This is a
*download-free* workflow: it generates chemically plausible candidate formulas
with **SMACT**, screens them for novelty, and can query live databases.

## Contents

1. Install
2. Quick start
3. Layout
4. Configuration (`config.yaml`)
5. Scripts
6. Data files
7. Source modules
8. How novelty is decided
9. Troubleshooting

---

## 1. Install

Python 3.11+ .

```bash
conda env create -f environment.yml && conda activate crystalfno-discovery
# or
python -m pip install -r requirements.txt
pip install -e ".[test]"
```

Run the tests with `python -m pytest`.

## 2. Quick start

```bash
# generate candidates, screen them, keep 1000 novel
./run_formula_search.sh            # full config universe
# or run the pieces manually:
python scripts/generate_formulas.py --out /tmp/cand.csv
python scripts/find_novel_formulas.py --candidates /tmp/cand.csv --known data/formulas/known_formulas.csv \
    --out data/formulas/novel_formulas.csv --novel-only
```

A ready-made set is committed at `data/formulas/novel_1000.csv`.

## 3. Layout

```
config.yaml            settings (elements, SMACT, network)
pyproject.toml         packaging
environment.yml        conda env
requirements.txt       pip deps
run_formula_search.sh  one-shot novelty search
manual.md  tutorial.md  README.md  FAQ.md
src/                   importable modules
  config.py formula_utils.py structure_utils.py database_utils.py common.py
  smactlib/  generator.py filters.py      # SMACT wrapper
  novelty/   canonicalize.py filter.py
  validation/ dataset.py
scripts/               generate_formulas.py find_novel_formulas.py search_known_formulas.py
data/formulas/         candidates.csv known_formulas.csv novel_1000.csv
```

## 4. Configuration (`config.yaml`)

| Key | Meaning |
|-----|---------|
| `elements` | element universe for generation |
| `smact.num_elements` | combo sizes, e.g. `[2,3,4]` |
| `smact.max_stoich` | max per-element coefficient |
| `smact.oxidation_states` | oxidation-state set (`icsd24`, ...) |
| `smact.electronegativity` | enable Pauling EN test |
| `smact.include_alloys` | treat all-metal combos as alloys |
| `filters.min/max_elements` | keep formulas within element count |
| `databases.*` | which DBs the live lookup queries |
| `network.*` | HTTP retries/timeouts used by the live lookup |

All scripts read config via `src/config.py`; override the file with the
`CRYSTALFNO_DISCOVERY_CONFIG` env var.

## 5. Scripts

### `scripts/generate_formulas.py`
Enumerate SMACT-valid candidate formulas.

```bash
python scripts/generate_formulas.py --elements Al Ti Nb --num-elements 2 3 --max-stoich 8 --out /tmp/cand.csv
python scripts/generate_formulas.py --no-alloys                 # disable alloy pathway
python scripts/generate_formulas.py --no-electronegativity-filter
```

### `scripts/find_novel_formulas.py`
Screen candidates against a known list.

```bash
python scripts/find_novel_formulas.py --candidates /tmp/cand.csv --known data/formulas/known_formulas.csv \
    --out /tmp/novel.csv --novel-only --with-distance
```

### `scripts/search_known_formulas.py`
Ask live databases whether a formula exists (no download).

```bash
python scripts/search_known_formulas.py /tmp/novel.csv --output /tmp/live.csv
```

## 6. Data files

| File | Role |
|------|------|
| `data/formulas/known_formulas.csv` | reference formulas to subtract (your known set) |
| `data/formulas/candidates.csv` | SMACT output |
| `data/formulas/novel_1000.csv` | example: 1000 formulas absent from `known_formulas.csv` |

## 7. Source modules

- `formula_utils.py` - canonicalisation, elements, composition vectors
- `novelty/` - `KnownIndex`, novelty flags, distance ranking
- `smactlib/` - thin SMACT wrapper (namespaced to avoid shadowing the `smact` lib)
- `database_utils.py` - CSV I/O + `NetClient` HTTP helper
- `structure_utils.py` - structure metadata/comparison helpers
- `validation/` - CrystalFNO interface
- `config.py`, `common.py` - config loading and logging

## 8. How novelty is decided

A candidate is *database-unseen* relative to your reference list:

- `novel_exact_formula` - reduced formula absent from every known source;
- `novel_element_combination` - that element set never co-occurred;
- `novel_anonymous_formula` - that anonymous stoichiometry never appeared.

These are conservative labels, **not** claims of new material.

## 9. Troubleshooting

- **`ModuleNotFoundError`** - not installed: `pip install -e ".[test]"` (scripts add
  `src/` to `sys.path`, so running from the repo root works too).
- **Slow/big runs** - pass a small element universe and `--max-stoich`.
- **Live lookup errors** - COD/OQMD/JARVIS throttle sustained use; retry later or
  disable a source. Per-formula results use `error` instead of failing the run.
- **No novelty flagged** - make sure `known_formulas.csv` exists and has a
  `formula` column.

## 10. Building a known-formula DB from local datasets

`scripts/extract_known_formulas.py` turns the locally downloaded datasets into a
known-formula reference **without any database server** (a custom streaming
extractor). Symlink or point at the datasets via `datasets/COD-Dataset` and
`datasets/OQMD-Dataset` (both ignored by git).

```bash
# COD: walk the CIF mirror, read _chemical_formula_sum from each entry
python scripts/extract_known_formulas.py cod \
    --cod-dir datasets/COD-Dataset/cif --output data/formulas/known_cod.csv

# OQMD: stream the MySQL dump and pull compositions.formula (no MySQL needed)
python scripts/extract_known_formulas.py oqmd \
    --dump datasets/OQMD-Dataset/qmdb__v1_8__022026.sql/qmdb_20260129_release.sql \
    --output data/formulas/known_oqmd.csv

# merge into the reference used by find_novel_formulas.py
python scripts/extract_known_formulas.py merge \
    --inputs data/formulas/known_cod.csv data/formulas/known_oqmd.csv \
    --output data/formulas/known_formulas.csv
```

Both extractors canonicalise to reduced formulas and write
`formula,database,database_id,anonymous_formula,elements`.

- COD reads only each CIF's header (`--header-bytes`, default 16 KB) so a full
  112 GB mirror can be scanned without parsing structures.
- OQMD streams the dump line-by-line via a small tuple parser, so the ~100 GB
  MySQL import (and a server) is not required; `.sql` and `.sql.gz` both work.
- `--limit` scans only N records (fast testing); `--log-every` controls progress.

After extraction, re-run `find_novel_formulas.py` to get novelty relative to the
real COD/OQMD formulas:

```bash
python scripts/find_novel_formulas.py --candidates data/formulas/candidates.csv \
    --known data/formulas/known_formulas.csv --out data/formulas/novel_formulas.csv --novel-only
```
