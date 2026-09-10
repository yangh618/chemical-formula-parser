# Tutorial: generate novel formulas

A hands-on, end-to-end walkthrough of generating candidate formulas and finding
which ones are not present in a known set.

> **Before you start** - install dependencies and make sure you are in the
> project root:
>
> ```bash
> conda env create -f environment.yml && conda activate crystalfno-discovery
> # or: python -m pip install -r requirements.txt && pip install -e ".[test]"
> ```

## Step 0 - understand the two inputs

1. **Element universe** - which elements to combine (in `config.yaml`; default is
   refractory metals `Al Cr Co Mo Nb Ta Ti V W Zr`).
2. **Known formulas** - the reference list to subtract (`data/formulas/known_formulas.csv`).
   A formula is called *novel* only when it is absent from this list.

## Step 1 - generate candidate formulas (SMACT)

Generate candidates for a small example first so it runs in seconds:

```bash
python scripts/generate_formulas.py --elements Al Ti --num-elements 2 --max-stoich 6 --out /tmp/cand.csv
```

This enumerates charge-neutral / alloy compositions of Al and Ti and writes
`/tmp/cand.csv`. Preview it:

```bash
column -s, -t /tmp/cand.csv | head
```

## Step 2 - screen for novelty

Compare the candidates against the known list:

```bash
python scripts/find_novel_formulas.py --candidates /tmp/cand.csv \
    --known data/formulas/known_formulas.csv \
    --out /tmp/novel.csv --novel-only
```

The important column is `novel_exact_formula` (`True` = absent from the known list).

## Step 3 - scale up to many formulas

Now generate over the whole element universe and a few combination sizes:

```bash
python scripts/generate_formulas.py --num-elements 2 3 --max-stoich 8 --out /tmp/cand_all.csv
python scripts/find_novel_formulas.py --candidates /tmp/cand_all.csv \
    --known data/formulas/known_formulas.csv \
    --out /tmp/novel_all.csv --novel-only
```

Keep the first 1000 novel rows (this is how `data/formulas/novel_1000.csv` was made):

```bash
head -n 1001 /tmp/novel_all.csv > data/formulas/novel_1000.csv
```

## Step 4 - one-shot convenience

Everything above is bundled in `run_formula_search.sh`:

```bash
./run_formula_search.sh "Al Ti Nb W" 3     # elements + combo size
./run_formula_search.sh                      # config.yaml universe
CHECK_LIVE=1 ./run_formula_search.sh         # also live-check vs COD/OQMD/JARVIS
```

## Step 5 - check a formula against live databases

To ask whether a proposed formula already exists in COD/OQMD/JARVIS (no bulk
download; needs network and can be slow):

```bash
python scripts/search_known_formulas.py data/formulas/novel_1000.csv --output /tmp/live.csv
```

## What the output means

`novel_exact_formula = True` means **not in the reference list you supplied**.
It is *not* a claim of a new material - a formula can exist in a database you did
not search, and many different crystal structures can share one formula.

## Keep it fast

- Fewer elements and smaller `--max-stoich` = much faster.
- 3-element combos over 10 metals at `max_stoich 8` already yields tens of
  thousands of candidates (about 12 s to generate here).

## Step 6 - Use your real COD/OQMD datasets as the "known" set

If you downloaded the COD CIF mirror and the OQMD dump, build the reference from
them (no MySQL needed) and re-screen:

```bash
python scripts/extract_known_formulas.py cod \
    --cod-dir datasets/COD-Dataset/cif --output data/formulas/known_cod.csv
python scripts/extract_known_formulas.py oqmd \
    --output data/formulas/known_oqmd.csv
python scripts/extract_known_formulas.py merge \
    --inputs data/formulas/known_cod.csv data/formulas/known_oqmd.csv \
    --output data/formulas/known_formulas.csv

python scripts/find_novel_formulas.py --candidates data/formulas/candidates.csv \
    --known data/formulas/known_formulas.csv --out data/formulas/novel_formulas.csv --novel-only
head -n 1001 data/formulas/novel_formulas.csv > data/formulas/novel_1000.csv
```

The full COD/OQMD extraction reads hundreds of GB, so run it from a normal
terminal (not a short-lived session). `--limit` scans only a sample for a quick
test first.
