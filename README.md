# Chemical Formula Parser

A Python toolkit for **parsing, generating and screening chemical formulas**.
It enumerates chemically plausible candidate compositions, screens them for
validity with SMACT, filters against known databases (COD, MP, OQMD) to find
database-unseen (novel) formulas, and selects maximally diverse subsets.

This is a *download-free* workflow. The heavy bulk-database download scripts
(COD/OQMD/JARVIS/AFLOW full collections) were removed; instead:

- **`generate_formulas.py`** — enumerate SMACT-valid candidate compositions over an
  element universe and write `data/formulas/candidates.csv`;
- **`merge_known_formulas.py`** — merge the local COD/MP/OQMD formula lists into a
  single `datasets/known_formulas.csv` holding every unique formula exactly once;
- **`find_novel_formulas.py`** — flag which candidates are absent from a supplied
  *known-formulas* list (`data/formulas/known_formulas.csv`) and rank them;
- **`search_known_formulas.py`** — ask the live databases (COD / OQMD / JARVIS)
  whether a proposed formula already exists, without downloading anything.

> **Terminology.** Candidates produced here are **database-unseen compositions** or
> **database- and literature-unreported candidates**. This code never claims "new
> material", "stable phase", or "discovered compound" — those require experiments
> and/or convex-hull/DFT validation (see *Scientific caveats*).

---

## 1. Pipeline overview

```
    element universe + SMACT settings (config.yaml)
                       |
                       v
        generate_formulas.py  ---------------> data/formulas/candidates.csv
                       |                                  |
                       v                                  |
   datasets/known_formulas.csv (unique COD+MP+OQMD         |
   formulas, one row each - build it with                  |
   scripts/merge_known_formulas.py)                        |
                       |                                  |
                       v                                  v
        find_novel_formulas.py  ---------------->  search_known_formulas.py
                       |                          (live COD/OQMD/JARVIS lookup)
                       v
        data/formulas/novel_formulas.csv
```

**What is the *known* list?** There are two ways to supply it:

1. **Merged local datasets (default for the generator scripts).** If you have the
   local COD / MP / OQMD formula extracts, `scripts/merge_known_formulas.py`
   combines `datasets/{cod,mp,oqmd}_formulas.csv` into
   `datasets/known_formulas.csv`, which contains **every unique formula exactly
   once** (columns `formula, databases`). The `databases` column keeps track of
   which sources contain the formula, so nothing is lost by merging. Candidate
   generation with `--filter-known` then drops every formula present in that file.
2. **Your own reference list.** `--known` / `--known-csv` accept any CSV with a
   `formula` column — e.g. an export from COD/OQMD/JARVIS or a curated file. The
   legacy per-database schema (`formula, database, anonymous_formula, elements`)
   is understood too, and is what gives the `found_in_<db>` /
   `novel_anonymous_formula` flags their per-database meaning.

Without a known list, novelty flags still report element/anonymous novelty but
there is no reference to subtract from.

---

## 2. Installation

Python 3.11 is the target (tested against 3.11/3.13).

**conda / mamba**

```bash
conda env create -f environment.yml
conda activate chemical-formula-parser
```

**pip**

```bash
python -m pip install -r requirements.txt
pip install -e ".[test]"        # from the project root, editable install
```

Run the tests:

```bash
python -m pytest
```

## 3. Configuration

`config.yaml` holds the element universe, SMACT parameters, filters and network
settings. Nothing is hard-coded in source; the top-level `config` module loads it
(`CHEMICAL_FORMULA_PARSER_CONFIG` overrides the path).

```yaml
elements: [Al, Cr, Co, Mo, Nb, Ta, Ti, V, W, Zr]    # default refractory space
smact: { num_elements: [2, 3, 4], max_stoich: 8, oxidation_states: icsd24,
         charge_neutrality: true, electronegativity: false, include_alloys: true }
filters: { min_elements: 2, max_elements: 6, max_atoms: 200 }
databases: { cod: true, oqmd: true, jarvis: true }  # used only by the live lookup
network: { retries: 5, backoff_factor: 1.5, timeout_seconds: 60 }
```

`databases.*` and `network.endpoints.*` are used by `search_known_formulas.py`
(the live lookup). They do **not** trigger any bulk download.

## 4. Generate candidates (SMACT)

```bash
python scripts/generate_formulas.py --num-elements 3 --max-stoich 8
python scripts/generate_formulas.py --num-elements 2 3 4 --max-stoich 8
python scripts/generate_formulas.py --elements Al Cr Co Mo Nb Ta Ti V W Zr --num-elements 3 4 --max-stoich 8
python scripts/generate_formulas.py --no-electronegativity-filter
```

Writes `data/formulas/candidates.csv` with one row per candidate reduced formula
plus SMACT filter flags (`smact_valid`, `charge_neutral`, `electronegativity_valid`).

Two enumeration regimes (driven by SMACT building blocks):

- **all-metal combinations** (`include_alloys`, default on): SMACT treats any
  all-metal composition as an alloy, so integer stoichiometry vectors up to
  `max_stoich` are enumerated (this is why the default refractory universe yields
  candidates — pure cation pairs have no charge-neutral oxidation balance);
- **combinations containing an anion**: charge-neutral enumeration over the
  chosen oxidation-state set, with the optional Pauling test.

> Enumerating *every* 3-4 element combination of the default universe can be
> large. Use `--elements`, `--num-elements` and `--max-stoich` to scope a run.

### 4.1 Generate only database-unseen candidates (`--filter-known`)

The example scripts (`examples/<system>/generate_candidates.py` and
`examples/generate_candidates.py`) wrap generation and novelty screening into one
command: they enumerate SMACT-valid candidates, drop every formula that already
appears in the merged known-formula list and return a maximally diverse subset.

```bash
cd examples/double_perovskites
python generate_candidates.py -n 1000 --filter-known -o novel_1000.txt
python generate_candidates.py -n 500 --filter-known --known-csv ../../datasets/known_formulas.csv
python generate_candidates.py -n 500 --filter-known --db cod mp   # per-database subset
```

- `--filter-known` — exclude formulas already in the known reference; this picks
  up `datasets/known_formulas.csv` automatically when it exists;
- `--known-csv PATH` — filter against a specific known-formulas CSV instead
  (any file with a `formula` column; takes precedence over `--db`);
- `--db cod mp` — fall back to the individual `datasets/{cod,mp,oqmd}_formulas.csv`
  files and use only those databases;
- `--no-diverse` — keep enumeration order instead of a spread-out subset.

## 5. Screen for novelty

### 5.1 Build the unique-formula reference (once)

```bash
python scripts/merge_known_formulas.py                # all: cod + mp + oqmd
python scripts/merge_known_formulas.py --db cod mp    # subset
```

This reads `datasets/{cod,mp,oqmd}_formulas.csv` (produced by
`scripts/extract_database_formulas.py`) and writes `datasets/known_formulas.csv`
with **every unique formula exactly once**:

```
formula,databases,anonymous_formula,elements,num_elements
Ac,cod;mp;oqmd,A,Ac,1
AcAg2As2,cod;mp,ABC3,Ac Ag As,3
NbAl2Cr,cod;mp;oqmd,ABC2,Al Cr Nb,3
```

- `databases` records which sources contain the formula, so merging loses no
  provenance (a formula present in all three databases is still a single row);
- `anonymous_formula` / `elements` are derived once (pymatgen) and make the
  element-combination and anonymous-stoichiometry novelty levels meaningful; pass
  `--no-metadata` for a faster formula-only build;
- formulas are stored in the flattened notation the generator emits
  (`Ac(AgAs)2` -> `AcAg2As2`), so membership tests compare like with like instead
  of matching strings literally;
- the file lives in the git-ignored `datasets/` folder and is rebuilt with one
  command.

### 5.2 Screen a candidates file

```bash
python scripts/find_novel_formulas.py                                    # datasets/known_formulas.csv
python scripts/find_novel_formulas.py --known datasets/known_formulas.csv
python scripts/find_novel_formulas.py --novel-only
python scripts/find_novel_formulas.py --with-distance
```

Input: SMACT candidates + a known-formulas CSV (`--known` accepts either the merged
`formula, databases` file or the per-database schema). Output:
`data/formulas/novel_formulas.csv` with per-row flags, e.g.
`found_in_cod`, `novel_exact_formula`, `novel_element_combination`,
`novel_anonymous_formula`, plus SMACT flags and (with `--with-distance`) the
nearest-known composition distance.

Three distinct novelty **levels** are reported separately and never conflated:

1. **`novel_exact_formula`** — reduced formula absent from the known list;
2. **`novel_element_combination`** — that *set of elements* never co-occurred;
3. **`novel_anonymous_formula`** — that anonymous stoichiometry (e.g. `ABC2`)
   never appeared.

`novel_exact_formula` is the key flag: **true only if the canonical reduced
formula is absent from every known source.** Labels are conservative — see the
caveats.

## 6. Live database lookup (no download)

```bash
python scripts/search_known_formulas.py candidates.txt
python scripts/search_known_formulas.py candidates.csv --output results.csv
python scripts/search_known_formulas.py in.csv --databases cod,oqmd
```

Reads proposed formulas from a text file (one per line) or a CSV (auto-detects a
`formula`/`reduced_formula`/`candidate` column; override with `--formula-column`),
canonicalises each, then asks the enabled databases whether that reduced formula
already exists via the OPTIMADE `chemical_formula_reduced` filter. Reports per
formula `KNOWN in <db>`, `not found`, or per-database `error` if a source's query
failed — one flaky database never aborts the batch. Optionally writes a CSV.

This needs **no bulk download**. Note it only reflects *reachable databases* at
query time and does not prove novelty.

## 7. Formula canonicalisation

`src/formula_utils.py` canonicalises formulas through pymatgen. Because
`reduced_formula` does not reduce fractional compositions in recent pymatgen,
canonicalisation uses `get_integer_formula_and_factor`:

```python
normalize_formula("Al2CrNb")                    # 'NbAl2Cr'
normalize_formula("CrAl2Nb")                    # 'NbAl2Cr'
normalize_formula("Cr0.5Al1.5")                 # 'Al3Cr'
get_anonymous_formula("Al2CrNb") == get_anonymous_formula("Ti2VTa") == "ABC2"
```

Malformed formulas raise `FormulaError` (or return `None` via `safe_*` helpers);
an invalid formula is never silently converted into a valid one. 118-D
composition vectors, element sets and the `composition_distance` /
`element_set_distance` baselines (for "near-known but unseen" screening) are
provided here too.

## 8. Structure-prediction interface

The package ships a small abstract interface for downstream structure-prediction
models. Implement it and feed it your novel candidates:

```python
from validation import StructurePredictorInterface, load_formulas

class MyModel(StructurePredictorInterface):
    def predict(self, formula):
        ...

model = MyModel()
for formula in load_formulas("data/formulas/novel_formulas.csv", as_dataframe=False):
    predicted = model.predict(formula)
```

`validation/dataset.py` also exposes `to_formula_list`.

## 9. Project structure (abridged)

```
config.yaml  pyproject.toml  requirements.txt  environment.yml  README.md
src/                          # top-level modules + subpackages (import directly)
  config.py formula_utils.py structure_utils.py database_utils.py common.py
  smactlib/    generator.py filters.py      # SMACT wrapper (avoids shadowing the smact lib)
  novelty/     canonicalize.py filter.py
  validation/  dataset.py
scripts/     generate_formulas.py find_novel_formulas.py search_known_formulas.py
             extract_database_formulas.py merge_known_formulas.py (unique-formula CSV)
src/databases/ (COD/OQMD/JARVIS/AFLOW clients) and bulk-download scripts were removed.
tests/       test_formula_utils.py test_structure_utils.py test_smact.py test_novelty.py
data/formulas/{candidates,known_formulas,novel_formulas}.csv
datasets/{cod,mp,oqmd}_formulas.csv  ->  datasets/known_formulas.csv (unique union)
logs/
```

## 10. Scientific caveats (important)

The code deliberately **never** asserts any of these:

1. "Not in MP" / "not in COD" ⇒ a new material — it only means *database-unseen*
   within the searched sources.
2. "New formula" ⇒ a new crystal structure (many structures share a formula).
3. "SMACT-valid" ⇒ thermodynamically stable — SMACT is a chemical-plausibility
   filter only.
4. "Low predicted energy" ⇒ stable — requires a valid reference / convex-hull
   calculation against the correct competing phases.

Use the vocabulary *database-unseen composition*, *candidate composition*,
*predicted structure*, and *candidate material*. The pipeline stops at candidate
generation; relaxation and hull-based stability screening are downstream,
out of scope, and deliberately not faked here.

## 11. Known limitations

- Novelty flags are only as good as the supplied known-formula reference (e.g. the
  merged `datasets/known_formulas.csv`) and, for the live lookup, the databases
  reachable at query time.
- Live per-formula lookups hit COD / OQMD / JARVIS one query at a time; these
  servers throttle sustained use and can be slow or unreachable from some networks.
- Composition-distance ranking uses the simple 118-D Euclidean baseline — treat it
  as a heuristic, not a learned descriptor.
