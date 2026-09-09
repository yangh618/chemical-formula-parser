# CrystalFNO Novel-Composition Discovery (SMACT)

A Python toolkit for **generating chemically plausible candidate compositions** and
screening them as **database-unseen (novel) formulas**, ready to feed a frozen
(pretrained) **CrystalFNO** crystal-structure generation model.

This is a *download-free* workflow. The heavy bulk-database download scripts
(COD/OQMD/JARVIS/AFLOW full collections) were removed; instead:

- **`generate_formulas.py`** — enumerate SMACT-valid candidate compositions over an
element universe and write `data/formulas/candidates.csv`;
- **`find_novel_formulas.py`** — flag which candidates are absent from a supplied
*known-formulas* list (`data/formulas/known_formulas.csv`) and rank them;
- **`search_known_formulas.py`** — ask the live databases (COD / OQMD / JARVIS)
whether a proposed formula already exists, without downloading anything.

> **This package does NOT train CrystalFNO or any other model.** CrystalFNO is only
> ever *consumed* through a small abstract interface (`CrystalFNOInterface.predict`)
> that you connect your own trained model to later.

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
                       |
                       v
   known_formulas.csv (your reference list)     search_known_formulas.py
                       |                                  (live COD/OQMD/JARVIS lookup)
                       v
        find_novel_formulas.py
                       |
                       v
        data/formulas/novel_formulas.csv
                       |
                       v
    frozen CrystalFNO (your model) ---> predicted structures
```

**What is the *known* list?** By default the tool reads `data/formulas/known_formulas.csv`
(columns `formula, database, database_id, anonymous_formula, elements`). The original
pipeline built this from bulk database metadata downloads; those downloaders are no
longer shipped here. You can instead supply your own reference list of formulas you
want to screen against — e.g. an export from COD/OQMD/JARVIS or a curated file. Without
a known list, novelty flags still report element/anonymous novelty but there is no
reference to subtract from.

---

## 2. Installation

Python 3.11 is the target (tested against 3.11/3.13).

**conda / mamba**

```bash
conda env create -f environment.yml
conda activate crystalfno-discovery
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
(`CRYSTALFNO_DISCOVERY_CONFIG` overrides the path).

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

## 5. Screen for novelty

```bash
python scripts/find_novel_formulas.py
python scripts/find_novel_formulas.py --novel-only
python scripts/find_novel_formulas.py --with-distance
```

Input: SMACT candidates + a known-formulas CSV. Output:
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

## 8. Connecting CrystalFNO

The package never trains a model. Implement the small interface and feed it
your novel candidates:

```python
from validation import CrystalFNOInterface, load_formulas

class MyFNO(CrystalFNOInterface):
    def predict(self, formula):
        ...

model = MyFNO()
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
src/databases/ (COD/OQMD/JARVIS/AFLOW clients) and bulk-download scripts were removed.
tests/       test_formula_utils.py test_structure_utils.py test_smact.py test_novelty.py
data/formulas/{candidates,known_formulas,novel_formulas}.csv
logs/
```

## 10. Scientific caveats (important)

The code deliberately **never** asserts any of these:

1. "Not in MP" / "not in COD" ⇒ a new material — it only means *database-unseen*
   within the searched sources.
2. "New formula" ⇒ a new crystal structure (many structures share a formula).
3. "SMACT-valid" ⇒ thermodynamically stable — SMACT is a chemical-plausibility
   filter only.
4. "CrystalFNO predicts a structure" ⇒ experimentally realisable.
5. "Low predicted energy" ⇒ stable — requires a valid reference / convex-hull
   calculation against the correct competing phases.

Use the vocabulary *database-unseen composition*, *candidate composition*,
*predicted structure*, and *candidate material*. The pipeline stops at candidate +
predicted structures; relaxation and hull-based stability screening are downstream,
out of scope, and deliberately not faked here.

## 11. Known limitations

- Novelty flags are only as good as the supplied `known_formulas.csv` and, for the
  live lookup, the databases reachable at query time.
- Live per-formula lookups hit COD / OQMD / JARVIS one query at a time; these
  servers throttle sustained use and can be slow or unreachable from some networks.
- Composition-distance ranking uses the simple 118-D Euclidean baseline — treat it
  as a heuristic, not a learned descriptor.
