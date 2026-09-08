# CrystalFNO External Validation & Novel-Composition Discovery

A production-oriented Python package that wraps two independent tasks around a
**frozen (pretrained) CrystalFNO** crystal-structure generation model:

- **Task A — External validation.** Collect crystal structures from databases
  *other than* the Materials Project (COD, OQMD, JARVIS, and optionally
  AFLOW), normalise them into one representation, and assemble
  *formula-unseen external tests* for CrystalFNO.
- **Task B — Novel-formula generation.** Use **SMACT** to generate chemically
  plausible candidate compositions, canonicalise them, and remove everything
  already present in the collected databases.

> **This package does NOT train CrystalFNO or any other machine-learning
> model.** CrystalFNO is only ever *consumed* through a small abstract
> interface (`CrystalFNOInterface.predict`) that you connect your trained
> model to later.

> **Terminology.** Candidates produced here are **database-unseen
> compositions** or **database- and literature-unreported candidates**. This
> code never claims "new material", "stable phase", or "discovered compound"
> — those require experiments and/or convex-hull/DFT validation (see
> *Scientific caveats*).

---

## 1. Pipeline overview

```
              EXISTING DATABASES (COD, OQMD, JARVIS, AFLOW*)
                               |
                               v
                      Unified structure DB          SMACT generation
                      (CIF + metadata CSV)                |
                               |                          v
                               v                    Chemical filtering
                    Known-formula database                 |
                               |                          v
                               +-------------------> Novelty filter
                               |                          |
                               v                          v
            Formula-unseen external tests        Novel candidate formulas
                               |                          |
                               +------------> frozen CrystalFNO ----> predicted
                               structures ---> relaxation / stability screening
```

*AFLOW is used as a **formula/metadata source only**; see section 3.*

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

`pymatgen >= 2026` changed several APIs (`Composition.reduced_formula` no
longer reduces fractional inputs; `SymmetryFinder` -> `SpacegroupAnalyzer`;
the RDF module was removed). The code here was written and tested against the
*installed* `pymatgen 2026.x`, `SMACT 4.0`, `ASE 3.29`, `spglib 2.7`
(versions checked on 2026-09-08). It uses `get_integer_formula_and_factor`
for canonicalisation so equivalent fractional formulas collapse correctly.

Run the tests:

```bash
python -m pytest
```

---

## 3. Database sources, access methods & licensing

All endpoints below were **verified against the live services** on 2026-09-08
and are configurable in `config.yaml` under `network.endpoints`.

| Database | Access used | Base endpoint (current) | Notes / licence |
|----------|-------------|--------------------------|-----------------|
| **COD** | per-id CIF + OPTIMADE metadata | `https://www.crystallography.net/cod/<id>.cif`, `.../cod/optimade/v1/structures` | Structure files are the per-record CIFs COD mirrors. Full `cod.zip`/`cod-cifs-mysql.txz` archives on `.../cod/archives/` are tens of GB; ingest a local copy with `--archive`. Open data; record provenance. |
| **OQMD** | OPTIMADE v1.2.0 | `https://oqmd.org/optimade/v1/structures` | Full structures incl. coordinates. Data CC-BY-4.0 (attribute OQMD). |
| **JARVIS-DFT** | OPTIMADE v1.1.0 | `https://jarvis.nist.gov/optimade/jarvisdft/v1/structures` | Full structures; **requires an `elements` filter** (unfiltered queries return 0 rows), so the client enumerates by element combination. NIST public-data policy. |
| **AFLOW** | AFLUX (public, no credentials) | `https://aflow.org/API/aflux/` | Used for **composition/metadata only**. Non-commercial licence; structure resolution intentionally not implemented (kept as a clean plugin point). |

**Materials Project (MP).** CrystalFNO was trained on MP. Point the pipeline at
a local MP export (`materials_project.metadata_file`) for training-overlap
removal and formula-unseen tests. The pipeline never downloads MP and works
without it.

**Licensing / redistribution (section 28).** Restricted datasets are never
redistributed. This package stores derived metadata (IDs, formulas, lattice
parameters, space groups) and per-record CIFs under the open COD/OPTIMADE
terms. `data/metadata/provenance.json` records source, URL, retrieval date,
version and counts for every database used.

## 4. Configuration

`config.yaml` holds the element universe, SMACT parameters, filters, database
toggles, MP pointer, network/timeouts and API endpoints. Nothing is
hard-coded in source; `crystalfno_discovery.config` loads it
(`CRYSTALFNO_DISCOVERY_CONFIG` overrides the path).

```yaml
elements: [Al, Cr, Co, Mo, Nb, Ta, Ti, V, W, Zr]    # default refractory space
smact: { max_stoich: 8, oxidation_states: icsd24, charge_neutrality: true, electronegativity: false }
filters: { min_elements: 2, max_elements: 6, max_atoms: 200,
           remove_partial_occupancy: true, remove_disordered: true }
databases: { cod: true, oqmd: true, jarvis: true, aflow: false }
materials_project: { metadata_file: null }
```

## 5. Downloading data

Every downloader is **resumable**: it skips records whose structure file
already exists and logs failures instead of aborting a large run.

```bash
python scripts/download_cod.py --test              # quick COD parser/connectivity test
python scripts/download_cod.py --limit 1000        # sample from COD OPTIMADE
python scripts/download_cod.py --archive data/raw/cod   # use a local COD mirror
python scripts/download_cod.py --all               # full enumeration (huge)

python scripts/download_oqmd.py --max-pages 10
python scripts/download_oqmd.py --all

python scripts/download_jarvis.py --max-pages 10   # enumerated by element combo
python scripts/download_jarvis.py --all

python scripts/download_aflow.py --max-pages 2     # composition metadata only
```

Downloads write raw payloads to `data/raw/<db>/` and normalised CIFs to
`data/structures/<db>/<id>.cif`, then rebuild `data/metadata/<db>.csv`.

## 6. Parsing structures / metadata

`python scripts/parse_structures.py --db all` scans every structure CIF and
re-writes the unified metadata CSV (idempotent, resumable-friendly). Each row
has `database, database_id, formula, reduced_formula, anonymous_formula,
elements, nelements, nsites, lattice_a..c, alpha, beta, gamma, volume,
space_group, space_group_number, structure_file` plus a **status** that tags —
never deletes — problem records:

```
valid | partial_occupancy | disordered | missing_coordinates |
invalid_lattice | parse_error | duplicate
```

Extra per-database fields (e.g. JARVIS `_jarvis_jid`) are preserved.

## 7. Formula canonicalisation

`src/crystalfno_discovery/formula_utils.py` canonicalises formulas through
pymatgen. Because `reduced_formula` does not reduce fractional compositions in
recent pymatgen, canonicalisation uses `get_integer_formula_and_factor`:

```python
normalize_formula("Al2CrNb")                 # 'NbAl2Cr'
normalize_formula("CrAl2Nb")                 # 'NbAl2Cr'
normalize_formula("Al0.6667Cr0.3333Nb0.3333")  # 'NbAl2Cr'
get_anonymous_formula("Al2CrNb") == get_anonymous_formula("Ti2VTa") == "ABC2"
```

Malformed formulas raise `FormulaError` (or return `None` via `safe_*`
helpers); an invalid formula is never silently converted into a valid one.
118-D composition vectors, element sets and the `composition_distance` /
`element_set_distance` baselines (for "near-known but unseen" screening) are
provided here too.


## 8. Building the known-formula database

`python scripts/build_known_formula_database.py [--mp data/metadata/mp.csv]`
reads `data/metadata/{cod,oqmd,jarvis,aflow}.csv` (whichever exist) plus an
optional MP export and writes `data/formulas/known_formulas.csv` with
`formula, database, database_id, anonymous_formula, elements`.

Provenance is preserved — a formula in two databases yields two rows (never
collapsed), because "found in X" must stay attributable per database.

## 9. SMACT generation

```bash
python scripts/generate_formulas.py --num-elements 3 --max-stoich 8
python scripts/generate_formulas.py --num-elements 2 3 4 --max-stoich 8
python scripts/generate_formulas.py --elements Al Cr Co Mo Nb Ta Ti V W Zr --num-elements 3 4 --max-stoich 8
python scripts/generate_formulas.py --no-electronegativity-filter
```

The generator is built on SMACT's building blocks (`smact.neutral_ratios`,
`smact.screening.pauling_test`, SMACT oxidation-state sets). Two regimes:

- **all-metal combinations** (`include_alloys`, default on): SMACT treats any
  all-metal composition as an alloy, so integer stoichiometry vectors up to
  `max_stoich` are enumerated (this is why the default refractory element
  universe produces candidates — pure cation pairs have no charge-neutral
  oxidation balance);
- **combinations containing an anion**: charge-neutral enumeration over the
  chosen oxidation-state set, with the optional Pauling test.

Output: `data/formulas/candidates.csv` with `smact_valid`,
`charge_neutral`, `electronegativity_valid` and which filters were applied.
`--no-electronegativity-filter` and the config `smact.electronegativity`
toggle these independently (not all filters are mandatory).

> Enumerating *every* 3-4 element combination of the default universe can be
> large. Use `--elements`, `--num-elements` and `--max-stoich` to scope a run.

## 10. Novelty filtering

```bash
python scripts/find_novel_formulas.py
python scripts/find_novel_formulas.py --novel-only
python scripts/find_novel_formulas.py --with-distance
```

Input: SMACT candidates + `known_formulas.csv`. Output:
`data/formulas/novel_formulas.csv` with per-row flags, e.g.
`found_in_cod`, `found_in_mp`, `novel_exact_formula`,
`novel_element_combination`, `novel_anonymous_formula`, plus SMACT filter
flags and (with `--with-distance`) the nearest-known composition distance.

Three distinct novelty **levels** (section 16) are reported separately and
never conflated:

1. **`novel_exact_formula`** — reduced formula absent from *all* searched DBs;
2. **`novel_element_combination`** — that *set of elements* never co-occurred;
3. **`novel_anonymous_formula`** — that anonymous stoichiometry (e.g.
   `ABC2`) never appeared.

`novel_exact_formula` is the key flag: **true only if the canonical reduced
formula is absent from every searched database.** Labels are conservative —
see the caveats.

## 11. External benchmark generation

```bash
python scripts/build_external_benchmark.py \
    --exclude-mp --exclude-disordered --min-elements 2 --max-elements 6
```

Produces `data/benchmark/{cod,oqmd,jarvis}_external_test.csv` and a
`summary.csv`. With `--exclude-mp` (using the MP set from the known database
or `--mp FILE`) every retained structure has a formula absent from the MP
training set — a true *formula-unseen* external test for the frozen model.
Disordered/partially-occupied/parse-failed structures can be excluded; no
structure is deleted from the source metadata.

## 12. Connecting CrystalFNO

The package never trains a model. Implement the small interface and feed it
candidates:

```python
from crystalfno_discovery.validation import CrystalFNOInterface, load_formulas

class MyFNO(CrystalFNOInterface):
    def predict(self, formula):
        # ... run your frozen CrystalFNO on the canonical formula ...
        return structure   # e.g. a pymatgen Structure

model = MyFNO()
for formula in load_formulas("data/formulas/novel_formulas.csv", as_dataframe=False):
    predicted = model.predict(formula)
    # validate / compare with structure_utils
```

`validation/dataset.py` also exposes `to_formula_list`.

## 13. Validation & structure utilities

`structure_utils` provides `validate_structure` (lattice/volume/coordinates/
composition/min-distance/overlap/symmetry) returning a structured dict, plus
`compare_structures` for predicted-vs-reference comparisons (lattice-parameter,
volume, density and space-group errors, composition match, RDF similarity and a
`StructureMatcher` match). The old pymatgen RDF module was removed upstream, so
a small periodic RDF is implemented here.

## 14. Logging & resumability

Every script logs to console **and** to `logs/`:
`download_cod.log`, `download_oqmd.log`, `download_jarvis.log`,
`generate_formulas.log`, `novelty_filter.log` (etc). Network errors are retried
with back-off and never abort a multi-record run unless you ask for it.

Resumability is file-based: a record whose `data/structures/<db>/<id>.cif`
already exists is skipped on restart; metadata CSVs are rebuilt by scanning the
structure directory, so re-running `parse_structures.py` is safe and
idempotent.

## 15. Provenance

Every downloader appends a `data/metadata/provenance.json` entry recording
source name, URL/endpoint, retrieval date (UTC), database version and record
count. Preserve database IDs in metadata; never redistribute restricted data.

## 16. Project structure (abridged)

```
config.yaml  pyproject.toml  requirements.txt  environment.yml  README.md
src/crystalfno_discovery/
  config.py formula_utils.py structure_utils.py database_utils.py common.py
  databases/   cod.py oqmd.py jarvis.py aflow.py
  smact/       generator.py filters.py
  novelty/     canonicalize.py filter.py
  validation/  dataset.py
scripts/     (9 CLI scripts)
tests/       test_formula_utils.py test_structure_utils.py test_smact.py
             test_novelty.py test_database_utils.py
data/{raw,structures,metadata,formulas,benchmark}
logs/
```

## 17. Scientific caveats (important)

The code deliberately **never** asserts any of these:

1. "Not in MP" or "not in COD" ⇒ a new material — it only means
   *database-unseen* within the searched sources.
2. "New formula" ⇒ a new crystal structure (polymorphs / many structures share
   a formula).
3. "SMACT-valid" ⇒ thermodynamically stable — SMACT is a chemical-plausibility
   filter only.
4. "CrystalFNO predicts a structure" ⇒ experimentally realisable.
5. "Low predicted energy" ⇒ stable — this requires a valid reference /
   convex-hull calculation against the correct competing phases.

Use the vocabulary *database-unseen composition*, *candidate composition*,
*predicted structure*, and *candidate material*. The pipeline stops at
candidate + predicted structures; relaxation and hull-based stability
screening are downstream, out of scope, and deliberately not faked here.

## 18. Known limitations

- **AFLOW** is a formula/metadata source only (structure resolution not
  implemented; see `databases/aflow.py`). It is disabled by default.
- **JARVIS full extraction** requires an element-by-element OPTIMADE sweep
  over the chosen element universe (unfiltered queries return nothing).
- **COD full mirror** is best done out-of-band from the multi-GB archives and
  ingested via `--archive`.
- Composition-distance ranking uses the simple 118-D Euclidean baseline from
  section 18 — treat it as a heuristic, not a learned descriptor.

