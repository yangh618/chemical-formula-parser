# Changelog

## 0.1.0 - novel-composition discovery (download-free)

- Flattened the package: modules importable directly under `src/`; internal
  `smact` wrapper renamed `smactlib` so it does not shadow the SMACT library.
- Removed bulk database download scripts and the COD/OQMD/JARVIS/AFLOW client
  package (full collections are not needed for local novelty screening).
- Added `scripts/search_known_formulas.py`: per-formula live lookup
  (COD/OQMD/JARVIS), no bulk download.
- Made `NetClient` retry mid-stream `ChunkedEncodingError` connection drops.
- Added `run_formula_search.sh`, `manual.md`, `tutorial.md`, `FAQ.md` and a
  generated example set `data/formulas/novel_1000.csv`.
