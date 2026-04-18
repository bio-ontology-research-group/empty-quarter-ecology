# `eq` — Python package for the Empty-Quarter reproducibility pipeline

Small, dependency-light package providing the core numerical and
bookkeeping primitives used by every notebook. Installed via
`make env` (`uv pip install -e .`).

## Module map

- `loader.py` — canonical readers for the feature table, taxonomy,
  metadata, XRF, climate, and iNaturalist tables. Returns
  pandas/polars frames with stable column names across every
  downstream notebook.
- `sample_id.py` — sample-ID parser that maps QIIME/Ampliseq
  identifiers (e.g. `S123_T5_rhizosphere`) to `(SiteNum, Trip,
  Compartment)` tuples; enforces integer-only site IDs in the
  range 1–60, discarding the four special-interest sites 61–64.
- `diversity.py` — Shannon, Simpson, Chao1, observed ASV, evenness.
  Operates on raw counts or relative abundance.
- `beta.py` — Bray–Curtis, Aitchison (CLR), Sørensen distance
  helpers, Procrustes, Mantel.
- `bnti.py` — vectorised β-NTI (Stegen 2012, 2013). Computes
  abundance-weighted β-MNTD on a precomputed FastTree phylogeny,
  compares to 999 tip-shuffle nulls, returns per-compartment β-NTI
  matrices.
- `core_transient.py` — Shade–Handelsman 2012 / Dini-Andreote 2015
  prevalence–abundance classifier (persistent core vs opportunist vs
  conditionally-rare). Replaces the naive detection-based rule
  flagged in the audit memory.
- `network.py` — CLR-Spearman compositional co-occurrence network;
  pseudo-FDR control; centrality ranking for keystone identification.
  Avoids SPIEC-EASI because of its memory and solver instability on
  this data size.
- `climate.py` — rolling rainfall-lag windows (7–14 d), pulse
  detection, BIE event labelling.
- `inaturalist.py` — site-to-iNaturalist observation matching (5 km
  radius, trip-window date filter); returns host-plant genus
  diversity.

## Design principles

- **Stable signatures**. Every public function returns the same
  columns regardless of which cache version is on disk; upstream
  refactors do not break downstream notebooks.
- **Vectorised where possible**. β-NTI, network edges, and mediation
  bootstrap all run end-to-end in pandas/numpy without Python-level
  loops over samples.
- **Test fixtures in `tests/`**. Run `pytest` from the repo root to
  execute the unit tests.
- **No hidden global state**. Modules read from `cache/*.parquet`
  only when called explicitly; notebooks drive the caching.
