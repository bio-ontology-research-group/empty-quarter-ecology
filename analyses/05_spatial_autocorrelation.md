# 05. Spatial autocorrelation tests

**Question.** Is community composition spatially autocorrelated beyond what within-site replicates would generate?

**Method.** Moran's I on PCoA-axes; per-compartment. `scripts/run_spatial_autocorrelation.py`, `cache/spatial_autocorrelation_tests.tsv`.

**Key results.**
- Significant positive Moran's I on dominant PCoA axes for all three compartments.
- Spatial structure is real — not artifact of pseudo-replication.
- Drives the choice to use partial Mantels (#7–8) instead of plain correlations for wind/climate analyses.

**Status.** solid.

**Cross-refs.** 07–08 (wind partial Mantels), 31 (temporal stability).
