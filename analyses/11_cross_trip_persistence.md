# 11. Cross-trip persistence (Tests 6 + 6A–6D)

**Question.** What fraction of OTUs persist across multiple sampling trips at the same site?

**Method.** Per (site, comp), count how many of the 5 trips each OTU appears in. Histogram. `scripts/test6_cross_trip_persistence.py`, output `cache/test6_persistence/`.

**Disconfirmation suite (6A–6D).** Persistence might be an artifact of detection thresholds, abundance bias, contamination by relic, or OTU-resolution. Tested all four: `scripts/test6a_abundance_stratified.py`, `scripts/test6b_read_floor_sensitivity.py`, `scripts/test6c_pma_validation.py`, `scripts/test6d_otu_persistence.py`. Output `cache/test6_disconfirmation/`.

**Key results.**
- **67% of OTUs are 1-trip ephemeral** (appear in only one trip at a site).
- Drops to **~57%** after collapsing ASVs to 99%-OTUs (clustering).
- Survives abundance stratification (6A).
- Survives detection-floor sensitivity (6B).
- PMA-validation (6C) inconclusive due to limited site coverage.
- OTU clustering at 97% (6D) gives an intermediate drop.

**Interpretation.** Ephemerality is **real**, not a detection or clustering artifact — but the magnitude is somewhat over-stated by ASV-level analysis (57% after OTU clustering). Consistent with dispersal-driven assembly: most lineages briefly appear, are detected, then are gone.

**Status.** solid (with magnitude caveat).

**Outputs.**
- `cache/test6_persistence/`
- `cache/test6_disconfirmation/`

**Cross-refs.** 09 (cosmopolitanism), 6 (iCAMP dispersal), 29 (alive-only refines).
