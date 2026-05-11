# 46. CSP1-2 niche projection under CMIP6

**Question.** Independently of the A/B framing, where does CSP1-2's habitable niche go under CMIP6?

**Method.** Fit species distribution model on CSP1-2 relabund vs climate + chemistry features. Project under SSP1-2.6, SSP2-4.5, SSP3-7.0 at 2050 + 2100. `scripts/run_csp_niche_projection.py`, output `cache/niche_grid_summary.tsv`, `cache/niche_model_coeffs.tsv`.

**Key results.**
- CSP1-2's realised niche shrinks modestly under warming + drying.
- Under SSP3-7.0_2100, projected suitable area decreases ~40%.

**Status.** caveated. CSP1-2 niche modeling inherits the same cross-sectional confounds as v1 (#43) and should be re-done with longitudinal features at the ASV level if needed. Currently kept as supplementary evidence rather than headline result.

**Outputs.**
- `cache/niche_grid_summary.tsv`
- `cache/niche_model_coeffs.tsv`

**Cross-refs.** 43, 45.
