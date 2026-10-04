# Conditional spatial-covariance sensitivity

This analysis evaluates the central route-composition test on the existing
60 sites, 630 campaign–site–position groups and 200 selected genera. The original
40.0734% route fit remains descriptive. The new probabilities condition on
specified spatial covariance models and use a propagated mean design.

## Reproduce the released analysis

Run from the ecology repository root with its scientific Python environment:

```sh
python analysis/v3/spatial_covariance_sensitivity.py \
  --output-dir analysis/v3/spatial_covariance_sensitivity_20260909
python analysis/v3/render_spatial_covariance_tex.py \
  --results analysis/v3/spatial_covariance_sensitivity_20260909 \
  --output /path/to/ecology-manuscript/generated/spatial_covariance_tables.tex
python -m pytest -q tests/test_spatial_covariance_sensitivity.py
```

The first command recomputes numerical outputs in the existing result directory.
To preserve the release, create a separate output directory and copy only
`frozen_design.json` and `frozen_design_addendum.json` into it, then pass that
directory to `--output-dir`. Inputs and numerical parameters must match the
frozen design exactly. The large input is
`analysis/v2/review/cache/genus_counts.tsv`; the corrected coordinate input is
`analysis/v3/spatial_turnover_rescue/results/site_coordinates.tsv`.

The released design was frozen before observed test probabilities. An initial
positive-definiteness guard stopped because centring imposes a known rank-59
support. The pre-test addendum records support whitening and the propagated
route design. No jitter, response-based covariance selection or omitted model
was introduced. `--plan-only` creates a fresh geometry/coverage plan for a new
study; it does not replace the released design or manufacture a support addendum.

## Model and calibration

For group CLR responses G, A averages groups within sites, B removes
campaign-by-position means, and Z assigns groups to sites. The observed response
is ABG, with H=ABZ and independent-group covariance V=ABB'A'. A site covariance
K contributes HKH'. Both components are scaled to mean marginal variance one
before mixing. The 72 spatial models combine three kernels, six geometry-defined
scales and four group-noise fractions; V alone is the 73rd model.

The test whitens the known covariance support and compares the two-dimensional
propagated route subspace with 9,999 uniform random subspaces. Gaussian separable
row/column covariance permits joint cross-genus dependence, including singular
column covariance. CLR transformation and genus selection are held fixed.
Each result is conditional on its covariance assumption. The finite-family
maximum is a sensitivity envelope, with no claim of validity for arbitrary
spatial errors or evidence that a particular covariance model fits the data.

Calibration uses 500 independent datasets per generator and 199 randomizations
per dataset, with independent random draws between datasets. Matched Gaussian,
off-grid, nonstationary, response-specific and heavy-tailed nulls are all retained.
The source TSV reports five reference procedures where defined, individual
replicate probabilities and exact binomial Monte Carlo intervals. The supplement
displays unrestricted site permutation, known-generating-covariance rotation and
the finite-family maximum for every generator. The heavy-tailed case has finite
covariance, but its Gaussian generating-covariance reference is omitted.

`summary.json` records the original descriptive fit and full probability range;
`run_manifest.json` records input, code and output hashes and software versions;
`SHA256SUMS` covers the analysis bundle. The renderer writes a separate provenance
manifest. The source date is a release identifier, not a preregistration claim.

Method references: Langsrud (2005), https://doi.org/10.1007/s11222-005-4789-5;
Hettegger et al. (2021), https://doi.org/10.1093/bioinformatics/btab063.
