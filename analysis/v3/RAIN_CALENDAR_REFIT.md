# Calendar-preserving rainfall search and complete nuisance refits

This analysis replaces the current manuscript's circular-lag search and
frozen-nuisance-residual bootstrap. Historical outputs remain archived.

Inputs are the final diversity table, corrected campaign geodata and the
two frozen daily rainfall products. No new weather retrieval is performed.
The complete input hashes and Python, NumPy, pandas, SciPy and Patsy versions
are recorded in rain_calendar_refit_20260909/manifest.json.

## Reproduce

From the ecology analysis root, use the analysis environment:

~~~bash
python analysis/v3/rain_calendar_refit.py \
  --geodata /path/to/data/metadata/geodata \
  --nasa /path/to/data/metadata/climate/nasa_power_daily_precipitation.tsv.gz \
  --open-meteo /path/to/data/metadata/climate/daily_weather_canonical.tsv \
  --output analysis/v3/rain_calendar_refit_20260909 \
  --bootstraps 9999 --simulations 300 --seed 20260910
python -m pytest -q tests/test_rain_calendar_refit.py
~~~

The root option defaults to the ecology root and supplies
analysis/v2/review/cache/alpha.tsv. All other scientific inputs are explicit.
For manuscript rendering, use the data repository's pinned figure environment
(Matplotlib 3.9.4 and FreeType 2.14.3):

~~~bash
python analysis/v3/render_rain_calendar_refit.py \
  --source analysis/v3/rain_calendar_refit_20260909 \
  --paper /path/to/authoritative/ecology-manuscript
~~~

This writes the vector figure, two numerical tables and an input/output
manifest. It uses the result tables directly.

## Statistical scope

- All six permutations of the three eligible annual rainfall fields enter the
  exact reference distribution. Collection month and day are preserved, and
  the 60 preceding actual calendar days supply each exposure.
- Every site and both products share a year map. Campaigns sharing a year
  remain coupled. The primary maximum covers both products, three endpoints
  and 60 peaks; product-specific maxima are additional summaries.
- The test conditions on the communities and collection design and assumes
  joint exchangeability of annual rainfall fields. The Monte Carlo calibration
  tests the ranking algorithm under a uniformly randomized year assignment.
- Whole-site multinomial resampling retains each site's observation history.
  Every draw refits the nuisance model. Weighted Frisch–Waugh–Lovell reduction
  gives the same estimates as the complete weighted design-matrix regression;
  the focused tests compare both calculations.
- Effect and explained-variation ranges retain the observed richness peak.
  Peak ranges reselect within richness; full-search summaries reselect both
  endpoint and peak separately within each product. These are percentile
  summaries conditional on the frozen weather and whole-site resampling.
- Wet-spell deletions remove maximal runs of regionally positive rainfall
  jointly across sites within a campaign lookback. They locate influential
  calendar periods. The code makes no independent-storm assumption.

The frozen annual fields permit probabilities in multiples of 1/6. The
current joint positive and absolute searches both give 2/6.
