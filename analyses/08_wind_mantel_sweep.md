# 08. Wind-Mantel sensitivity sweep

**Question.** How robust is the wind-Mantel result to wind-score definition, distance-stratification, asymmetric vs symmetric edges, threshold choices, and time window?

**Method.** Comprehensive grid: 5 time windows × asymmetric/symmetric × distance-stratified vs pooled × multiple threshold cutoffs × multiple wind-score definitions. **11,520 partial Mantel tests in total.** `scripts/run_wind_mantel_sweep.py`, render in `scripts/render_wind_sweep.py`.

**Inputs.** Same as #07.

**Key results.**
- Effect **monotonically grows with window length** (1d → 365d) — annual memory.
- Asymmetric (true wind direction) > symmetric, but both significant.
- Distance-stratified Mantels confirm wind effect is **not just a geographic-proxy artifact** — it persists within distance bins.
- Threshold-robust across all reasonable wind-distance cutoffs.

**Outputs.** Per-cell test results in `cache/wind_dispersal/` (large grid).

**Interpretation.** The wind-dispersal signal is robust against any single parameter choice. Mechanism is reified.

**Status.** solid. AOD validation pending (NASA AppEEARS EULA — see #55).

**Cross-refs.** 07, 56.
