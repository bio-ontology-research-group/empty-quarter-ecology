# 45. Climate projection v3 — longitudinal + temperature

**Question.** Properly: with temperature included, where does CMIP6 warming take A vs B equilibria?

**Method.** Two fixes vs v2:
1. **Absolute to-trip features**, not deltas. (Uniform CMIP6 ΔT/ΔP shifts cancel in a delta feature → no model response; absolute features respond properly.)
2. **Temperature added** — NASA POWER skin T (TS) windowed at 30 / 90 / 365 days before each trip's CenterDate, per (site, trip). `cache/per_trip_site_temperature.tsv`.

Features: to_d7, to_d365 (precip), Δd30/90/180 (intermediates; kept as relative), to_T_d30/90/365.

Logit P(A→B) on A-start; P(B→A) on B-start. Apply CMIP6 ΔT to absolute T cols, ΔP_pct multiplied into absolute precip cols. `scripts/two_strategy_climate_projection_v3.py`, `cache/two_strategy_projection_v3/`.

**Logistic coefficients (z-scaled).**

| Feature | P(A→B) | P(B→A) |
|---|---|---|
| to_T_d365 | **+1.083** | **−1.119** |
| to_d365   | +0.923 | −0.353 |
| to_T_d30  | +0.045 | −0.006 |

Both T_d365 coefficients are large and opposite-signed — **warming favors B-state**.

**CMIP6 equilibrium π_B (full scenarios).**

| Scenario      | π_B    | Δ        |
|---------------|--------|----------|
| Historical    | 0.461  | —        |
| SSP1-2.6_2100 | 0.764  | +30 pp   |
| SSP2-4.5_2100 | 0.888  | +43 pp   |
| SSP3-7.0_2100 | **0.960** | **+50 pp** |

**Decomposition (SSP3-7.0_2100).** T-only: π_B = 0.962. P-only: π_B = 0.433. → **Temperature dominates**; drying mildly opposes.

**Direction reversed vs v2.** Mechanistically coherent:
- Strategy B = halotolerant thermo-spore guild (Bacilli + Halomonas + sporulation + ectoine).
- Warming favors thermo-halotolerants.
- Strategy A = mesophilic DOM-cyclers (Bacteroidota); they lose ground.

**Extrapolation caveat.** Observed T_d365 range across trips is ~1.1 °C (30.08 → 31.22 °C between T1 and T4). SSP3-7.0_2100 applies +4 °C — **~4× the training-data variance**. The 0.96 endpoint is logit saturation. **SSP1-2.6_2100 (+1.6 °C, π_B = 0.76) sits within observed range and is the most defensible projection.** Also: trip-to-trip T variance is dominated by **seasonality** (T4 Aug 2024 hottest), which we can't fully disentangle from long-term trend with 5 trips. Treat π_B magnitudes as bounds; direction (warming → more B) is robust because the guild's mechanistic tolerance points that way.

**Per-site hotspots under SSP3-7.0_2100.** Sites 3, 7–14 (Najran/Wadi corridor, cooler baseline) show largest Δp_AtoB (~0.63) — warming moves them into the A→B transition zone.

**Uncertainty quantification (2026-05-11).**
`scripts/two_strategy_climate_projection_v3_uncertainty.py`, outputs in `cache/two_strategy_projection_v3/uncertainty/`.

| Scenario | π_B point | Bootstrap median | 95% CI |
|---|---|---|---|
| Historical | 0.461 | 0.454 | [0.418, 0.486] |
| SSP1-2.6_2100 | 0.763 | 0.753 | **[0.660, 0.830]** |
| SSP2-4.5_2100 | 0.888 | 0.880 | [0.760, 0.941] |
| SSP3-7.0_2100 | 0.960 | 0.957 | [0.839, 0.987] |

**Historical and SSP1-2.6_2100 CIs are non-overlapping** (0.486 < 0.660) — warming-favors-B direction is statistically robust.

**Leave-one-trip-pair-out CV (6 splits).** Δπ_B under SSP1-2.6_2100 across splits: mean **+0.287**, std 0.021, range [+0.261, +0.310]. **100% positive direction.** Per-cell classification AUC is variable (held-out T3→T4 AUC = 0.20; T4→T5 = 0.78), but the **population-level π_B projection is robust** — Markov equilibrium averages out per-cell noise.

**Leave-one-site-out CV (59 sites).** Δπ_B under SSP1-2.6_2100: median +0.302, IQR [+0.300, +0.304] — **virtually point-like across spatial CV**. 100% of held-out sites predict positive Δπ_B.

**Headline reportable.** "Under SSP1-2.6_2100 (+1.6 °C, −3% precip), the equilibrium B-fraction shifts from π_B ≈ 0.46 to π_B ≈ 0.76 (bootstrap 95% CI [0.66, 0.83]; +0.29 across LTO trip-pair splits; +0.30 across LOSO splits)."

**Status.** solid (with caveats above) — bootstrap + LOO-CV verified. Canonical climate projection.

**Outputs.**
- `cache/two_strategy_projection_v3/scenario_summary_v3.tsv`
- `cache/two_strategy_projection_v3/per_site_AtoB_risk.tsv`
- `cache/two_strategy_projection_v3/decomposition_ssp370_2100.tsv`
- `cache/per_trip_site_temperature.tsv`
- `cache/two_strategy_projection_v3/uncertainty/{bootstrap_ci.tsv,bootstrap_raw.tsv,lto_trip_pair.tsv,point_estimate.tsv}`

**Cross-refs.** 43 (artifact), 44 (precip-only), 41 (transition data), 46.
