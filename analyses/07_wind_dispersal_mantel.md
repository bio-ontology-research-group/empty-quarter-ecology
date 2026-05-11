# 07. Wind-dispersal Mantel

**Question.** If iCAMP says ~67% homogenizing dispersal, what is the physical vector? Test wind.

**Method.** For each (compartment, trip): build a wind-connectivity matrix (asymmetric, directional, from NASA POWER hourly wind over a defined backward window) and run **partial Mantel** of BC-distance against wind-distance, controlling for geographic distance. `scripts/run_wind_dispersal_mantel.py`. Per-trip per-compartment with 1, 7, 30, 90, 365-day windows.

**Inputs.**
- `cache/nasa_power_daily.parquet` (WS2M wind)
- `cache/distance_bray.parquet`
- `cache/pairwise_geometry.tsv`

**Key results.**
- Partial Mantel **r = 0.05–0.34** for BC ~ wind | geographic distance, per (compartment, trip).
- Effect **strongest at 365-day window** → annual-scale wind memory.
- Per-compartment: deep shows the largest effect, rhizosphere the smallest.

**Outputs.**
- `cache/mantel_per_trip_compartment.tsv`
- `cache/wind_dispersal/*`

**Interpretation.** Wind names the mechanism behind iCAMP's homogenizing dispersal — it operates over annual timescales (consistent with cumulative dust deposition).

**Status.** solid. Bookended with full sweep (#08).

**Cross-refs.** 06 (iCAMP gives the macroscopic statistic), 08 (sweep gives robustness).
