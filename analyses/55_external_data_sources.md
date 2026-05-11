# 55. External data sources used

| Source | Purpose | Cache / status |
|---|---|---|
| **NASA POWER** (daily, per-site) | TS (skin T), UV index, shortwave radiation, evapotranspiration, GWETTOP soil moisture, surface pressure, WS2M wind | `cache/nasa_power_daily.parquet` (60 sites × 2022-01-01 to 2025-09-30); `scripts/fetch_nasa_power.py`, `scripts/fetch_power_historical.py` |
| **NASA POWER historical** (1995–2024) | Long-term climate baseline + trends | `cache/climate_historical_1995_2024.parquet`, `cache/climate_trends_per_site.tsv` |
| **SMAP** | Surface soil moisture (microwave) | `cache/smap_daily.parquet`, `cache/smap_by_sample.parquet`; `scripts/fetch_smap_moisture.py` |
| **SRTM** | Topography (elevation, slope) | `cache/srtm_topo.parquet`, `cache/topography.tsv`; `scripts/fetch_srtm_topography.py` |
| **NASA AppEEARS — MODIS NDVI** | Per-sample NDVI (vegetation proxy) | `cache/appeears_ndvi/`, `cache/modis_ndvi.parquet`, `cache/ndvi_by_sample.parquet`, `cache/ndvi_per_compartment.tsv`; `scripts/fetch_appeears_ndvi.py` |
| **NASA AppEEARS — MODIS / MERRA-2 AOD** | Dust-event identification | `cache/appeears_aod/`, `cache/merra2_aod/`; `scripts/fetch_appeears_aod.py`, `scripts/fetch_merra2_aod.py` — **pending EULA acceptance** (Task #56) |
| **Open-Meteo (legacy)** | Was used briefly for climate pulls; **deprecated** for rate limits | `cache/climate_extended_tmp/`; `scripts/fetch_openmeteo_*.py` — **do not use** |
| **iNaturalist** | Plant species per site (RA1 plant-microbe coupling) | `cache/inat_observations.tsv`, `cache/inat_plant_microbe_coupling.tsv`, `cache/plant_diversity_by_site.tsv`, `cache/plants_by_site_trip.tsv` |
| **EMP (Earth Microbiome Project)** | Global cosmopolitanism reference (#9) | `cache/emp_cosmopolitanism/` |
| **CMIP6 SSP scenarios** | ΔT, ΔP_pct per SSP for projection | `cache/cmip6_interventions.tsv`; values applied in #43, #44, #45 |
| **GeoData (T1–T5)** | Per-(site, trip) CenterDate, lat/lon, AnnualMeanTemp, AnnualTotalPrecip | `data/geodata/trip*_geodata.tsv` |
| **XRF (T1–T5)** | Soil geochemistry (725 measurements) | `data/geochemistry/xrf_lab_table_all_trips.tsv` |

**Note.** All external data lives under `cache/` (per-fetcher subdirectories) or `data/` (committed raw). NDVI mediation, lagged compartment coupling, and pulse-reserve all consume the climate stack.

**Status.** All sources solid except MERRA-2 AOD (pending EULA acceptance by user).

**Cross-refs.** 7, 8, 15, 17, 30, 45.
