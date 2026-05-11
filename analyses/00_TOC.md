# Empty Quarter Amplicon — Analysis Inventory

Complete catalogue of every experiment run on this dataset, organised so they can be threaded into a narrative. Each file: **Question · Method · Inputs · Key results · Interpretation · Status · Files**.

`Status` values: **solid** | **caveated** | **superseded by X** | **retracted** | **negative result** | **in-progress**.

---

## A. Foundation (the dataset itself)

| # | File | Pithy headline |
|---|---|---|
| 01 | [01_dataset_baseline.md](01_dataset_baseline.md) | 1,237 samples · 60 sites · 5 trips · 3 compartments |
| 02 | [02_xrf_geochemistry.md](02_xrf_geochemistry.md) | 725 XRF measurements; sabkha vs sandy chemistry |
| 03 | [03_permanova_habitat.md](03_permanova_habitat.md) | Compartment > trip > site for community structure |
| 04 | [04_distance_decay.md](04_distance_decay.md) | Tax vs functional decay; functional decays slower |
| 05 | [05_spatial_autocorrelation.md](05_spatial_autocorrelation.md) | Moran's I confirms spatial structure |

## B. Community assembly

| # | File | Pithy headline |
|---|---|---|
| 06 | [06_icamp_process_attribution.md](06_icamp_process_attribution.md) | ~67% homogenizing dispersal across compartments |
| 07 | [07_wind_dispersal_mantel.md](07_wind_dispersal_mantel.md) | Wind names iCAMP's homogenizing-dispersal vector |
| 08 | [08_wind_mantel_sweep.md](08_wind_mantel_sweep.md) | 11,520 partial Mantels; effect grows with window |
| 09 | [09_emp_cosmopolitanism.md](09_emp_cosmopolitanism.md) | >50% EQ taxa are EMP-cosmopolitan |
| 10 | [10_within_replicate_variance.md](10_within_replicate_variance.md) | Replicate < between-site variance (Test 3) |
| 11 | [11_cross_trip_persistence.md](11_cross_trip_persistence.md) | 67% of OTUs are 1-trip ephemeral (Test 6 + 6A–D) |

## C. Functional ecology

| # | File | Pithy headline |
|---|---|---|
| 12 | [12_functional_redundancy.md](12_functional_redundancy.md) | Allison-Martiny slope 0.21 → high redundancy (Test 4) |
| 13 | [13_osmolyte_blackqueen.md](13_osmolyte_blackqueen.md) | Uptake/biosynth ratio ~230× → Black Queen (Test 5) |
| 14 | [14_functional_icamp.md](14_functional_icamp.md) | Functional iCAMP near-uniform vs taxonomic 67% (Test 1) |
| 15 | [15_pulse_reserve.md](15_pulse_reserve.md) | Pulse-reserve precip alignment per site |
| 16 | [16_time_for_space.md](16_time_for_space.md) | Time-for-space substitution prep |
| 17 | [17_lagged_compartment_coupling.md](17_lagged_compartment_coupling.md) | Surface→deep coupling at trip-lag |
| 18 | [18_within_site_hill.md](18_within_site_hill.md) | Within-site Hill profile fits |
| 19 | [19_thermal_performance_curve.md](19_thermal_performance_curve.md) | Thermal calibration + EQ-specific curves |

## D. Causal inference (the retracted thread)

| # | File | Pithy headline |
|---|---|---|
| 20 | [20_mediation_tier1_retracted.md](20_mediation_tier1_retracted.md) | "88% mediated" retracted; recast non-linear Hill |
| 21 | [21_causal_inference_tiers.md](21_causal_inference_tiers.md) | DML + PC/LiNGAM/FCI structure-learning |
| 22 | [22_jsdm_knockout_alternatives.md](22_jsdm_knockout_alternatives.md) | JSDM perturbation: alternative knockouts |
| 23 | [23_stoichiometric_supply_demand.md](23_stoichiometric_supply_demand.md) | Element supply/demand mismatch |

## E. Relic-DNA detection (the pivot)

| # | File | Pithy headline |
|---|---|---|
| 24 | [24_pma_viability_qc.md](24_pma_viability_qc.md) | PMA per-pair ICC = 0.53 (workable) |
| 25 | [25_relic_indicator_composite.md](25_relic_indicator_composite.md) | Path B logistic + GB; AUC 0.70 → 0.785 w/ damage proxies |
| 26 | [26_mapdamage_pilot.md](26_mapdamage_pilot.md) | NEGATIVE: EQ DNA biologically intact, not aDNA |
| 27 | [27_relic_indicator_priors.md](27_relic_indicator_priors.md) | Bayesian taxonomic + MAG-presence priors |
| 28 | [28_relic_indicator_sensitivity.md](28_relic_indicator_sensitivity.md) | Cross-site AUC 0.69; false-relic ≈80% for rare specialists |

## F. Alive-only re-analyses (after relic filter)

| # | File | Pithy headline |
|---|---|---|
| 29 | [29_alive_population_analyses.md](29_alive_population_analyses.md) | 84% relic; alive = halotolerant Bacilli + Halomonas |
| 30 | [30_alive_climate_response.md](30_alive_climate_response.md) | Shannon~MAT is RELIC: alive ρ≈0 vs all ρ=−0.40 |
| 31 | [31_alive_temporal_stability.md](31_alive_temporal_stability.md) | All-pairwise BC; alive 2× more temporally stable |
| 32 | [32_alive_csp12_collapse_and_correction.md](32_alive_csp12_collapse_and_correction.md) | CSP1-2 "collapse" overturned by MAG prior; rank 6 not 1 |
| 33 | [33_alive_remaining_analyses.md](33_alive_remaining_analyses.md) | iCAMP, mediation, wind-Mantel, distance-decay on alive |

## G. Keystone hunt and Nibribacter

| # | File | Pithy headline |
|---|---|---|
| 34 | [34_keystone_hunt_alive.md](34_keystone_hunt_alive.md) | Nibribacter top cross-compartment in alive net |
| 35 | [35_nibribacter_mag_function.md](35_nibribacter_mag_function.md) | 16 MAGs; KEGG functional profile |
| 36 | [36_nibribacter_kegg_corrected.md](36_nibribacter_kegg_corrected.md) | Regex artifact corrected: trehalose, not betaine |
| 37 | [37_nibribacter_xrf_climate.md](37_nibribacter_xrf_climate.md) | Prefers sandy non-saline sites |
| 38 | [38_keystone_vs_guild_knockout.md](38_keystone_vs_guild_knockout.md) | 3-guild beats single-keystone: 7% vs 46% knockout |

## H. Two-strategy architecture

| # | File | Pithy headline |
|---|---|---|
| 39 | [39_two_strategy_temporal.md](39_two_strategy_temporal.md) | Anti-correlated A (DOM) vs B (halotolerant) strategies |
| 40 | [40_two_strategy_scrutiny.md](40_two_strategy_scrutiny.md) | CLR ρ=−0.47, p=97th pct null; B mech = sporulation+ectoine |
| 41 | [41_transition_asymmetry.md](41_transition_asymmetry.md) | A→A 88% (resilient), B→B 52% (transient) |
| 42 | [42_a_vs_b_network.md](42_a_vs_b_network.md) | B network 2.4× denser, fewer modules |

## I. Climate projection

| # | File | Pithy headline |
|---|---|---|
| 43 | [43_climate_projection_v1_artifact.md](43_climate_projection_v1_artifact.md) | Cross-sectional ARTIFACT; geographic confounding |
| 44 | [44_climate_projection_v2_longitudinal.md](44_climate_projection_v2_longitudinal.md) | Precip-only transitions; −5 pp π_B |
| 45 | [45_climate_projection_v3_with_temperature.md](45_climate_projection_v3_with_temperature.md) | T added; +50 pp π_B; warming → more B |
| 46 | [46_csp_niche_projection_cmip6.md](46_csp_niche_projection_cmip6.md) | CSP1-2 niche shrinkage under CMIP6 |

## J. Cross-cohort / external benchmarks

| # | File | Pithy headline |
|---|---|---|
| 47 | [47_cross_desert_and_atacama.md](47_cross_desert_and_atacama.md) | Gurbantunggut / Namib / McMurdo / Atacama comparisons |
| 48 | [48_atacama_within_desert.md](48_atacama_within_desert.md) | Atacama CSP1-2 is Altiplano-only |
| 49 | [49_hill_cross_cohort.md](49_hill_cross_cohort.md) | EQ Hill does NOT transfer (n_Atacama=1.16, linear) |

## K. Metagenomic guild (betA)

| # | File | Pithy headline |
|---|---|---|
| 50 | [50_betA_guild_metagenomic.md](50_betA_guild_metagenomic.md) | 296 metagenomes, 9,229 bins, 152 producer bins |
| 51 | [51_betA_leak_asymmetry.md](51_betA_leak_asymmetry.md) | FAILS: 54.8% producers also have uptake |
| 52 | [52_betA_xcomparator.md](52_betA_xcomparator.md) | Rubellimicrobium ALSO has betA (single-keystone retracted) |

## L. Referee-defense & sensitivity

| # | File | Pithy headline |
|---|---|---|
| 53 | [53_adversarial_keystone.md](53_adversarial_keystone.md) | CSP1-2 rank 5/15 by knockout magnitude |
| 54 | [54_referee_defense_suite.md](54_referee_defense_suite.md) | 13 sensitivity analyses; thermal-bound retracted |

## M. External data sources used

| # | File | Pithy headline |
|---|---|---|
| 55 | [55_external_data_sources.md](55_external_data_sources.md) | NASA POWER, SMAP, SRTM, NDVI, AOD, iNat plants |
| 56 | [56_data_inventory.md](56_data_inventory.md) | Where the canonical tables live |

## Cross-cutting narrative threads

- **The relic-DNA pivot**: 24 → 25 → 26 → 27 → 28 → 29 → 32 (CSP1-2 keystone collapse and partial rescue)
- **The two-strategy discovery**: 32 → 34 → 38 → 39 → 40 → 41 → 42 → 45
- **The climate projection saga**: 43 (artifact) → 44 (longitudinal, precip-only) → 45 (T added, direction reversed)
- **The betA arc**: 38 (3-guild beats 1-keystone) → 50 (guild census) → 51 (leak-asymmetry FAILS) → 52 (multi-genus betA)
