# 56. Canonical data inventory

Where the load-bearing tables live, in order of how often downstream analyses touch them.

## Core feature tables

| Path | Description |
|---|---|
| `cache/feature_table.parquet` | ASV × sample counts (canonical, pre-filter) |
| `cache/feature_table_alive.parquet` | ASV × sample counts, alive-only (score ≤ 0.3) |
| `cache/feature_table_relic.parquet` | ASV × sample counts, relic-only |
| `cache/feature_table_bnti.tsv` | Subset used for iCAMP βNTI |
| `cache/taxonomy.parquet` | ASV → SILVA taxonomy |
| `cache/metadata.parquet` / `cache/metadata.tsv` | Sample-level metadata |
| `cache/metadata_with_rainfall.parquet` | Metadata + recent-precip windows |

## Per-ASV scores

| Path | Description |
|---|---|
| `cache/relic_priors/relic_score_with_mag_prior.tsv` | **Canonical** per-ASV relic-likelihood with MAG prior (#27) |
| `cache/relic_priors/relic_score.tsv` | Initial Path B composite (#25) |

## Geometry / distance

| Path | Description |
|---|---|
| `cache/distance_bray.parquet` | Pairwise BC distances |
| `cache/distance_aitchison.parquet` | Pairwise Aitchison |
| `cache/pairwise_geometry.tsv` | Geographic distance, bearing per pair |

## XRF geochemistry

| Path | Description |
|---|---|
| `cache/xrf_summary_all_trips.tsv` | Aggregate XRF table (all 5 trips) |
| `cache/xrf_lithology_pca.tsv` | Lithology PCA loadings |
| `cache/xrf_chemodiversity.tsv` | Per-element Shannon |
| `cache/xrf_per_compartment.tsv` | Compartment-level summaries |
| `cache/xrf_site_compartment_panel.tsv` | (site, comp) panel |
| `cache/per_element_shannon.tsv` | Element-level Shannon |

## Two-strategy (the late narrative)

| Path | Description |
|---|---|
| `cache/two_strategy_temporal/per_sample_strategy_with_precip.tsv` | Per-sample log₂(A/B) + precip windows |
| `cache/two_strategy_scrutiny/` | CLR analyses, mechanism counts |
| `cache/transition_asymmetry/all_transitions.tsv` | 453 trip-to-trip transitions |
| `cache/transition_asymmetry/per_cell_sequences.tsv` | Per-(site, comp) trip sequences |
| `cache/network_A_vs_B/` | A-only and B-only stratified networks |

## Climate projection

| Path | Description |
|---|---|
| `cache/two_strategy_projection/scenario_summary.tsv` | v1 (artifact) |
| `cache/two_strategy_projection_v2/scenario_summary_v2.tsv` | v2 (longitudinal, precip-only) |
| `cache/two_strategy_projection_v3/scenario_summary_v3.tsv` | **v3 (canonical, longitudinal + T)** |
| `cache/two_strategy_projection_v3/decomposition_ssp370_2100.tsv` | T-only / P-only / combined |
| `cache/per_trip_site_temperature.tsv` | Per-(site, trip) T_d30/90/365 from NASA POWER |
| `cache/cmip6_interventions.tsv` | CMIP6 ΔT/ΔP_pct deltas |

## Metagenomic guild (betA + Nibribacter)

| Path | Description |
|---|---|
| `cache/betA_guild_census.tsv` | 152 producer bins |
| `cache/betA_producers_per_bin.tsv` | Per-bin producer flag + uptake flag |
| `cache/betA_per_sample_summary.tsv` | Per-sample producer counts |
| `cache/leak_asymmetry_per_bin.tsv` | Producer ∩ uptake test (54.8% fail) |
| `cache/nibribacter_mags/per_mag_ko_assignments.tsv` | Per-MAG KO calls |
| `cache/nibribacter_mags/corrected_function_summary.tsv` | Corrected KEGG counts (#36) |
| `cache/xcomparator_betA_summary.tsv` | Multi-genus betA carriage |

## Causal frame (mostly retracted)

| Path | Description |
|---|---|
| `cache/causal_frame_tier1.parquet` | Tier-1 mediation frame (88% retracted) |
| `cache/causal_tier2_*.tsv` | PC/LiNGAM/FCI edge tables |
| `cache/causal_tier3_interventions*.tsv` | Counterfactual intervention tables |

## iCAMP

| Path | Description |
|---|---|
| `cache/icamp/process_summary_all.tsv` | Across-compartment process mix |
| `cache/icamp/RCbray_*.parquet` | Per-compartment RCbray |
| `cache/bnti/*` | βNTI intermediate |

## Sample naming

| Helper | Description |
|---|---|
| `scripts/_sample_parse.py` | Parses sample names with trip prefixes T/F/S/V; **use everywhere** |

**Status.** This is the file-finder for new authors / re-runs.

**Cross-refs.** Every other entry.
