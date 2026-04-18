# Notebooks

Each notebook is self-contained and reproducible — given the inputs
listed, it produces the named outputs via `quarto render`. Notebooks
are numbered in execution order; later numbers depend on caches
written by earlier numbers. Run the full chain with `make figures`
from the repository root (after `make env && make cache`).

Inputs referenced as `cache/*.parquet` are written by notebook 00 and
consumed by every downstream notebook. Inputs under `data/` are
small, shipped with the repository. Inputs under `cache/` that
require external fetches are documented in `scripts/README.md`.

---

## `00_load_and_qc.qmd` — Loading and quality control

**Purpose.** Loads the 1,237-sample Ampliseq output, applies QC
filters, decontam-style removal of lab contaminants, and writes the
analysis-ready feature table.

- **Inputs:** `data/taxonomy/feature-table-trips1-5.tsv` (regenerate
  from raw reads via nf-core/Ampliseq; see Data section of the main
  README), `data/taxonomy/taxonomy-trips1-5.tsv`,
  `data/geodata/trip{1–5}_geodata.tsv`,
  `src/eq/loader.py`, `src/eq/sample_id.py`.
- **Outputs:** `cache/feature_table.parquet`,
  `cache/taxonomy.parquet`, `cache/metadata.parquet`,
  `cache/decontam_flagged.tsv`.
- **Figures produced:** none (diagnostic plots only).
- **Compute:** local, ~2 min.

## `01_scale_and_phyla.qmd` — Spatial scale and dominant phyla

**Purpose.** Reports the 1,000 km transect geometry and per-phylum
compartment signatures; produces Fig. 1 panels and the per-element
XRF–Shannon correlations (Table S1).

- **Inputs:** `cache/feature_table.parquet`,
  `cache/taxonomy.parquet`, `cache/metadata.parquet`,
  `data/geochemistry/xrf_lab_table_filtered.tsv`.
- **Outputs:** `figures/fig1a_sites.pdf`, `figures/fig1c_phyla.pdf`,
  `figures/fig1_xrf_pca_compartment.pdf`,
  `figures/fig1_shannon_surface_deep.pdf`.
- **Manuscript:** Results §1 "Scale, dominant phyla, and
  geochemical drivers of diversity".
- **Compute:** local, <1 min.

## `02_assembly.qmd` — Null-model assembly decomposition

**Purpose.** Computes β-NTI on the 500 most prevalent ASVs per
compartment using 999 tip-shuffle nulls; produces Fig. 2 panels
for compartment-specific assembly regimes.

- **Inputs:** `cache/feature_table.parquet`,
  `cache/taxonomy.parquet`, `src/eq/bnti.py` (vectorised
  implementation of Stegen 2013).
- **Outputs:** `figures/fig2ab_pcoa.pdf`,
  `figures/fig2cd_bnti.pdf`, `cache/bnti_{surface,deep,rhiz}.tsv`.
- **Manuscript:** Results §2 "Compartment-specific assembly
  regimes"; Suppl S2.
- **Compute:** local, ~5 min (999 permutations).

## `03_temporal.qmd` — Rainfall-lag dynamics

**Purpose.** Fits 7–14 d lagged rainfall→Shannon/community Bray-Curtis
panel regression; reports the core/opportunist/conditionally-rare
classification and the 408-genus conditionally rare bloomer pool
(Shade-Handelsman/Dini-Andreote framework).

- **Inputs:** `cache/feature_table.parquet`,
  `data/climate/daily_weather.tsv`.
- **Outputs:** `figures/fig3a_rainfall_lag.pdf`,
  `figures/fig3b_core_transient.pdf`,
  `cache/conditionally_rare_pool.tsv`.
- **Manuscript:** Results §3 "Pulse-reserve dynamics"; Suppl S3.
- **Compute:** local, <2 min.

## `04_depth.qmd` — Depth refugium comparison

**Purpose.** Paired surface/deep Shannon contrasts; tests the
depth-refugium hypothesis for hyperarid sand.

- **Inputs:** `cache/feature_table.parquet`,
  `cache/metadata.parquet`.
- **Outputs:** `figures/fig4a_shannon_by_trip.pdf`,
  `figures/fig4b_paired_surface_deep.pdf`,
  `figures/fig4d_host_plant.pdf`.
- **Manuscript:** Results §4 "Depth refugia"; Suppl S4.
- **Compute:** local, <1 min.

## `05_distance_decay.qmd` — Distance decay of similarity

**Purpose.** Bray–Curtis and Sørensen distance-decay slopes by
compartment and trip; Mantel tests against Haversine geographic
distance.

- **Inputs:** `cache/feature_table.parquet`,
  `cache/metadata.parquet`.
- **Outputs:** `figures/fig5a_distance_decay.pdf`,
  `figures/fig5b_mantel.pdf`, `cache/mantel_results.tsv`.
- **Manuscript:** Results §5 "Distance-decay"; Suppl S5.
- **Compute:** local, ~1 min.

## `06_function.qmd` — Functional redundancy

**Purpose.** PICRUSt2 pathway-level functional distance vs taxonomic
distance; Procrustes superposition.

- **Inputs:** `data/functional/picrust2/path_abun_unstrat.tsv`,
  `cache/feature_table.parquet`.
- **Outputs:** `figures/fig6a_tax_vs_func.pdf`,
  `figures/fig6b_decay_tax_vs_func.pdf`.
- **Manuscript:** Results §6 "Functional redundancy"; Suppl S10.
- **Compute:** local, <2 min.

## `07_network.qmd` — Co-occurrence network and keystone ranking

**Purpose.** Builds the compositional co-occurrence network via
CLR-Spearman + pseudo-FDR per compartment; ranks genus keystones by
centrality; recovers MND1–*Nitrospira* and the CSP1-2 keystone hub.

- **Inputs:** `cache/feature_table.parquet`,
  `cache/taxonomy.parquet`.
- **Outputs:** `figures/fig7_network.pdf`,
  `cache/network_{surface,deep,rhiz}_edges.tsv`,
  `cache/keystone_ranks.tsv`.
- **Manuscript:** Results §7 "Co-occurrence network and keystones".
- **Compute:** local, ~3 min.

## `08_csp_mag.qmd` — CSP1-2 MAGs and dark-matter annotation

**Purpose.** Parses CheckM2 + GTDB-Tk + gapseq outputs for the four
co-assembled CSP1-2 MAGs; reports nitrogenase, trehalose, DNA-repair
inventories; exports `cache/csp1-2_asvs.fasta` for the cross-desert
VSEARCH search.

- **Inputs:** `cache/mag_checkm2.tsv`,
  `cache/mag_gtdbtk.tsv`,
  `cache/mag_gapseq_pathways.tsv` (produced by cluster job; see
  `scripts/gspa_run.sh`).
- **Outputs:** `figures/fig8a_checkm2.pdf`,
  `figures/fig8b_functional.pdf`,
  `figures/fig8c_darkmatter.pdf`,
  `cache/csp1-2_asvs.fasta` (input for notebook 16).
- **Manuscript:** Results §8 "CSP1-2 candidate-phylum MAGs".
- **Compute:** local notebook (annotation happens on cluster);
  ~1 min.

## `09_causal_tier1.qmd` — Panel FE + DML + mediation

**Purpose.** Two-way within-transformation panel regression with
cluster-robust SEs, Chernozhukov 2018 Double/Debiased ML for XRF
ATEs on Shannon and CSP1-2, and bootstrap causal mediation
(Imai 2010) of the S → CSP1-2 → Shannon path. Reports the
88% indirect-effect share used in the abstract.

- **Inputs:** `cache/feature_table.parquet`,
  `data/geochemistry/xrf_lab_table_filtered.tsv`,
  `data/climate/daily_weather.tsv`.
- **Outputs:** `figures/fig9ab_causal_tier1.pdf`,
  `figures/figS4_mediation_sensitivity.pdf`,
  `cache/dml_ate_estimates.tsv`,
  `cache/mediation_bootstrap.tsv`.
- **Manuscript:** Results causal subsection.
- **Compute:** local, ~3 min (bootstrap with N=1,000).

## `10_causal_tier2.qmd` — PC + FCI + LiNGAM structure learning

**Purpose.** Learns causal graph structure from the geochemistry +
diversity + CSP1-2 panel. Runs 200 bootstrap iterations of PC / FCI
/ LiNGAM; reports edge-presence frequencies and consensus PAG.

- **Inputs:** `cache/mediation_input.parquet` (from notebook 09).
- **Outputs:** `figures/fig9c_causal_pag.pdf`,
  `cache/fci_edge_freq.tsv`, `cache/lingam_edge_freq.tsv`.
- **Manuscript:** Results causal subsection.
- **Compute:** local, ~10 min (200 bootstraps).

## `11_causal_tier3.qmd` — Bayesian hierarchical state-space twin

**Purpose.** NumPyro/NUTS hierarchical state-space model for
compartment × trip × site. Posterior interventional forecasts:
do(rain +10 mm), do(S −1 SD), do(S −2 SD, P +1 SD).

- **Inputs:** `cache/feature_table.parquet`,
  `data/climate/daily_weather.tsv`,
  `data/geochemistry/xrf_lab_table_filtered.tsv`.
- **Outputs:** `figures/fig9d_digital_twin.pdf`,
  `cache/numpyro_posterior.parquet`.
- **Manuscript:** Results causal subsection; abstract.
- **Compute:** local, ~12 min (4 chains × 2,000 warmup + 2,000
  draws; fits on a laptop with jax CPU backend, GPU accelerates).

## `12_causal_nonlinear.qmd` — Hill / spline / GP dose-response

**Purpose.** Non-linear rainfall→Shannon dose-response under three
functional forms: Hill/Monod (BIE-threshold, Reynolds 2004), Bayesian
cubic B-spline with shrinkage, and Matérn-3/2 Gaussian process.

- **Inputs:** `cache/feature_table.parquet`,
  `data/climate/daily_weather.tsv`.
- **Outputs:** `figures/fig9e_nonlinear_doseresp.pdf`,
  `cache/hill_params.tsv`, `cache/spline_posterior.parquet`,
  `cache/gp_posterior.parquet`.
- **Manuscript:** Results causal subsection; Suppl S19.
- **Compute:** local, ~5 min.

## `13_vegetation_mediation.qmd` — NDVI mediation + iNat plants

**Purpose.** Tests the rainfall → NDVI → Shannon mediation path
(47.5% of total effect); reports cross-kingdom plant-diversity
(iNaturalist) × microbial Shannon correlation in the rhizosphere.

- **Inputs:** `cache/metadata.parquet`,
  `cache/feature_table.parquet`,
  `cache/inat_plants.parquet`,
  `cache/modis_ndvi.parquet` (from `scripts/fetch_appeears_ndvi.py`).
- **Outputs:** `figures/fig10_vegetation_mediation.pdf`,
  `cache/ndvi_mediation.tsv`.
- **Manuscript:** Results vegetation subsection; Suppl S22, S25.
- **Compute:** local, ~2 min.

## `14_keystone_knockout.qmd` — JSDM latent-factor knockout

**Purpose.** Gradient-boosted baseline + conditional-MVN joint species
distribution model with 3 latent factors; computationally removes
CSP1-2 and reports propagated drop in Shannon and the identity of
the four most-affected genera.

- **Inputs:** `cache/feature_table.parquet`,
  `cache/taxonomy.parquet`.
- **Outputs:** `figures/fig14_keystone_knockout.pdf`,
  `cache/jsdm_knockout.tsv`.
- **Manuscript:** Suppl S20, S26.
- **Compute:** local, ~4 min.

## `15_intervention_scenarios.qmd` — Reclamation + lagged coupling

**Purpose.** Posterior draws under the do(S−2 SD, P+1 SD) reclamation
scenario; lagged cross-compartment coupling regression; topographic
integration.

- **Inputs:** `cache/numpyro_posterior.parquet`,
  `cache/feature_table.parquet`,
  `cache/srtm_topo.parquet` (from `scripts/fetch_srtm_topography.py`).
- **Outputs:** `figures/fig15_reclamation.pdf`,
  `cache/reclamation_posterior.tsv`,
  `cache/lagged_compartment_coupling.tsv`.
- **Manuscript:** Suppl S27, S28.
- **Compute:** local, ~2 min.

## `17_xrf_compartment.qmd` — Per-compartment XRF × Shannon (Fig S30)

**Purpose.** Decomposes the site-mean per-element XRF × Shannon
correlations reported in main Fig. 1C by compartment. Identifies
which salinity ions (S, Cl, Na) act in which soil compartment and
surfaces the P sign flip between deep (+) and rhizosphere (−).

- **Inputs:** `data/geochemistry/xrf_lab_table_filtered.tsv`,
  `cache/feature_table.parquet`, `cache/metadata.parquet`.
- **Outputs:** `figures/figS30_xrf_per_compartment.pdf`,
  `cache/xrf_per_compartment.tsv`.
- **Manuscript:** Suppl S30.
- **Compute:** local, <1 min.

## `16_cross_desert.qmd` — Cross-desert generalisation (Fig 11)

**Purpose.** Summarises the cross-desert comparison against published
16S surveys of Atacama (PRJEB17617), Namib (PRJNA628615), and
McMurdo (PRJNA721735). Produces Fig. 11 (Shannon distributions,
CSP1-2 prevalence, Namib abundance–Shannon, Atacama EC–Shannon).

- **Inputs:** `cache/crossdesert/per_sample.tsv`,
  `cache/crossdesert/comparison_summary.tsv`,
  `cache/crossdesert/eq_shannon_reference.tsv` (produced by
  `stage4_analyze.py` / `stage5_genus_compare.py` on the Unimatrix
  cluster; see `scripts/README.md`).
- **Outputs:** `figures/fig11_crossdesert.pdf`,
  `figures/fig11_crossdesert.png`.
- **Manuscript:** Results cross-desert subsection; Suppl S29.
- **Compute:** local, <1 min (summary only; the upstream SINTAX +
  VSEARCH clustering happens on the cluster at ~90 min of 8-thread
  wall-clock).

## `99_audit.qmd` — Cross-notebook sanity checks

**Purpose.** Re-reads the cached outputs of every earlier notebook
and checks the published numbers against the caches. Catches
regressions when upstream data, parameters, or code change.

- **Inputs:** every `cache/*.parquet` / `cache/*.tsv`.
- **Outputs:** `figures/audit_contamination.pdf`,
  `figures/audit_per_element_shannon.pdf`, stdout summary table.
- **Compute:** local, ~1 min.
