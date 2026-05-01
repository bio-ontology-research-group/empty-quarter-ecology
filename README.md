# Empty Quarter microbiome — reproducibility repository

Computational analyses and manuscript figures for:

> **A longitudinal, causally-grounded microbiome atlas of Earth's largest sand desert.**
> Tawfiq, R., Martinez Arbas, S., Abdelhakim, M., Niu, K., Lopez Velazquez, A.,
> Rueping, M., Hoehndorf, R. (in preparation).

This repository contains every notebook, script, and (small) data table
needed to reproduce the figures, tables, and statistics of the manuscript.
It is a computational companion to the semantic data descriptor submitted
in parallel:

> **A semantic knowledge graph of the Rub al-Khali metagenomic expedition.**
> (Companion data paper; in preparation.) Provides the formal knowledge
> representation of samples, sites, geochemistry, taxonomy, and
> provenance (~602,186 triples, 147 sites, 1,607 soil samples) and the
> SPARQL endpoint used as the upstream data source for the analyses
> below.

## Scope and design

The analyses span two compute contexts:

1. **Local (laptop / workstation).** QC, ordination, diversity,
   networks, PICRUSt2 interpretation, causal inference, manuscript
   figures. Python + Quarto. 16 GB RAM sufficient; no GPU required.
   Single command: `make env && make cache && make figures`.

2. **Cluster (Unimatrix / IBEX Slurm).** Compute-heavy steps —
   nf-core/Ampliseq, β-NTI null modelling, MAG recovery, GTDB-Tk
   classification, GSPA/gapseq dark-matter annotation, cross-desert 16S
   reprocessing. Pipeline scripts under `scripts/` are intended to be
   run from a cluster session; outputs are pulled back into `cache/`
   for the local notebooks.

## Notebook inventory

All notebooks are Quarto `.qmd` files under `notebooks/`; each renders
to HTML under `_output/notebooks/` via `make figures`.

| Notebook | Produces | Manuscript section |
| --- | --- | --- |
| `00_load_and_qc.qmd` | `cache/feature_table.parquet`, decontam | Methods |
| `01_scale_and_phyla.qmd` | Fig 1 (sites, phyla, XRF) | Results §1 |
| `02_assembly.qmd` | Fig 2 (PCoA, β-NTI) | Results §2 |
| `03_temporal.qmd` | Fig 3 (rainfall-lag) | Results §3 |
| `04_depth.qmd` | Fig 4 (depth refugium) | Results §4 |
| `05_distance_decay.qmd` | Fig 5 (distance decay) | Results §5 |
| `06_function.qmd` | Fig 6 (functional redundancy) | Results §6 |
| `07_network.qmd` | Fig 7 (co-occurrence network) | Results §7 |
| `08_csp_mag.qmd` | Fig 8 (CheckM2, functional, MAGs) | Results §8 |
| `09_causal_tier1.qmd` | Fig 9a,b (Panel FE + DML + mediation) | Results causal |
| `10_causal_tier2.qmd` | Fig 9c (FCI + LiNGAM structure learning) | Results causal |
| `11_causal_tier3.qmd` | Fig 9d (Bayesian state-space twin) | Results causal |
| `12_causal_nonlinear.qmd` | Fig 9e (Hill / spline / GP dose-response) | Results causal |
| `13_vegetation_mediation.qmd` | Fig 10a (NDVI mediation, iNat plants) | Results vegetation |
| `14_keystone_knockout.qmd` | Fig 14 (JSDM knockout) | Suppl S20 |
| `15_intervention_scenarios.qmd` | Fig 15 (reclamation, lagged coupling) | Suppl S27 |
| `16_cross_desert.qmd` | Fig 11 (Atacama / Namib / McMurdo) | Results §cross-desert |
| `99_audit.qmd` | Cross-notebook sanity checks | — |

## Repository layout

```
├── notebooks/         # 18 Quarto notebooks
├── src/eq/            # Python package: loading / diversity / β-NTI / network
├── scripts/           # Data fetchers (AppEEARS, SMAP, NASA POWER,
│                      # Copernicus DEM), Slurm job templates, figure builders
├── R/                 # R scripts for HMSC / joint species distribution models
├── analysis/          # Per-research-question analyses (RQ01–RQ24)
├── data/              # Small metadata and intermediate tables
├── figures/           # Rendered figure PDFs (main + supplement)
├── Makefile           # `make env | cache | figures | check | clean`
├── pyproject.toml     # Python dependency declarations
├── uv.lock            # Reproducible Python environment (uv)
├── environment.yml    # Conda-friendly alternative
├── _quarto.yml        # Quarto project config
└── README.md          # this file
```

## Data

Small metadata and derived tables are included under `data/`:

- `data/climate/daily_weather.tsv` — daily weather pulled from NASA POWER
- `data/climate/monthly_weather_averages.tsv` — monthly aggregates
- `data/geochemistry/xrf_lab_table_filtered.tsv` — Trip 5 XRF panel
  (158 samples × 33 elements)
- `data/geodata/trip{1–5}_geodata.tsv` — site coordinates + metadata
- `data/taxonomy/taxonomy-trips1-5.tsv` — ASV → lineage taxonomy
- `data/functional/picrust2/path_abun_unstrat*.tsv` — PICRUSt2 pathway
  abundances (~7 MB each)
- `data/functional/picrust2/metagenome_pred_metagenome_unstrat.tsv` —
  PICRUSt2 gene predictions (23 MB)

**Large data not included in the repository** (regenerate or fetch
externally):

- `data/taxonomy/feature-table-trips1-5.tsv` — 1,227-sample feature
  table (1.7 GB). Regenerate from raw reads via nf-core/Ampliseq or
  fetch from the companion KG SPARQL endpoint (see data paper).
- `data/taxonomy/ASV_seqs-trips1-5.fasta` — ASV representatives
  (155 MB). Same provenance as the feature table.
- `data/functional/picrust2/{KO,EC}_predicted.tsv`,
  `data/functional/picrust2/ko_pred_metagenome_unstrat.tsv` — large
  PICRUSt2 KO/EC tables (189 MB – 4.4 GB). Regenerate by running
  `picrust2_pipeline.py` on the ASV FASTA.

Raw sequencing reads are deposited at the European Nucleotide Archive
(ENA) under umbrella project `PRJEB104209`; amplicon reads under
`PRJEB106069`.

## Cross-desert comparison (Fig 11)

The cross-desert generalisation (Atacama / Namib / McMurdo) is
reproducible via the staged pipeline checked in at
`cache/crossdesert/` (paths are cluster-local; pipeline scripts are in
the `scripts/` of this repository):

```
stage1_download.sh        # ENA/SRA fastq retrieval
stage2_process.sh         # cutadapt + VSEARCH OTU-97
stage2b_reprocess.sh      # per-desert primer handling
stage3_classify.sh        # SINTAX taxonomy (drive5 SILVA v123)
stage4_analyze.py         # Shannon + direct CSP1-2 detection
stage5_genus_compare.py   # desert vs Empty-Quarter genus spectra
```

Inputs for the VSEARCH-based CSP1-2 detection are the four CSP1-2
MAG-associated ASVs recovered from Empty-Quarter co-assemblies
(`cache/csp1-2_asvs.fasta`, generated by `08_csp_mag.qmd`).
Raw fastqs (1.8 GB) are archived on Unimatrix at
`/data/emptyquarter/ecology-paper-runs/crossdesert/raw/`.

## Climate-impact and projection analyses (Figs 5–7)

Three complementary climate-impact analyses extend the digital twin
into present and future climate:

| Stage | Script | Output | Manuscript |
| --- | --- | --- | --- |
| Pull NASA POWER 1995–2024 daily T and P at 60 sites | `scripts/fetch_power_historical.py` | `cache/climate_historical_1995_2024.parquet` | Fig 7 |
| Per-site historical climate trends (Mann-Kendall + sign test) | `scripts/run_climate_trends.py` | `cache/climate_trends_per_site.tsv`, `figures/fig_climate_trends.pdf` | Fig 7 |
| CMIP6 SSP × horizon scenarios on the Tier-3 hierarchical twin | `scripts/run_cmip6_projections.py` | `cache/cmip6_interventions.tsv`, `figures/fig_cmip6_projections.pdf` | Fig 5 |
| CSP1-2 climate-niche projection on global 10-arcmin grid | `scripts/run_csp_niche_projection.py` | `cache/niche_model_coeffs.tsv`, `cache/niche_grid_summary.tsv`, `figures/fig_niche_projection.pdf` | Fig 6 |

**Reproduction order** (each script is self-contained and depends only
on the listed upstream artefacts):

```bash
# Per-site historical trends (Fig 7) — needs ~5 min of NASA POWER pulls
uv run python scripts/fetch_power_historical.py
uv run python scripts/run_climate_trends.py

# CMIP6 projections via the Tier-3 twin (Fig 5)
# Requires cache/causal_frame_tier1.parquet from notebook 11
uv run python scripts/run_cmip6_projections.py

# Global niche projection (Fig 6) — needs WorldClim historical + CMIP6 future rasters
# (~80 MB; download commands in .gitignore comment under data/worldclim/)
uv run python scripts/run_csp_niche_projection.py
```

`scripts/fetch_openmeteo_historical.py` is an alternative ERA5-archive
historical fetcher for higher-resolution data; in practice the
Open-Meteo free-tier rate limits made the NASA POWER fetcher more
reliable for the 60-site batch, so it is the one used for Fig 7.

CMIP6 ensemble Δ(T, P) values for the Arabian Peninsula are taken
from Almazroui et al. 2020 (Earth Syst. Environ. 4:611–630).
Climate-niche projections use WorldClim 2.1 historical 10-arcmin and
WorldClim 2.1 CMIP6 future bioclim grids (UKESM1-0-LL, 2081–2100,
SSP2-4.5 and SSP3-7.0); these are gitignored under `data/worldclim/`
(re-fetch commands documented at the top of `.gitignore`).

## Quickstart

```bash
# one-time setup
make env

# regenerate the cached feature-table + QC (runs 00_load_and_qc)
make cache

# render every figure notebook
make figures

# sanity-check outputs
make check
```

## Citing

If you use this code or its outputs, please cite the manuscript:

```bibtex
@article{tawfiq2026emptyquarter,
  title   = {A longitudinal, causally-grounded microbiome atlas of
             Earth's largest sand desert},
  author  = {Tawfiq, Rund and Martinez Arbas, Susana and Abdelhakim,
             Marwa and Niu, Kexin and Lopez Velazquez, Alejandra and
             Rueping, Magnus and Hoehndorf, Robert},
  year    = {2026},
  note    = {in preparation}
}
```

and the companion data paper describing the knowledge graph:

```bibtex
@article{rubalkhali2026kg,
  title   = {A semantic knowledge graph of the Rub al-Khali metagenomic
             expedition},
  author  = {{Bio-Ontology Research Group}},
  year    = {2026},
  note    = {companion data paper, in preparation}
}
```

## License

Code: MIT. Derived data tables retain their original source licences.
Raw reads: ENA accession terms.

## Contact

Bio-Ontology Research Group, King Abdullah University of Science and
Technology (KAUST). Issues and pull requests welcome.
