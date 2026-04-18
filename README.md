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
