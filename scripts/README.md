# `scripts/` — external data fetchers and cluster pipelines

These scripts pull remote sensing and climate data, submit
compute-heavy jobs to Slurm clusters, and build individual figure
panels. They are called either from notebooks or from the Makefile
`cache` and `figures` targets.

## Data fetchers

All fetchers are idempotent — if the target cache already exists
they exit quickly. Re-run with `--force` to refetch. Each writes
to `cache/<name>.parquet` (or equivalent TSV/netCDF).

### `fetch_nasa_power.py`
Pulls daily NASA POWER records (temperature, precipitation, humidity,
UV index, shortwave radiation, topsoil moisture) per site. No
credentials required. Writes
`data/climate/daily_weather.tsv` + `cache/nasa_power.parquet`.

**Usage:** `./fetch_nasa_power.py --sites data/geodata/trip5_geodata.tsv`

### `fetch_openmeteo_extended.py`
Legacy Open-Meteo fetcher (kept for cross-validation). Subject to
rate-limit 429s; prefer NASA POWER.

### `fetch_appeears_ndvi.py`
MODIS NDVI via NASA AppEEARS. Submits a point-sample request for
all sites. Requires Earthdata credentials in `.secrets/earthdata.netrc`
(gitignored). Writes `cache/modis_ndvi.parquet`.

**Usage:** `./fetch_appeears_ndvi.py --submit` (submits task) then
`./fetch_appeears_poll.py --task-id <id>` (polls until complete).

### `fetch_appeears_poll.py`
Resilient poller for long-running AppEEARS tasks with retry backoff
(60 s timeout, exponential backoff on 429). Downloads the CSV bundle
once the task completes.

### `fetch_smap_moisture.py`
SMAP Level-3 9 km daily soil moisture via the NASA Earthdata
`earthaccess` library. Requires Earthdata credentials. Writes
`cache/smap_moisture.parquet`. Used as a null-covariate test in
`15_intervention_scenarios.qmd`.

### `fetch_srtm_topography.py`
Copernicus GLO-30 DEM tiles from AWS Open Data (no credentials
required). Computes elevation, slope, aspect, and curvature per site;
writes `cache/srtm_topo.parquet`. Consumed by
`15_intervention_scenarios.qmd`.

## Figure builders

### `build_fig11_crossdesert.py`
Builds the four-panel Fig. 11 (Shannon distributions per desert,
CSP1-2 prevalence bars, Namib abundance–Shannon scatter, Atacama
EC–Shannon scatter). Reads `cache/crossdesert/per_sample.tsv`,
`comparison_summary.tsv`, `eq_shannon_reference.tsv`. Writes
`figures/fig11_crossdesert.pdf` and `.png`. Called directly by
`notebooks/16_cross_desert.qmd`.

## Analysis runners

### `run_bnti.py`
Command-line wrapper around `src/eq/bnti.py`. Computes β-NTI for one
compartment in parallel; used by `02_assembly.qmd`.

**Usage:** `./run_bnti.py --compartment surface --n-null 999 --threads 8`

## Slurm submission templates

### `scripts/slurm/`
Job scripts used on the Unimatrix and IBEX clusters:

- `ampliseq.sh` — nf-core/Ampliseq primary processing.
- `gspa_run.sh` — GSPA-DM / gapseq dark-matter annotation for the
  four CSP1-2 MAGs. Uses a wrapper script written into `/tmp` at
  job start to work around `noexec` on `/data` compute nodes.
- `picrust2.sh` — PICRUSt2 pipeline.
- `gtdbtk.sh` — GTDB-Tk classification for MAGs.

## Cross-desert pipeline (cluster-side, referenced by notebook 16)

These five scripts live on the Unimatrix cluster under
`/data/emptyquarter/ecology-paper-runs/crossdesert/` and produce the
`cache/crossdesert/*.tsv` files consumed by
`notebooks/16_cross_desert.qmd`. They are shipped here
(under `scripts/crossdesert/`) for auditability:

- `stage1_download.sh` — ENA/SRA fastq retrieval via filereport API.
- `stage2_process.sh` / `stage2b_reprocess.sh` — cutadapt primer
  trimming + VSEARCH dereplication + OTU-97 clustering. Atacama
  primer handling differs because Qiita-pipeline reads are
  pre-trimmed.
- `stage3_classify.sh` — SINTAX taxonomy against drive5 SILVA v123
  (does not include CSP1-2; see note below).
- `stage4_analyze.py` — per-desert Shannon, OTU richness, and CSP1-2
  VSEARCH `-usearch_global` at 97% / 85% V4 identity against
  `cache/csp1-2_asvs.fasta` produced by `08_csp_mag.qmd`.
- `stage5_genus_compare.py` — genus-level median relative abundance
  and desert-vs-EQ Spearman comparison.

**Note on SINTAX:** drive5 SILVA v123 lacks Dadabacteria taxonomy,
so CSP1-2 detection uses direct VSEARCH alignment (not SINTAX).
Replacing the reference with SILVA 138.2 would improve taxonomy
assignment cosmetically but is not required for the quantitative
results.

## Credentials

Scripts that need external credentials read them from `.secrets/`
(gitignored):

- `.secrets/earthdata.netrc` — Earthdata login for AppEEARS + SMAP
  (machine-readable `netrc` format).
- `.secrets/opentopography_api.key` — OpenTopography API key (used
  by an older SRTM fetcher; current Copernicus DEM path does not
  require it).

Obtain Earthdata credentials from
`https://urs.earthdata.nasa.gov/users/new`.
