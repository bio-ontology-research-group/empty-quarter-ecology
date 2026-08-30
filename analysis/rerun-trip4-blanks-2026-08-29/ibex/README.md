# Ibex re-denoising of the six Trip 4 extraction blanks (2026-08-29)

Decision (Robert, 29 Aug 2026): denoise the six Trip 4 extraction-blank libraries
with the canonical nf-core/ampliseq v2.14.0 settings and run the per-trip
contaminant screen with them (`../run_rerun.sh`).

## What is on Ibex

- FASTQs (all six present, raw reads, 1.9 to 24 MB per file):
  `/ibex/project/c2014/EmptyQuarter_Data/soil/amplicon_16S/novaseq_14_07_25/raw_reads/`
  `M-25-0684_SEB_SU0168`, `M-25-0770_EB1_SU0216`, `M-25-0771_EB2_SU0228`,
  `M-25-0772_EB3_SU0240`, `M-25-0773_EB4_SU0252`, `M-25-0774_EB5_SU0264`
  (`_L002_R{1,2}_001.fastq.gz`).
- Trip 4 PCR blank (index pair 275): no library. The only "275" files in the run
  are Trip 3 xGen-275 libraries (F31PRr3, F10Sr3) and `M-25-0495_50Dr1_UDP0275`;
  nothing with `SU0275`. It was not sequenced under this run and is not in the
  job. (The run does contain `M-25-0875_NTC-2_SU0312` and
  `M-25-0555_PCR-Ctrl-Trip1_UDP0374`, neither of which is the Trip 4 PCR blank.)
- Canonical provenance (this is what the 351,472-ASV table actually is):
  - Trips 1-4: `novaseq_14_07_25/raw_reads/ampliseq/` (Aug 2025), params in
    `results/pipeline_info/params_2025-08-11_11-05-34.json`; DADA2 output
    1,117,343 ASVs x 1,041 samples (`parallelize/ASV_table.tsv`), then QIIME 2
    filtering to `final_analysis/tables/table-decontam-filtered-clean.qza`.
  - Trip 5: `Trip5/analysis/ampliseq/results_trip5_fixed/` (Mar 2026,
    `run_ampliseq_trip5_fixed.sh`), 331,305 ASVs x 258 samples.
  - Merge: `Trip5/analysis/ampliseq/merge_qiime2.sh` -> `merged_qiime2/feature-table.tsv`
    (351,472 ASVs; ids are QIIME 2 md5-of-sequence; verified that the canonical
    fasta header equals md5 of its sequence).

## Parameters reused (identical in both canonical runs)

`-profile kaust` (slurm executor, singularity, image cache
`/ibex/user/hohndor/.singularity/nf_images/`), `-r 2.14.0`,
`--FW_primer CCTACGGGNGGCWGCAG --RV_primer GGACTACNVGGGTWTCTAAT`,
`--illumina_novaseq`, `--dada_ref_taxonomy silva=138.2`, `max_ee 2`, `min_len 50`,
`sample_inference independent`, `seed 100`, `trunc_qmin 25`, `trunc_rmin 0.75`,
no ASV length/SSU filter. Both runs auto-derived `truncLen = c(240, 238)`
(`dada2/args/filterAndTrim.args.txt`); here it is passed explicitly
(`--trunclenf 240 --trunclenr 238`) so that six low-depth blanks do not
re-derive a different truncation. Deviations: `--skip_qiime` (QIIME 2 steps do
not touch `dada2/ASV_table.tsv`), smaller resource requests
(`trip4_blanks.config`), and of course the input and outdir. The DADA2 error
model is learned from the six blanks only (the canonical runs learned it from
all libraries of their run), which is unavoidable when denoising a subset.

Note for the data paper: `data-paper/metadata/amplicon/README.md` gives the
canonical command with `--RV_primer GACTACHVGGGTATCTAATCC` (Bakt_785R), but
both Ibex runs whose merge is the canonical table recorded
`GGACTACNVGGGTWTCTAAT` (806RB) in their params.json. The 806RB sequence was
reused here because the ASV sequences (and md5 ids) must match the table as it
was actually produced.

## Files

- `samplesheet_trip4_blanks.tsv`: ampliseq TSV (`sampleID forwardReads reverseReads run`),
  sample ids `TRIP4_SEB`, `TRIP4_EB1`..`TRIP4_EB5` (what `screen_extended.py` expects).
- `params.json`, `trip4_blanks.config`, `run_trip4_blanks.sbatch`: the job.
- `make_md5_table.py`: rewrites `dada2/ASV_table.tsv` with md5-of-sequence ids
  (`#OTU ID` header) -> `trip4_blanks_md5.tsv`; run at the end of the sbatch job.
- `job.id`: SLURM job id. `finish.sh`: rsync + `run_rerun.sh` (per-trip) + compare.

## Ibex layout

Run directory `/ibex/user/hohndor/eq-trip4-blanks-2026-08-29/` (`/ibex/project`
is at 100%); symlinked next to the previous run as
`novaseq_14_07_25/raw_reads/ampliseq/results_trip4_blanks_2026-08-29`. Canonical
outputs untouched. Log: `/ibex/user/hohndor/eq-trip4-blanks-2026-08-29/slurm-<jobid>.out`;
`DONE` is created after the md5 table is written.
