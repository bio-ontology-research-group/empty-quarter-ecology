# Control-library re-denoising and per-trip contaminant screen (2026-08-30)

Decision (Robert Hoehndorf, 30 Aug 2026): denoise the seven Trips 1-3 control
libraries that Marwa Abdelhakim listed on 30 Aug together with the six Trip 4
extraction blanks, update the contaminant analysis, and report the consolidated
result in both papers. Supersedes `../rerun-trip4-blanks-2026-08-29/` for the
blank table (same settings; the DADA2 error model is now learned from 14
libraries instead of 6). The 29 Aug directory keeps its own outputs untouched.

## What was denoised

Run `novaseq_14_07_25` (250626_A01018_0375_AH2JLYDRX7), 14 libraries, nf-core/ampliseq
2.14.0 with the canonical settings (`ibex/params.json`, `ibex/controls.config`,
`ibex/run_controls.sbatch`; Ibex job 51027761, run directory
`/ibex/user/hohndor/eq-controls-2026-08-30/`):

- Trip 4 extraction blanks `TRIP4_SEB`, `TRIP4_EB1`..`TRIP4_EB5`
  (M-25-0684, M-25-0770..0774);
- Trips 1-3 controls `CTL_ExtractionCtrlPro_Trip1` (M-25-0929), `CTL_PCRCtrl_Trip1`
  (M-25-0555), `CTL_Ctrl1_Trip1` (M-25-0323), `CTL_Ctrl2` (M-25-0553), `CTL_Ctrl3`
  (M-25-0554), `CTL_NegCtrl1_Trip2` (M-25-0870), `CTL_NegCtrl2_Trip2` (M-25-0871),
  and `CTL_NTC2` (M-25-0875, on the same run, not in Marwa's list).

Correction found afterwards: the eight Trips 1-3 libraries had already been
denoised in the July 2025 run (`ibex_20250714_qiime2`, e-prefixed ids
`e0323_Ctrl_1_Trip1` etc.) and carry recorded roles in
`data-paper/evidence/controls/source_snapshots/ibex_20250714_qiime2/controls-metadata.tsv`
and `analysis/v3/control_audit/control_analysis_roles.tsv`; only the canonical merged
table excludes them (the Trips 1-4 export was taken after control removal). The
Trip 4 blanks were never in any denoising run before 29 Aug.

## Files

- `ibex/`: samplesheet, params, config, sbatch, `job.id`, `wait_and_fetch.sh`, and
  `results/` (ASV table and md5-id table `controls_md5.tsv`, sequences, SILVA 138.2
  taxonomy, DADA2 stats, args, pipeline_info without the HTML reports).
- `characterize_controls.py` -> `outputs/control_characterization.{tsv,md,json}`:
  per library reads, ASVs, top genera, ZymoBIOMICS mock-genus share, and read share
  in ASVs that also occur in the 1,247 canonical biological profiles.
- `screen/`: the per-trip contaminant screen of `../rerun-trip4-blanks-2026-08-29/run_rerun.sh`
  run with `TRIP4_BLANK_TABLE=ibex/results/controls_md5.tsv` (outputs and logs copied
  here; the 29 Aug outputs were restored from git afterwards).
- `sensitivity/`: patched copies of `scripts/controls/build_control_sensitivity_inputs.py`
  (`--expected-profiles`, per-screen candidate check, Trip 4 profile ids) and
  `run_control_ecology_sensitivity.sh` (inputs from `screen/outputs`), plus
  `outputs/` with the 25-conclusion sensitivity on the Trip 5 + Trip 4 filtered table.

## Results

DADA2 (non-chimeric reads / ASVs): Trip 4 blanks 15,021-180,892 reads, 150-465 ASVs;
`Ctrl1_Trip1` 281,716 / 24; `Ctrl2` 322,257 / 2,816; `Ctrl3` 201,379 / 2,811;
`ExtractionCtrlPro_Trip1` 408,918 / 1,165; `NegCtrl1_Trip2` 105,840 / 2,253;
`NegCtrl2_Trip2` 87,301 / 2,268; `NTC2` 76,084 / 460; `PCRCtrl_Trip1` 28,427 / 1,686.

Characterisation (`outputs/control_characterization.md`): `Ctrl1_Trip1` is 99.99 %
ZymoBIOMICS genera (recorded as a negative control: label under review); `Ctrl2`,
`Ctrl3` (recorded positive) and `ExtractionCtrlPro_Trip1` (recorded extraction blank)
are soil-like (63-84 % of reads in ASVs shared with biological profiles); `NTC2` is
98 % *Aeromonas*, the reagent signature that also dominates the Trip 4 blanks
(23-34 % *Aeromonas* in SEB, EB3-EB5).

Per-trip screen (`screen/outputs/comparison.md`): Trip 5 screen reproduces the
published 351 candidates exactly; Trip 4 screen (6 blanks vs 95 linked profiles):
7 candidate ASVs (4 unassigned at genus level, *Nesterenkonia*, *Brevibacterium*,
*Alcaligenes*), 31 of 95 profiles lose at least one read, median removal 0.00 %
(IQR 0.00-0.01 %), maximum 14.31 %, pooled 0.14 %, no profile below 25,000 reads.
Identical to the 29 Aug six-library run.

Sensitivity of the 25 tracked conclusions with both screens applied together (312
filtered profiles; `sensitivity/outputs/headline_result_sensitivity.tsv`): all 25
verdicts stable; largest shifts paired Shannon q 0.0256 -> 0.0312 (published Trip 5-only
run: 0.0305) and paired distance-decay omnibus p 0.0053 -> 0.0075 (0.0077).
