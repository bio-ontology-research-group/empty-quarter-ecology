# Re-run of the extraction-blank contaminant screen with the Trip 4 blanks

Date: 2026-08-29. Requested by Robert Hoehndorf after Marwa Abdelhakim's mail of
26 Aug 2026 (Message-ID <D5204C75-7A9B-470A-A883-E0CA75DB1749@KAUST.EDU.SA>,
`~/Mail/inbox/360188`), which attached the Trip 4 EB-sample map and asked whether
the Trip 4 extraction blanks were used in the contamination audit.

Nothing outside this directory was modified. No commit was made.

## What the published screen is

`scripts/controls/run_assay_aware_control_audit.py` (results in
`analysis/v3/control_audit/`, manuscript Methods "Assay controls and sensitivity
analyses"). It is a prevalence screen in the style of decontam's prevalence
method, implemented directly in Python on the canonical ASV table
(`data/processed/taxonomy/taxon-tables/feature-table-trips1-5.tsv`, 351,472 ASVs
x 1,271 profiles):

- training blanks: the 17 Trip 5 extraction blanks EB1-EB17 (columns of the
  canonical table), linked by extraction day to 220 Trip 5 profiles via
  `data/metadata/samples/Sequenced_Samples_by_EB_FifthTrip.xlsx` (217 of the
  220 are in the canonical table; V16Dr1, V16Sr1, V7Dr1 are not);
- per ASV: 2x2 table of presence in the 17 blanks vs presence in the 217 linked
  profiles, one-sided Fisher exact test (enrichment in blanks);
- candidate when p < 0.10 AND present in >= 2 blanks AND blank prevalence >
  biological prevalence (a 3x3 grid of thresholds 0.01/0.05/0.10 x 1/2/3 blanks
  is also reported); the frequency method was not run (no DNA concentrations);
- the 351 candidates are removed only from the 217 linked profiles in a
  sensitivity copy (`trip5_mapped_feature_table_control_filtered.tsv.gz`);
  `build_control_sensitivity_inputs.py` + `run_control_ecology_sensitivity.sh`
  then re-run 25 tracked conclusions on that copy (all stable).

## The finding that blocks a straight re-run

The six Trip 4 extraction blanks in Marwa's map (EB, EB1-EB5; xGen UDI index
pairs 168, 216, 228, 240, 252, 264) were sequenced on run `novaseq_14_07_25` as
libraries `M-25-0684_SEB_SU0168`, `M-25-0770_EB1_SU0216`, `M-25-0771_EB2_SU0228`,
`M-25-0772_EB3_SU0240`, `M-25-0773_EB4_SU0252`, `M-25-0774_EB5_SU0264`
(`data/metadata/samplesheets/additional_fastqs_v2.tsv`; FastQC in the working
repo shows 23k / 246k / 233k / 42k / 36k / 32k read pairs). They were **never
denoised**: they are absent from the July 2025 denoising samplesheet
(`ibex_20250714_16s_samplesheet.tsv`, 1,046 libraries), from the July 2025 QIIME 2
snapshot (1,041 samples), and from the canonical Trips 1-5 table. The raw and
trimmed FASTQs exist only on Ibex (`/ibex/project/c2014/EmptyQuarter_Data/soil/
amplicon_16S/novaseq_14_07_25/{raw_reads,trimmed_reads}/`). No FASTQ is on this
laptop.

The canonical table's `EB1`-`EB5` columns are the **Trip 5** blanks (index pairs
307/410/319/331/422, `ibex_trip5_16s_samplesheet.tsv`), not the Trip 4
libraries of the same label. The 24 control columns of the canonical table are
EB1-EB18 and Negative1/2/4/5/6/7, all from the Trip 5 sequencing run
(`control_ground_truth.tsv` CGT-007/CGT-008). Consequently the manuscript
sentence "6 other extraction blanks, from Trip 4, lacked a day mapping in the
analysis inputs and were characterized" is not what the data show: the seven
characterization-only profiles (EB18, Negative1/2/4-7) are Trip 5 libraries with
no extraction-day mapping, and no Trip 4 blank profile exists in any analysis
input. The answer to Marwa's question is therefore: no, the Trip 4 blanks were
not used, and could not have been.

Including them requires denoising the six libraries (plus, if wanted, the Trip 4
PCR blank, index pair 275, which has no FASTQ record in the samplesheets) with
the same nf-core/ampliseq v2.14.0 / DADA2 settings as the canonical run so that
the md5 ASV identifiers are comparable, on Ibex where the reads are. That is a
cluster job, not a laptop job, and it was not submitted (per instructions). Once a
Trip 4 blank ASV table exists, the screen below runs in about 5 minutes here.

## Extended inputs built here (step 1)

`build_trip4_batch_map.py` -> `inputs/`:

- `EB_Sample_Map_FourthTrip2 correct.xlsx`, `marwa_2026-08-26_mail.txt`: Marwa's
  attachment and mail text (saved from the mail; originals untouched).
- `extraction_batch_map_extended.tsv`: 318 rows (blank, sample ID, canonical
  column), Trip 4 and Trip 5 together.
- `blank_libraries.tsv`: one row per blank with index pair, library name, FASTQ
  path, and status.
- `mapping_report.json`: mismatches (below).

Mapping check against the canonical table columns:

- Trip 4: 6 blanks, 97 sample IDs listed, **95 match** canonical columns
  (`e####_S<site><PR|S|D>r1`). Unmatched: `S34Dr1`, `S40PRr1` (both EB4; not in
  the canonical table at all). Index pairs in the map agree with the FASTQ
  names for all six libraries. Internal inconsistencies in the workbook: EB
  declares 15 sequenced but lists 13 IDs (site 59: only S59Sr1 listed, while
  S59PRr1/S59Dr1 exist in the table and site 59 is also in the "not sequenced"
  note); EB4 declares 26 but lists 23; EB2's "Sites Sampled" says 19 but the IDs
  are S14*, and site 19 sits under EB5. 82 Trip 4 profiles (28 sites) have no
  sequenced blank, matching Marwa's "28 of 60 sites".
- Trip 5: unchanged from the published run (220 listed, 217 matched, V16Dr1,
  V16Sr1, V7Dr1 absent).

## What was launched (step 2)

`bash run_rerun.sh` under nohup (log `run.log`, PID in `run.pid`; Python from
`../data-paper/.venv`, the documented lighter environment). Because no Trip 4
blank profiles exist, this run is the **validation** configuration: the
generalised screen `screen_extended.py` with the same parameters as the published
run (p < 0.10, >= 2 blanks, blank prevalence > biological prevalence), Trip 5
blanks only, Trip 4 screen reported as `not_run`. It verifies that the extended
harness reproduces the published 351 candidates and the per-profile removal
fractions, so that the full re-run is a one-line change once the Trip 4 blank
table exists:

    TRIP4_BLANK_TABLE=/path/to/trip4_blanks.tsv bash run_rerun.sh        # per-trip screens
    MODE=pooled TRIP4_BLANK_TABLE=/path/to/trip4_blanks.tsv bash run_rerun.sh

`trip4_blanks.tsv` is a `#OTU ID` x {SEB, EB1..EB5} count table (or a BIOM file)
with md5 ASV identifiers. `--mode per-trip` (default) screens the Trip 4 blanks
against the 95 linked Trip 4 profiles and the Trip 5 blanks against the 217
Trip 5 profiles, each candidate set applied to its own trip (extraction days are
trip-specific, so this is the batch-aware design and matches Marwa's
"characterized separately"). `--mode pooled` uses all 23 blanks vs 312 profiles.

Outputs (`outputs/`): `primary_contaminant_calls.tsv` (with a `screen` column),
`filter_sensitivity.tsv`, `removal_fraction_by_{profile,campaign,compartment}.tsv`,
`extraction_batch_summary.tsv`, `trip5_mapped_feature_table_control_filtered.tsv.gz`
(published format), `mapped_feature_table_control_filtered.tsv.gz` (Trip 5 +
Trip 4 linked profiles), `summary.json`, and `comparison.md` from `compare.py`.

The downstream 25-conclusion stage was not re-run: `build_control_sensitivity_inputs.py`
hard-codes 217 profiles and V-prefixed (Trip 5) IDs, so extending it to Trip 4
profiles needs a patched copy (parse `e####_S..` IDs, drop the 217 assertion).

## Comparison with the published run (step 3)

See `outputs/comparison.md` (regenerate with `../../../data-paper/.venv/bin/python compare.py`).

Validation run of 2026-08-29 (run.log; finished in about 2.5 min on the laptop):

- Trip 5 screen: 351 candidate ASVs, identical to the published set (0 new, 0
  dropped); 217 filtered profiles; removed reads median 0.40% (IQR 0.16-0.99%),
  max 56.60%, pooled 2.19%; 0 profiles below 25,000 reads; per-profile removed
  read counts identical to the published run for all 217 profiles; the filtered
  table body is byte-identical to the published one apart from the provenance
  header line (sha256 of body: rerun 3fc68a16c87913c8..., published 3fc68a16c87913c8...).
- Trip 4 screen: not run (no Trip 4 blank profiles). The 95 linked Trip 4
  profiles are listed in outputs/removal_fraction_by_profile.tsv with role
  mapped_biological_profile_screen_not_run.

Net effect on the manuscript today: none of the reported numbers change. What
changes is the description of the blanks (see "The finding that blocks a
straight re-run"), and whether Robert wants the six Trip 4 blank libraries
denoised on Ibex so the Trip 4 screen can actually be run.
