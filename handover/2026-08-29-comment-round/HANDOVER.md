# Handover: Empty Quarter papers, co-author comment round (29-30 Aug 2026)

Written by Robert's assistant session of 29 Aug 2026. Robert takes it from here.
Ledgers with every comment and its status: ecology paper
`~/Documents/papers/empty-quarter-ecology-overleaf/comments-2026-08-29.md`;
KG/data paper `~/Documents/papers/empty-quarter-data-paper/comments-2026-08-29.md`.

## 1. What was done and where it is

### Ecology paper ("Landscape-scale bacterial biogeography across the Rub' al-Khali ...")
- Source of truth: Overleaf project 6a7422eef2624277cf100b3b (the older project
  686e50288b374b13ddfeef7e is retired). Local clone:
  `~/Documents/papers/empty-quarter-ecology-overleaf` (commits e3f5578, 4571938,
  abfbb40, 3ff4068, all pushed to Overleaf main).
- Applied: ORCIDs (comment block after the author list; class has no ORCID field);
  Funding section with Susana's FNR statement (KAUST grant numbers still TODO);
  Daniela's sequencing wording; Marwa's amplification split (TaKaRa 30-cycle first
  round, 1 % gel, QIAGEN 8-cycle indexing), "24 control profiles (23 extraction
  blanks and 1 PCR blank)" (now known to be wrong, see section 2), Trip 4
  limitation wording, new Methods subsection "Relic-DNA removal and viability
  assessment"; "soil position" -> "compartment" everywhere, root-adjacent sentence
  restored in Methods with the one-time explanation (rhizosphere targeted, plants
  desiccated, soil as close to the root as possible, sometimes root-attached),
  short term "root-adjacent soil" elsewhere; Susana's Eida et al. passage replaced
  by "Our near-root samples lie in a compartment Eida et al. did not sample; their
  lower evenness extends this local organization outward from the root without
  identifying a host filter"; intro results-summary paragraph kept; Figure 1 c/d
  swapped, Fig 1b legend/ticks, Fig 3 axis titles (regenerated, repro commit
  df29c75, pushed); Table S1 row "Quality-controlled profiles per expedition"
  325/24/478/177/233; Methods sentence on the 25,000-read depth (no rarefaction
  analysis exists); "Rub' al-Khali" throughout, "Empty Quarter" once (kept in the
  keyword list, decide); shotgun source stated (125 Trip 5 libraries, separate
  study) with a red \todo{BioProject accession}; Software subsection naming Claude
  with author verification; Susana's rain-section rewrite adopted (all numbers
  agree with S7, p = 0.0485 unchanged); typo fixes from the co-author pass.
- Not applied (Robert): flow-cell wording (S1 vs SP, asked Daniela + Marwa),
  Marwa's control numbers (see 2), grant numbers, "Empty Quarter" keyword.

### KG / data paper ("The Empty Quarter database and knowledge graph")
- Source of truth: Overleaf project 6a741d888df17653bc201f8a, pushed to 0cb2b9c.
  GitHub `bio-ontology-research-group/empty-quarter-data-paper` main 284577f: the
  diverged histories were merged; `paper/` equals the Overleaf tree; both local
  clones (`~/Documents/papers/empty-quarter-data-paper`, `~/Public/software/
  empty-quarter-paper-repos/data-paper`) are at 284577f.
- Applied: ORCIDs, FNR funding; Marwa's "Laboratory controls", "DNA extraction",
  "16S amplification" text as written (her "waiting for Rund's input" is a \todo);
  Daniela's sequencing sentences; Michel: digital-twin paragraph moved to future
  work, consistency checking stated as OWL 2 EL with ELK, abstract now "expressed
  in OWL 2 DL ... checked under the OWL 2 EL profile with ELK"; Rund: full-IRI
  prefix example, refersTo -> isMeasurementValueOf everywhere (verified against
  rdf/generators), typos; ChEBI 247 (3 Dec 2025) and PubChem snapshot 23 Jul 2026
  pinned; validation tool versions moved to Methods, 05_validation results-only;
  XRF reporting rule: Trips 1-4 workbooks repeat identical values under several
  sections so max = last; Trip 5 has real duplicate readings for one sample
  (V14Dr1), max vs last differ in 3 cells (Al 1.1/1.0, Na 0.35/0.30, Si 21.5/20.2),
  one rule would change 1 of 725 lab records; Methods say this in one paragraph;
  rank agreement reported as spearman_both_positive (12 analytes, median 0.21,
  range -0.40 to 0.48); new Methods subsection "Shotgun-derived companion inputs"
  (150 CoverM tables, eggNOG-mapper 2.1.12, 990-genome KO matrix; assembly/binning
  versions and shotgun accessions are not recorded anywhere, \todo).
- Not applied (Robert): Data availability accession (Rund now states PRJEB104209
  and dropped the PRJEB106069 caveat; 04_data_records still names PRJEB106069;
  verify before submission); Marwa's control numbers; Acknowledgements and CRediT.

### Repositories
- `empty-quarter-ecology` README: "private" -> public, title aligned; GitHub
  description set to the paper title. The data repo is the same GitHub repository
  renamed to `empty-quarter-ecology-data`.

## 2. Findings that change the manuscript text (controls)
1. The six Trip 4 extraction blanks (M-25-0684_SEB, M-25-0770_EB1 .. 0774_EB5,
   run novaseq_14_07_25) were sequenced but never denoised: absent from the July
   2025 samplesheet and the canonical ASV table. The EB1-EB5 columns in the table
   are Trip 5 blanks (label collision). All 24 control columns (EB1-EB18,
   Negative1/2/4-7) are Trip 5 libraries. So the sentences "6 other extraction
   blanks from Trip 4 lacked a day mapping and were characterized" and "24 control
   profiles (23 extraction blanks and 1 PCR blank)" are wrong and must be
   rewritten once Marwa answers (which Negative libraries are extraction blanks vs
   PCR-stage NTCs).
2. The Trip 4 PCR blank (index pair 275) has no library on Ibex; it was never
   sequenced.
3. The six Trip 4 blanks are now denoised (Ibex job 50996981, ampliseq 2.14.0,
   same primers/DB/settings as the canonical runs, truncLen 240/238, md5 ids):
   1,824 ASVs. Per-trip contaminant screen re-run with them:
   `analysis/rerun-trip4-blanks-2026-08-29/outputs/comparison.md`:
   Trip 5 screen reproduces the published 351 candidates exactly; Trip 4 screen:
   7 candidate ASVs (genera: 4 unassigned, Nesterenkonia, Brevibacterium,
   Alcaligenes), 95 linked profiles, 31 lose at least one read, median removal
   0.00 % (IQR 0.00-0.01 %), max 14.31 %, pooled 0.14 %, no profile below 25,000
   reads. The downstream 25-conclusion sensitivity stage was NOT rerun (it
   hard-codes the 217 Trip 5 profiles; needs a patched copy). Batch map: 95 of
   97 Trip 4 IDs matched; S34Dr1 and S40PRr1 (EB4) not in the table; 82 profiles
   / 28 sites without a sequenced blank (matches Marwa's 28 of 60).
4. Marwa's workbook inconsistencies (asked her): sheet EB declares 15 samples,
   lists 13 (site 59 missing); EB4 declares 26, lists 23; EB2 says 19 in the sites
   column but IDs are S14.
5. Marwa's numbers vs the audit (asked her): 9 batches / 28 of 60 sites vs our 6
   blanks / "23 sites"; extraction days trip-specific vs mixed; field negatives
   on Trips 1-3 and 5 vs no sterile-bag inventory in the records.
6. Primer discrepancy (Robert to decide): the data paper's
   `metadata/amplicon/README.md` and Methods say reverse primer Bakt_785R
   (GACTACHVGGGTATCTAATCC), but both canonical Ibex runs recorded 806RB
   (GGACTACNVGGGTWTCTAAT) in params.json; the new run reused 806RB.

## 3. Co-author status (29 Aug)
Responded with approval/ORCID: Daniela, Susana (FNR statement), Jood, Xiang,
Alejandra; Sulaiman approves (ORCID and explicit preprint line still to collect);
Rund and Marwa responded with edits/text. Reminder sent 29 Aug (reply asked by
4 Sep) to Maxat (coolmaksat@gmail.com; KAUST address bounces), Mohammed Alarawi,
Hind Aldakhil, Abderahmane Derouiche, Michel Dumontier (edited the KG paper, no
mail), Raik Gruenberg, Kexin Niu, Krishnakumar Sivakumar, Tiannyu Wang, Magnus
Rueping. Susana's affiliation is now "NIUM, Esch-sur-Alzette" (check the full
name).

## 4. Emails (texts in emails/)
- 29 Aug, to Daniela + Marwa, "NovaSeq flow cell" (S1 v1.5 vs SP), SENT.
- 29 Aug, to Marwa, "Empty Quarter controls: a few points to reconcile", with
  control_analysis_roles.tsv, trip5_extraction_batch_summary.tsv, summary.json
  attached, SENT.
- 29 Aug, to Marwa, follow-up (Trip 4 blanks never denoised, PCR blank 275,
  workbook checks, Negative libraries), SENT.
- 29 Aug, to the ten silent co-authors, approval + ORCID reminder, SENT.
- 29 Aug, to Sulaiman, thanks, final version to follow, SENT.
- Answers awaited from Daniela/Marwa (flow cell), Marwa (controls, two mails),
  the ten co-authors.

## 5. Open decisions for Robert
1. Controls text in both papers after Marwa's answers (section 2.1, 2.2, 2.5).
2. Whether to rerun the downstream 25-conclusion sensitivity with the Trip 4
   screen (needs a patched copy of scripts/controls/build_control_sensitivity_
   inputs.py + run_control_ecology_sensitivity.sh).
3. Primer wording 785R vs 806RB in the data paper.
4. Grant numbers in both Funding sections; KG Acknowledgements and CRediT.
5. Shotgun BioProject accession and assembly/binning/CoverM versions (nobody
   recorded them; likely Krishnakumar or Tiannyu).
6. Data availability accession (PRJEB104209 vs PRJEB106069).
7. "Empty Quarter" in the ecology keyword list.
8. Overleaf comment threads (not the tracked changes) were never read: the
   Chrome extension was not connected. Open the Review panels of both projects.
9. Workflow scheme for the repo README / Supplementary Figure S1 (Susana).

## 6. Ibex artefacts
Run dir `/ibex/user/hohndor/eq-trip4-blanks-2026-08-29/` (symlinked as
`.../novaseq_14_07_25/raw_reads/ampliseq/results_trip4_blanks_2026-08-29`;
/ibex/project/c2014 is 100 % full). Local copies: `analysis/rerun-trip4-blanks-
2026-08-29/ibex/results/` (md5 table, ASV table/seqs, DADA2 stats, args,
pipeline_info, slurm log). Re-run the screen any time with
`TRIP4_BLANK_TABLE=$PWD/ibex/results/trip4_blanks_md5.tsv bash run_rerun.sh`
(MODE=pooled for a single 23-blank screen).

## 7. Addendum, 30 Aug 2026 (assistant session)
- Mail checked 30 Aug morning: no replies yet from Daniela/Marwa (flow cell,
  controls) or the ten reminded co-authors; only Daniela's auto-reply.
- Chrome extension still not connected; Overleaf comment threads (item 5.8)
  remain unread.
- GitHub `validate` CI was red on both repos (since 5 Aug, not caused by the
  comment round). Ecology repo: setup-python `cache: pip` found no
  requirements.txt; fixed with `cache-dependency-path` (commit 501f802). The
  Trip 4 blank re-run had tracked three nextflow pipeline_info HTML reports,
  which verify_repository.py rejects (blanket .html rule); untracked and
  gitignored, files kept on disk. Local verify_repository.py passes.
  Data repo: FILE_MANIFEST.tsv was stale after paper/ followed the Overleaf
  tree (8 files added, 4 removed); rebuilt (commit 5a18615 in
  ~/Public/software/empty-quarter-paper-repos/data-paper; the
  ~/Documents/papers/empty-quarter-data-paper clone does not have it yet).
  Neither commit is pushed.
- Still red in the data repo: scripts/manuscript/test_manuscript_consistency.py
  fails 17 tests + 1 error against the Overleaf-authoritative paper/ tree. The
  tests expect statements the current manuscript no longer carries (nine
  paired PMA aliquots, the F46Dr2 reconciliation, metadata/DATA_DICTIONARY.tsv
  citation, pressure "multiplied by 100 and asserted in pascals", Field-XRF vs
  lab-XRF separation, "two separately configured acquisitions" for climate,
  value-to-quality relation directions, Zenodo staging row wording), the
  removed paper/PRE_SUBMISSION_CHECKLIST.md, a stale env_table.tex, and
  "Lopez Velazquez" where the co-author now writes "Lopez-Velazquez". Robert
  to decide per test whether the manuscript lost reconciled content during the
  Overleaf takeover (port it back) or the test is stale (update the test).
