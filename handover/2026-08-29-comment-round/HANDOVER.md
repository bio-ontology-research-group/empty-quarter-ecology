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

## 8. Addendum, 30 Aug 2026 (second assistant session, ~/Documents clone)
- Supersedes parts of section 7: the Chrome extension was connected in this
  session; the KG paper's Overleaf Review panel holds only Michel's two
  comments on 01_introduction.tex (applied 29 Aug, anchors confirmed). The
  ecology project's panel was not opened here.
- Ecology CI: the cache-path fix alone still fails at render_figures.py (the
  renderer refuses any runtime but Python 3.11.14 / matplotlib 3.9.4 /
  FreeType 2.14.3); validate.yml now checks out empty-quarter-ecology-data and
  builds that runtime with setup-micromamba from its conda lock. Committed on
  top of 501f802.
- Data repo: the ~/Documents/papers/empty-quarter-data-paper clone carries the
  manifest rebuild plus restored paper/README.md, PRE_SUBMISSION_CHECKLIST.md,
  rebuilt PDFs, regenerated bibliography custody and environmental audit, and
  ontology-identifier corrections to the manuscript listings (KG ledger,
  "Applied 30 Aug 2026"). The unpushed 5a18615 in
  ~/Public/software/empty-quarter-paper-repos/data-paper conflicts on
  FILE_MANIFEST.tsv and should be dropped or rebased.
- Public SPARQL endpoint serves older dna/xrf modules: all Section
  "Technical validation" queries return 0 rows there; redeploy before submission.

## 8. Addendum, 30 Aug 2026 evening (assistant session; Robert's decisions applied)
- Overleaf review panels cleared (Chrome connected): ecology 70 tracked changes accepted
  (all Susana's, reviewed, none needing new action); KG paper: Michel's 2 comments resolved,
  3 changes accepted. Abderahmane's 30 Aug Overleaf pass (rain -> rainfall, commas) kept;
  his slips fixed (\label{fig: landscape}, west-east, landscape-scale).
- Name policy applied to the ecology paper: Rub' al-Khali introduced once (Latin + Arabic
  via arabtex, pdflatex; LaTeX \begin/\end restored after loading), Empty Quarter thereafter.
  Data paper not yet harmonised.
- Replies applied: Daniela (SP flow cell), Abderahmane and Kexin (approval + ORCID),
  Marwa (controls; corrected workbook; S40PRr1 failed amplification, S34Dr1 dropped
  without mention). Both papers pushed to Overleaf and GitHub after every change.
- Controls: 14 control libraries of run novaseq_14_07_25 denoised on Ibex (job 51027761),
  analysis/rerun-controls-2026-08-30. Trip 4 screen 7 candidates / 0.14 % pooled; all 25
  conclusions stable with both screens. Reported in both papers. Finding: the Trips 1-3
  controls were already in the July 2025 run with recorded roles; Ctrl_1_Trip1 is a Zymo
  standard despite its negative label, Ctrl_2/Ctrl_3/Extraction_Ctrl_Pro are soil-like.
- Mail: Robert sent the Marwa reply and the Rund mail (13:53). A corrected follow-up to
  Rund (label discrepancies) is drafted in Gnus (*claude-mail-31*) for Robert to send.
- Open: Rund's answers (extraction days Trips 1-3, plate links, the four mislabelled
  libraries); grant numbers; shotgun BioProject; PRJEB accession; Sulaiman's ORCID;
  ten co-authors still silent (Maxat, Mohammed, Hind, Michel, Raik, Krishnakumar,
  Tiannyu, Magnus); data-paper consistency test (14 failures); public SPARQL endpoint
  serves old modules (KG ledger).

## 9. Addendum, 31 Aug 2026 (assistant session; Rund's reply)

Rund replied 31 Aug 11:47 (from ECCB; realistic target for her items: end of week,
not Wednesday). Robert's reply drafted in Gnus (*claude-mail-5*): accepts Friday,
pastes Marwa's Trip 4/5 answers (Rund cannot open Mattermost: previews only, no
open/reply; Marwa's 25 Aug questions to her unread), reformulates Q1 as kit
linkage, promises the index-728/run lookup.

### Facts from Rund that fold into the papers
- Trips 1-3 extraction blanks were per extraction KIT, not per extraction day:
  one blank for PowerSoil and one for PowerSoil Pro. Presumably
  Extraction-Ctrl-Pro-Trip1 (M-25-0929) is the Pro blank; whether a PowerSoil
  blank was sequenced, and its ID, asked in Robert's reply. Supplement S2
  ("Assay-aware control analysis", ~line 197-210) currently says one blank per
  extraction day (true for Trip 5 workbook) and hedges Trips 1-3; once Q3 is
  answered, S2 should state kit-level blanks for Trips 1-3 explicitly. Main text
  only claims day linkage for Trip 5 (main.tex ~748), so no change there.
- Trips 1 and 2 were extracted on the same days (cross-trip extraction batches);
  confirms the S2 hedge "may have combined samples from different trips".
- Nothing was sequenced on extraction days for Trips 1-3; no per-day
  extraction blank exists or can be reconstructed. Do not wait for notebook
  dates (Rund's notebooks not physically with her; Robert told her not to dig).
- Trip 3 PCR blank: exactly one prepared, index "xGen 10nt UDI Index Pair 728".
  Whether sequenced/library ID unknown; Robert owes the lookup (map index 728 in
  run sample sheets; determine which run(s) carry Trip 3 samples). If Trip 3
  samples sit on another run, Rund can find that run's blank index too.
- Q3 (roles/contents of the 8 July-2025 control libraries, esp. the four
  suspected mislabels Ctrl-1-Trip1, Ctrl-2, Ctrl-3, Extraction-Ctrl-Pro-Trip1):
  Rund still checking. Marwa's catalogue answer (D6322 HMW standard used for
  Trips 1-2) matches Ctrl-1-Trip1's Zymo profile, supporting a label swap.

### Manuscript coordination
- Rund is still editing the ecology manuscript on Overleaf and reconciling
  inconsistencies between the two papers (sample/compartment descriptions,
  methods wording). Her commit 651e276 (31 Aug) landed in Methods.
- Cleanup needed in 651e276 (fix after she finishes or coordinate): empty
  \ref{} main.tex:824; "campagin" typo twice (826-827); Unicode en-dash in
  "Benjamini–Hochberg" (827); "Each sampling campaign covered 60 sites"
  (~674) now sits oddly against the 4 Trip-1-only locations paragraph below.
- She rewrote the environmental-measurements methods paragraph (site+campaign
  averaging, campaign + linear/quadratic transect covariates, site-clustered
  SEs, BH over nine tests): verify this matches the analysis actually run
  before accepting.
- Comms with Rund by email only until her Mattermost access is fixed; the
  Marwa<->Rund Mattermost loop is broken in both directions.
