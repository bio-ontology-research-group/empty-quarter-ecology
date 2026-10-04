# Trip 4 control-removal diversity summary

This additive reporting check fills the supplementary T4 before/after Shannon TODO. It reuses the frozen seven candidate ASVs and 95 linked biological profiles from `analysis/rerun-controls-2026-08-30/screen/outputs/` in the ecology analysis repository. It changes neither the screen nor its cohort.

Run `python reproduce.py` with NumPy, pandas and SciPy (the paper environment). `before.tsv` contains full-precision pre-removal Shannon, raw richness and depth; `removed_asv_counts.tsv` contains the exact seven removed-ASV counts. For count total N and natural-log Shannon H, sum(c log c) = N(log N - H). Subtract the removed-ASV sum(c log c), reduce N, and calculate post-removal Shannon. Raw richness decreases by the number of removed ASVs with positive count. Spearman correlations compare the 95 before/after profiles.

The full-source derivation in `extract_and_validate.py` additionally reconstructs original entropy, depth and raw richness directly from every canonical count row. It validates all 95 depths, richness values, removed-read counts and removed-ASV counts against the retained screen ledger. Direct-count Shannon and the cached full-precision values agree within 3.6e-15. `input-sha256.json` identifies the exact full inputs in the ecology repository; `summary.json` records the independent validation and unrounded correlations; `t4_before_after.tsv` records each pair.

Extraction usage: `python extract_and_validate.py /path/to/ecology /path/to/output`. The full count table is approximately 1.7 GB and is not duplicated here. The compact inputs suffice to reproduce the reported correlations. A high rank correlation summarizes ordering; individual Shannon values can change (maximum absolute change 0.3094).

On 4 October 2026, the source extraction was replayed against the 30 August screen, whose Trip 4 workbook is byte-identical to data commit 5a17782. The biological cohort, seven candidate ASVs and reported diversity correlations are unchanged; the input manifest now identifies the current screen.
