# Contaminant screen rerun vs published (generated 2026-08-30T04:22:08.223807+00:00)

Mode: per-trip; Trip 4 screen: run

## Candidate ASVs

- Published (Trip 5 blanks EB1-EB17 vs 217 profiles): 351
- Rerun screen Trip4: 7 candidates; overlap with published set 0; new 7; dropped 351
- Rerun screen Trip5: 351 candidates; overlap with published set 351; new 0; dropped 0
- Trip 5 screen reproduces the published candidate set exactly: YES
- Top genera among Trip 4 / pooled candidates: NA (4), Nesterenkonia (1), Brevibacterium (1), Alcaligenes (1)
- Top genera among candidates NOT in the published set: NA (4), Nesterenkonia (1), Brevibacterium (1), Alcaligenes (1)

## Profiles

- Published: 217 filtered profiles; reads removed median 0.40% (IQR 0.16-0.99%), max 56.60%; pooled 2.19%
- Rerun Trip5: 217 filtered profiles, 217 lose >= 1 read, reads removed median 0.40% (IQR 0.16-0.99%), max 56.60%; pooled 2.19%; profiles below 25,000 reads after filtering: 0
- Rerun Trip4: 95 filtered profiles, 31 lose >= 1 read, reads removed median 0.00% (IQR 0.00-0.01%), max 14.31%; pooled 0.14%; profiles below 25,000 reads after filtering: 0
- Trip 5 profiles whose removed-read count differs from the published run: 0

## Headline numbers the manuscript reports (Methods, 'Assay controls')

| quantity | manuscript / published | rerun |
|---|---|---|
| candidate ASVs | 351 | {'Trip4': 7, 'Trip5': 351} |
| linked profiles filtered | 217 (Trip 5) | Trip 5 217, Trip 4 95 |
| training blanks | 17 | Trip 5 17, Trip 4 6 |

## Downstream (25 tracked conclusions)

Not rerun by this script. The published downstream stage (scripts/controls/build_control_sensitivity_inputs.py + run_control_ecology_sensitivity.sh) hard-codes 217 Trip 5 profiles and V-prefixed IDs; extending it to Trip 4 needs a patched copy. Compare analysis/v3/control_sensitivity/headline_result_sensitivity.tsv once that is run.
