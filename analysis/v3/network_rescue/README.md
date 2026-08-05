# Conservative network-claim rescue

This directory replaces the thresholded-Spearman network analysis with a
reproducible, compositional sensitivity analysis. It does not modify either
manuscript.

The primary analysis:

1. parses only core sites 1--60;
2. sums sequencing replicates within campaign × site × compartment;
3. takes the exact campaign/site intersection shared by surface, deep, and
   rhizosphere samples;
4. retains taxa meeting the same prevalence rule in every compartment and
   ranks them by pooled mean relative abundance;
5. applies a fixed-pseudocount CLR, centers taxa within campaign, and fits
   GraphicalLasso models with one common regularization parameter;
6. calibrates regularization against null data made by independently permuting
   every taxon within campaign;
7. requires full-data selection, site-cluster bootstrap stability, sign
   consistency, and low edge-wise null frequency;
8. repeats the density comparison across taxon counts and regularization
   multipliers.

The numerical nonzero threshold is an absolute partial correlation of
`1e-3`, safely above the recorded `2e-4` solver tolerance, so solver residuals
are not counted as edges.

Campaign-specific models are not forced when the matched sample sizes fail the
predeclared observation rule. Campaigns are never combined into wet/dry bins.
All reported edges are descriptive conditional associations.

Run from any directory:

```bash
uv run \
  --with 'numpy==2.1.*' \
  --with 'pandas==2.2.*' \
  --with 'scikit-learn==1.5.*' \
  python /path/to/empty-quarter/analysis/v3/network_rescue/run_network_rescue.py \
  --project-root /path/to/empty-quarter \
  --output-dir /path/to/results
```

The CLI takes explicit `--project-root`, `--input-table`, and `--output-dir`
arguments and is safe to call from a Nextflow task directory. A minimal DSL2
entry point is provided in `main.nf`; the repository workflow environment is
used for its Python dependencies:

```bash
nextflow run analysis/v3/network_rescue/main.nf -with-conda \
  --project_root /path/to/empty-quarter \
  --outdir /path/to/results
```

Outputs:

- `cohort_accounting.tsv`: sample/group attrition and exact matching;
- `cohort_groups.tsv`: group-level read depth, QC, and matched inclusion;
- `taxa_selection.tsv`: prevalence, abundance rank, and primary inclusion;
- `alpha_calibration.tsv`: observed versus independently permuted null models;
- `network_edges.tsv`: every tested pair and each stability gate;
- `network_metrics.tsv`: primary compartment summaries;
- `bootstrap_metrics.tsv`: paired site-cluster bootstrap summaries;
- `sensitivity_metrics.tsv`: taxon-count and regularization sensitivity;
- `compartment_comparisons.tsv`: matched uncertainty and robustness gates;
- `campaign_feasibility.tsv`: why campaign models were or were not estimated;
- `network_rescue.tsv` and `claim_verdict.json`: machine-readable disposition;
- `parameters.json`: checksums, seeds, thresholds, cohort counts, and software;
- `README.md`: generated human-readable result summary.

For fixed input bytes, parameters, and dependency versions, outputs contain no
timestamps and are deterministic.
