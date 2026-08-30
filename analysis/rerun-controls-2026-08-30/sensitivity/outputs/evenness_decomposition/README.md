# Post-hoc evenness decomposition

**Status:** `post_hoc_evenness_decomposition_supported`

The lower campaign-averaged Shannon entropy in root-adjacent soil is accompanied by a much clearer lower-evenness signal. Hurlbert expected richness at 25,000 reads shows no root-adjacent versus surface difference and mixed evidence for root-adjacent versus shallow subsurface (the bootstrap interval for the mean excludes zero, but the paired Wilcoxon q is 0.0878). The evenness direction persists in a campaign- and log-depth-adjusted GEE.

## Primary paired results

- **Deep-Surface**: expected-richness difference 142.350 (q=0.06616); evenness-sensitivity difference 0.0120 (95% bootstrap CI 0.0004 to 0.0245; q=0.01338).
- **Rhizosphere-Surface**: expected-richness difference -7.951 (q=0.8309); evenness-sensitivity difference -0.0309 (95% bootstrap CI -0.0440 to -0.0184; q=2.056e-06).
- **Rhizosphere-Deep**: expected-richness difference -150.300 (q=0.08923); evenness-sensitivity difference -0.0429 (95% bootstrap CI -0.0561 to -0.0303; q=3.33e-07).

## Interpretation boundary

The source column named pielou is exactly Shannon divided by log Hurlbert expected richness. Because expected standardized richness replaces observed richness in the denominator, the analysis calls it an evenness sensitivity rather than conventional Pielou evenness.

Limitation: Hurlbert expected richness is unavailable for 55 core-frame profiles below the 25,000-read standard, and normalized evenness is additionally undefined for one single-ASV profile. The corresponding GEE therefore uses 617 of 633 site-campaign-position blocks. All 60 sites contribute to the campaign-averaged paired contrasts, but some site means are based on fewer profiles or campaigns. The decomposition is post hoc.

Permitted wording: In a post-hoc decomposition, the lower paired Shannon distribution in root-adjacent samples was accompanied by a clearer lower normalized-evenness signal than expected-richness signal: expected richness did not differ from surface, and its root-adjacent--shallow evidence was mixed across summaries. This describes the diversity profile and does not identify a root-mediated mechanism.

Prohibited wording: Do not call H/log(E[S_25k]) conventional Pielou evenness, do not describe a causal rhizosphere filter, and do not treat post-hoc decomposition as a preregistered primary endpoint.

## Reproduction

```bash
uv run --python .venv/bin/python analysis/v3/evenness_decomposition_analysis.py
```
