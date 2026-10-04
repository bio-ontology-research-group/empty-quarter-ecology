# ASV-resolution and neighbour-graph sensitivity

The transect association is not an artefact of genus aggregation: at amplicon-sequence-variant resolution the same model gives partial R2 0.3455 to 0.3728 against 0.4007 at genus level (800 ASVs 0.3728 [0.3208, 0.4249]; 2000 ASVs 0.3538 [0.3043, 0.4034]; 2642 ASVs 0.3455 [0.2962, 0.3948]). Both resolutions use the same 630 site-campaign-compartment groups at 60 sites and the corrected site coordinates; the genus figure is the primary fit. 1 of the 630 groups (campaign 3, site 54, Rhizosphere) carries fewer than 2,000 reads on prevalence-filtered ASVs while passing the genus-assigned read threshold; it is retained because the cohort is defined once, by the primary quality control. The residual Moran diagnostic declines with the neighbour count: it is detected at k = [3, 4, 5, 6] and not at k = [8, 10] (Moran I 0.1101 down to -0.0098). The fixed-k residual autocorrelation statement is bounded to the short-neighbourhood scale.

Do not promote the ASV-resolution fit to the primary result; it is a supplementary resolution sensitivity on the same design and inherits every design limit of the primary model, including the collection-order alias. Do not state unqualified residual spatial autocorrelation without naming the neighbour count.

Inputs: `analysis/v2/review/cache/genus_counts.tsv`,
`analysis/v2/review/cache/asv_filt_counts.tsv` and the corrected
site table `analysis/v3/spatial_turnover_rescue/results/site_coordinates.tsv` (checksums in `claim_verdict.json`).

Reproduce from the repository root:

```sh
python analysis/v3/spatial_resolution_sensitivity.py \
  --output-dir analysis/v3/spatial_resolution_sensitivity
```

Evidence files: `asv_resolution_sensitivity.tsv`,
`moran_k_sensitivity.tsv`, `claim_verdict.json`.
