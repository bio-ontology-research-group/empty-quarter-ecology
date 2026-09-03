# Taxon and predicted-pathway context

- Status: `descriptive_context_complete`
- Cohort: 1227 core-site profiles (sites 1-60), 351472 ASVs, 1604 named genera
- Genus sums recomputed from the canonical feature table reproduce `analysis/v2/review/cache/genus_counts.tsv` (max |difference| 0 reads)

## Leading phyla (mean share of total reads)

- Pseudomonadota: 26.9 % (prevalence 100 %)
- Bacillota: 22.5 % (prevalence 100 %)
- Actinomycetota: 18.7 % (prevalence 99 %)
- Chloroflexota: 8.6 % (prevalence 99 %)
- Bacteroidota: 7.1 % (prevalence 99 %)
- Planctomycetota: 3.2 % (prevalence 98 %)
- Gemmatimonadota: 2.7 % (prevalence 98 %)
- Acidobacteriota: 2.3 % (prevalence 95 %)

## Leading genera (all core profiles)

- Domibacillus (Bacillota): 6.22 %
- Massilia (Pseudomonadota): 3.55 %
- Bacillus (Bacillota): 2.97 %
- Microvirga (Pseudomonadota): 2.41 %
- Flavisolibacter (Bacteroidota): 2.05 %
- Noviherbaspirillum (Pseudomonadota): 1.92 %
- Streptomyces (Actinomycetota): 1.80 %
- Herpetosiphon (Chloroflexota): 1.52 %
- Paenibacillus (Bacillota): 1.38 %
- Lysinibacillus (Bacillota): 1.32 %

## Transect replacement (200 primary genera, site-level CLR)

- Supported route correlations (BH q < 0.05 over 200 tests): 124 (47 increase eastward, 77 decrease eastward)
- Strongest eastward decreases: Cellulomonas, Blastocatella, Marmoricola, Ellin6055, Roseisolibacter
- Strongest eastward increases: Halalkalibacter, Polygonibacillus, Ammoniphilus, Bacillus, Pseudalkalibacillus

## Genus-set overlap (rarefied presence)

- Empty Quarter west third vs Empty Quarter east third: 368 vs 317 genera, 267 shared, Jaccard 0.639; 50/50 leading genera of A detected in B
- Empty Quarter west third vs Empty Quarter central third: 368 vs 364 genera, 331 shared, Jaccard 0.825; 50/50 leading genera of A detected in B
- Empty Quarter central third vs Empty Quarter east third: 364 vs 317 genera, 275 shared, Jaccard 0.677; 50/50 leading genera of A detected in B
- Empty Quarter all core profiles vs Atacama pit all depths (PRJEB39249): 364 vs 76 genera, 43 shared, Jaccard 0.108; 17/50 leading genera of A detected in B
- Empty Quarter surface vs Atacama pit 2.5-10 cm (PRJEB39249): 336 vs 61 genera, 33 shared, Jaccard 0.091; 15/50 leading genera of A detected in B

## Environmental gradients along the route (Spearman rho with route position)

- 49-month mean air temperature: rho = +0.98 (n = 60; west third 27.89, east third 30.07)
- 49-month mean monthly rainfall: rho = +0.55 (n = 60; west third 2.16, east third 3.67)
- 49-month mean relative humidity: rho = +0.92 (n = 60; west third 24.72, east third 31.82)
- archived-soil pH (admitted measurements, site mean): rho = +0.82 (n = 60; west third 7.85, east third 8.23)

## Predicted pathway classes (share of predicted pathway abundance)

- biosynthesis: 66.8 % (231 pathways)
- energy_central_metabolism: 17.2 % (55 pathways)
- degradation_utilization: 10.7 % (135 pathways)
- other: 5.3 % (41 pathways)

## Autotrophy-related pathways

- rank 44: CALVIN-PWY Calvin-Benson-Bassham cycle: 0.59 %
- rank 122: P23-PWY reductive TCA cycle I: 0.41 %
- rank 284: CODH-PWY reductive acetyl coenzyme A pathway: 0.04 %
- rank 304: PWY-5392 reductive TCA cycle II: 0.02 %
- rank 360: PWY-5743 3-hydroxypropanoate cycle: 0.00 %
- rank 362: PWY-5744 glyoxylate assimilation: 0.00 %
- rank 460: PWY-5789 3-hydroxypropanoate/4-hydroxybutanate cycle: 0.00 %

## Permitted wording

- Relative abundances are shares of reads and describe the marker-gene profile, not cell counts.
- The pathway classes are keyword groupings of MetaCyc descriptions defined in `run_manifest.json`; they are not MetaCyc ontology classes.
- The Atacama comparison is a presence-based overlap under a shared detection rule; it does not merge feature tables and does not correct for primer or DNA-fraction differences.
- Route correlations describe which genera change along the transect; they do not identify a cause.
