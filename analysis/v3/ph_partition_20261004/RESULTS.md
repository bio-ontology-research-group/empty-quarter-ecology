# pH partition results and author interpretation

The variance-denominator problem is fixable using the existing processed data. A coherent partition shows a small unique pH association and a larger component shared with geography. The observational sampling design still limits separation of causal pH effects from spatially covarying environmental processes.

| Estimand | Unique pH | Unique geography | Shared association | Unexplained |
|---|---:|---:|---:|---:|
| 560 groups, equal total weight per site | 0.312% | 6.778% | 4.763% | 88.147% |
| 60 site-average adjusted profiles | 1.382% | 13.564% | 23.856% | 61.197% |
| 560 groups, equal weight per group (sensitivity) | 0.301% | 6.592% | 4.168% | 88.939% |

All entries are raw in-sample fractions after the specified nuisance adjustment. Each row has its own common denominator and sums to 100%. The grouped row includes within-site variation; the site-average row describes spatial variation in average adjusted profiles. The latter is the relevant companion for the previous site-average geography discussion. Its combined pH/geography association is 38.803%, with geography alone accounting for 37.420%. The shared fraction is a covariance component; assigning that component to a pH-mediated pathway would require independent causal evidence.

The original 37.54% and 19.78% values use two different residual community responses and denominators. Subtracting them or computing their relative reduction does not identify an explained fraction of the same response. The new site-average calculation also adjusts both predictor sets against the same campaign-by-compartment nuisance design before averaging.

pH and geography remain strongly aligned between sites: the adjusted geographic columns account for 84.39% of adjusted site-average pH variation (VIF 6.41). At the grouped scale this association is 53.44% (VIF 2.15). The scale distinction explains why the grouped and site-average results differ. The site-fixed sensitivity gives pH a within-site partial R2 of 0.307%, equivalent to 0.215% of the grouped primary denominator.

Deleting any single site leaves unique pH at 0.276–0.385% in the grouped model and 1.248–1.626% in the site-average model. Corresponding geographic fractions are 6.303–7.074% and 12.469–14.347%. Deleting only the highest-pH group (trip 4, site 60, rhizosphere, pH 9.98, one admitted pH specimen) gives grouped fractions 0.406% unique pH, 6.272% unique geography, 5.363% shared and 87.958% unexplained.

Whole-site resampling and contiguous spatial-block resampling quantify sensitivity to the sampled sites. Each uses 1,999 draws, retaining all repeated groups and re-estimating nuisance fits. At the site-average scale, the shared fraction's percentile range is 15.13–32.14% under independent-site sampling, and 0.22–31.38% under 15-site blocks. The same 15-site blocks span a median 253 km and involve only four drawn blocks. Spatial sensitivity and edge effects are therefore material, and these ranges have no certified confidence coverage.

Raw unique-pH fractions are biased upward in resampling: for the grouped model the independent-site bootstrap median is 0.521% against a point estimate of 0.312%; for site averages it is 2.322% against 1.382%. Their basic centered diagnostic ranges are −0.297–0.291% and −1.523–1.432%, respectively. Negative endpoints are retained as bias diagnostics. Percentile bounds above zero must not be turned into significance evidence. No p-value, mediation percentage, adjusted-R2 or causal pH estimate is reported.

Sparse nuisance categories disappear in 176, 327, 740 and 923 grouped draws for site clusters and block lengths 5, 10 and 15 respectively. The pH/geographic increments remain identifiable in all 7,996 draws for each estimand. All 15,992 complete bootstrap partitions are saved. The larger loss of campaign categories and unequal endpoint sampling under longer blocks are part of the sensitivity limitation, not hidden fit failures.

The package reconstructs the exact 560-group, 200-genus response and ordered taxon selection used by the archived pH code from actual pinned processed inputs. Six tests pass. Clean-environment regeneration produces byte-identical outputs, including every bootstrap draw. The scope begins at processed pH/count/coordinate inputs; laboratory specimen identity, raw-read generation and causal identification remain separate scientific questions.

For the manuscript, report the between-site partition with its defined response and the grouped sensitivity alongside it if useful. Retain the conclusion that pH covaries with geography and has a small conditional association. A replicated intervention or geographically crossed pH sampling would add the independent variation needed to estimate causal pH effects.
