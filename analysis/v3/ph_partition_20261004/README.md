# pH and geography association partition, 4 October 2026

This self-contained analysis partitions the association of archived-soil pH and geography with amplicon community composition. It regenerates from frozen **processed** pH measurements, genus counts, profile metadata and corrected site coordinates included here. It leaves the earlier pH analysis unchanged. The output is an observational, raw descriptive variance partition; each reported fraction has an explicitly shared denominator.

## Reproduce

From any working directory, with [uv](https://docs.astral.sh/uv/) available:

```bash
bash /path/to/analysis/v3/ph_partition_20261004/reproduce.sh /tmp/ph-partition-new-results
```

The destination must not exist. The wrapper creates an isolated temporary CPython 3.11.14 environment, installs the six exact versions in `requirements.lock.txt`, runs six tests, and regenerates both partitions and 1,999 bootstrap replicates for each of four resampling designs. All numerical work runs with one BLAS thread and seed 20261004 (resampling seed plus block length). Installation requires access to Python/package downloads or an existing uv cache; analysis itself has no network or external data dependency. The temporary environment is removed on exit. `results/` contains the checked reference outputs; compare a new run with `diff -r results /tmp/ph-partition-new-results` from this directory. The clean reproduction log and comparison result are in `validation/`.

With the pinned dependencies already installed, the numerical command is:

```bash
python analyze.py --output /tmp/ph-partition-new-results --replicates 1999
```

## Frozen input boundary and provenance

`input_manifest.json` identifies every input/provenance file with a SHA256, its original source-relative path, the original repository HEAD and its complete dirty status. The analysis verifies packaged hashes before reading inputs. The package includes actual bytes, and requires none of the absolute source paths from earlier manifests. Genus counts are deterministically gzip-compressed; the manifest also gives their uncompressed source hash.

The upstream source is the ecology repository at b5cefb2 plus existing local September corrections. HEAD alone is insufficient to reproduce those corrections. The exact existing `ph_ecology_analysis.py`, earlier pH summary, pH ingest summary, input manifest and repository reproduction guide are archived in `provenance/`; they document origin and are not imported or executed. The historical `DATA_REPOSITORY.lock` is archival context, not the dependency pin of this standalone package. Every operational input is packaged and pinned independently.

The executable chain starts with admitted `EQ-PH-SHARED-v1.0.0` measurements, processed per-profile genus counts/metadata and the corrected 60-site coordinate table. It joins pH to ecological profiles at campaign, site and compartment, averages admitted pH values within each group, sums count profiles within the same group, and retains groups with at least 2,000 genus-assigned reads. The frozen cohort contains 560 groups at 60 sites and 1,089 selected count profiles. The code checks group pH values against the preserved September group table. All profiles in an eligible group contribute; matching labels are not asserted to establish physical-specimen pairing. Archived pH source rows and original laboratory provenance remain those described by the retained pH ingest report.

Within this fixed cohort, genera detected in at least 20% of groups are ranked by mean relative abundance, retaining the first 200 with stable tie ordering. Each response coordinate is `log(count + 0.5)` minus the group mean log count across these 200 genera. Features, filtering, transformation and response values remain fixed in every fit and resample. The output cohort and ordered taxon list are saved. This is reproducibility from processed matrices; raw sequencing, genus inference, pH workbook admission, laboratory provenance and knowledge-graph construction are upstream boundaries.

## Grouped primary estimand

Each group has weight `1 / number_of_retained_groups_at_its_site`, so each of 60 sites has total weight one. The response is the same 560 by 200 CLR matrix in all models. The nuisance design N contains an intercept and all observed campaign-by-compartment categories (rank 15). P is standardized group pH; G contains the standardized transect coordinate and its square. Standardization is fixed from the original cohort, and does not change polynomial column spaces. The four nested weighted least-squares fits use N, N+P, N+G and N+P+G. Their weighted residual sums of squares are S_N, S_P, S_G and S_PG. The fractions are:

- unique pH: `(S_G - S_PG) / S_N`;
- unique geography: `(S_P - S_PG) / S_N`;
- shared association: `(S_N - S_P - S_G + S_PG) / S_N`;
- unexplained: `S_PG / S_N`.

They sum to one. Shared fractions remain signed, including under suppression; no clipping is applied. All four models retain identical rows, columns, weighting and nuisance treatment. The denominator contains both within-site and between-site composition variation after nuisance adjustment. Fractions are raw in-sample associations. No adjusted-R2 or unbiased population-explained-variance claim uses a fabricated effective sample size. The same-cohort equal-group-weight sensitivity is reported separately.

## Between-site companion estimand

To address the geographic pattern of site-average composition, the companion first residualizes the group CLR response, pH, linear transect and quadratic transect against the **same** weighted N design. It then averages each residual object within each site. Four regressions of the fixed, centered 60 by 200 site-average adjusted response use an intercept, adjusted site-average pH, adjusted site-average geographic columns, or both. The partition uses their common intercept-only residual sum of squares. Thus the predictors are averaged nuisance-adjusted pH/geographic columns, and the numerator and denominator belong to the same site response. Nuisance adjustment is re-estimated within every resample.

The grouped denominator and site-average denominator answer different questions. Their fractions cannot be subtracted or compared as a reduction. The legacy 37.54% and 19.78% statistics fitted geography to two different residual responses and used two different total sums of squares. Their relative change is not a partition of one community response. This package supplies the coherent grouped and between-site partitions directly.

## Resampling and sensitivity

Every resample selects whole sites with all their groups and fixed within-site weights. Repeated selection multiplies a site's sufficient statistics by its draw count. Campaign-position nuisance coefficients are re-estimated for each replicate; selected taxa and transformations stay fixed. Tests verify equivalence to explicit repeated-row fitting.

Four designs each use 1,999 draws: independent site clusters, and overlapping non-circular blocks of 5, 10 or 15 consecutive sites ordered by corrected transect coordinate. For a block length L, starting indices are sampled uniformly from `0..60-L`; blocks are concatenated and trimmed to 60 site draws. The transect endpoints are never made adjacent. Median spans are about 66, 155 and 253 km, with exact ranges in the summary. These block sizes are declared sensitivity settings, not estimated optimal correlation ranges. Block sampling gives edge sites lower expected inclusion and provides only 12, 6 or 4 drawn blocks; sparse campaign-position categories may disappear. Missing nuisance columns use a pseudoinverse, and diagnostics verify that pH and geographic increments retain ranks 1, 2 and 3 above the supported nuisance space. Missing-category counts and unidentified draws are reported.

The saved percentile ranges describe resampling stability. Raw explained fractions are nonnegative and biased upward; their percentile ranges can exclude zero even when evidence for a conditional association is weak. The additional basic ranges, `2*point - opposite_percentile`, are untruncated bias diagnostics and may be negative. Neither range is certified as a coverage-calibrated confidence interval or a significance test. Site-cluster sampling assumes sites can represent independent draws; contiguous blocks explore local spatial dependence under approximately homogeneous ordering assumptions. A single environmentally structured transect, irregular spacing and campaign coverage limit generalization under either scheme. Every bootstrap draw is retained, permitting alternatives without hiding this uncertainty.

Further outputs report geographic collinearity with pH, all 60 leave-one-site-out fits, removal of the highest-pH group (with its remaining site groups reweighted), and a site-fixed-effect pH fit. The latter targets within-site pH association and has its own residual denominator, also reported on the grouped primary denominator. None of these steps resolves causal direction, soil/pH covariation, measurement error, missing-campaign selection, laboratory/campaign confounding or uncertain original-specimen matching.

## Validation and interpretation

Six tests cover partition arithmetic and signed shared fractions, explicit clustered-bootstrap multiplicity, nuisance/weight invariance, non-circular block sampling and fixed-seed draws, input/cohort/feature/weight/design consistency, and re-estimation of common nuisance adjustment for site-average resampling. The main run separately compares sufficient-statistic fits against direct weighted least squares. Clean-environment regeneration is compared byte-for-byte at the output boundary.

`RESULTS.md` reports the numerical interpretation and limitations. `MANUSCRIPT-CANDIDATE.md` contains proposed prose for author review, separately from the manuscripts. The analysis can fix denominator arithmetic and estimate conditional associations. Separating causal pH effects from spatially covarying environmental processes would require additional independent variation, such as replicated pH manipulations or sampling that crosses pH levels with geographic position.

Methodological context: the [vegan variation-partitioning documentation](https://vegandevs.github.io/vegan/reference/varpart.html) describes partitions, signed shared components and the bias of raw fractions. Its ordinary adjusted-R2 default is not adopted for this clustered weighted design. The [R boot documentation](https://stat.ethz.ch/R-manual/R-devel/library/boot/html/tsboot.html) describes block resampling and distinguishes endpoint wrapping; this implementation uses explicitly non-circular blocks. Consulted 4 October 2026.
