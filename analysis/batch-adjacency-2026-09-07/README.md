# Library-order neighbour test for processing-batch cross-contamination (2026-09-07)

Decision (Robert Hoehndorf, 7 Sep 2026): the one sequenced Trips 1-3 extraction
blank (M-25-0929, PowerSoil Pro kit blank) carries a soil profile and cannot be
used; test directly whether co-processed Trips 1-3 profiles share DNA, with
Trips 4 and 5 as comparators.

## Method (`batch_adjacency_test.py`)

- Input: canonical `feature-table-trips1-5.tsv` (gzipped copy on Ibex), 1,237
  ecological profiles, relative abundance, Bray-Curtis over all ASVs
  (`outputs/braycurtis_1237.npy`, row order `outputs/braycurtis_profiles.txt`).
- `batch_meta.tsv`: profile, trip, site, compartment, depth, site coordinates
  (data repo `data/metadata/geodata/trip*_geodata.tsv`), library series and
  number (M-number from the July 2025 and Trip 5 samplesheets), DNA
  concentration and kit (`data/release/sample_ledger.tsv`). Trip 2 (24 profiles,
  15 without M-number) and the 29 Trip 3 profiles re-sequenced in the July run
  fall below the 50-profile minimum and are not tested.
- Within each series, profiles are ranked by library number; a pair counts as
  prepared together when ranks differ by at most 2. Different-site pairs only.
  Strata = geographic-distance decile x same-compartment. Statistic = stratum-
  weighted mean of (BC adjacent - BC non-adjacent). Null: 999 permutations of
  library order; one-sided p for adjacent pairs being more similar.
- Per profile: neighbour excess similarity = mean BC to matched non-neighbours
  minus mean BC to neighbours; Spearman against log10 DNA yield (contamination
  would give low-yield profiles a larger excess, i.e. a negative correlation).
- Ibex job 51426151, python/3.9.16 (numpy 1.26.2, scipy 1.11.4), 394 s.

## Result (`outputs/batch_adjacency_results.json`)

| trip | series | profiles | adjacent pairs | delta | p (more similar) | Spearman yield vs excess |
|---|---|---|---|---|---|---|
| 1 | July 2025 M-25 | 325 | 318 | -0.003 | 0.145 | +0.11 (p 0.095, n 238) |
| 3 | Trip 3 run M-23 | 449 | 716 | +0.007 | 0.997 | +0.08 (p 0.12, n 401) |
| 4 | July 2025 M-25 | 177 | 177 | +0.016 | 1.000 | no yields in ledger |
| 5 | Trip 5 run M-25 | 233 | 453 | -0.022 | 0.001 | n 11, not computed |

Trips 1 and 3: no excess similarity between co-prepared profiles; low-yield
profiles are, if anything, less similar to their neighbours than high-yield
ones. Trip 5 shows a small batch signal, which is the trip with per-day blanks
and the extraction-blank screen. e0917_46Dr1 (the profile nearest the blank)
is not unusually close to its own library neighbours (see JSON).

Reported in the ecology supplement S2 and the data descriptor Methods.
