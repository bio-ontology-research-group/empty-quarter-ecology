# Physical well-to-well analysis

This analysis tests physical positions on preparation plates. It is distinct
from the library-number adjacency analysis in the neighbouring directory.
The user supplied `well_analysis.zip` on 5 October 2026. Its original bytes
are retained as `author_source.zip`; the extracted plate map is
`well_positions.csv`. The archive contains the original scripts and outputs.

## Reproduce

Use the ecology repository's locked Python environment and the data repository
commit specified by `DATA_REPOSITORY.lock`. Install its checksum-pinned bulk
canonical feature table, then run from the ecology repository root:

```sh
make well-analysis DATA_REPO=/path/to/pinned/data-repository PYTHON=/path/to/python
```

The two Python scripts also accept explicit input paths (`--help`).
`01_build_bc.py` uses the deposited plate map, current batch metadata and all
ASVs in the canonical feature table. `02_neighbor_permtest.py` reproduces the
supplied statistical procedure with explicit paths and a JSON result added.
The source archive, plate map and all inputs are identified by checksums in
`verification.json` and `outputs/input_provenance.json`. Python and dependency
versions are in `verification.json`; they match the existing environment.

## Design and selection

The plate map contains 344 occupied positions on six 96-well plates and 298
unique derived profile identifiers. As in the supplied analysis, profiles
listed more than once retain their first listed position (46 additional rows
are discarded). Ten profile identifiers do not match the biological,
quality-filtered metadata. Both the repeated profiles and unmapped positions
are exported for inspection. The remaining 288 profiles comprise 259 T1,
one T2 and 28 T3 profiles. A campaign needs at least ten profiles to be tested.

Bray–Curtis distances use relative abundance over all ASVs. Pairs must belong
to the same plate and different sites. Matching strata combine within-campaign
geographic-distance decile with same versus different soil compartment (20
strata). The primary adjacency definition includes an edge or corner; integer
well coordinates with radius 1.5 give exactly the same neighbours as radius
sqrt(2). The secondary definition includes edges only. The statistic is the
mean adjacent-minus-nonadjacent Bray–Curtis contrast within strata, weighted
by adjacent-pair count. Each of 999 permutations shuffles sample-to-position
assignments within plates. The one-sided p-value counts null statistics no
greater than observed, with the add-one correction. The shared random stream
starts at seed 42 and processes campaign and adjacency definitions in the
supplied order. Optional per-well DNA-yield diagnostics retain the supplied
matching fallback and pseudocount rules; they do not determine the global test.

## Current results

| Campaign | Neighbours | Profiles | Adjacent pairs | Contrast | One-sided p |
|---|---|---:|---:|---:|---:|
| T1 | Edge or corner | 259 | 441 | +0.00514797 | 0.866 |
| T1 | Edge only | 259 | 190 | +0.00330777 | 0.740 |
| T3 | Edge or corner | 28 | 69 | -0.01102049 | 0.173 |
| T3 | Edge only | 28 | 34 | -0.01096703 | 0.203 |

The T1 headline remains +0.005 and p=0.87 at the manuscript's precision.
The supplement now specifies same/different-compartment matching, the
first-position rule and within-plate permutations with seed 42. The main
paper's values and conclusion are unchanged.

## Source verification

The original archive's summary and all numerical result tables were reproduced
using its deposited distances and well metadata. Its T1 p-value was 0.867.
The current replay uses the corrected coordinates from the pinned metadata;
two mapped profiles (T1 and T3 at site 52) have updated coordinates. This
changes geographic matching, but not the regenerated Bray–Curtis matrix,
which is exactly identical to the supplied one (maximum difference zero).
`original_replay/` records that verification; `outputs/` contains the current
results. Historical coordinates are retained only inside the source archive.

`tests/test_physical_well_analysis.py` verifies provenance, coordinate joins,
selection counts and the full 999-permutation replay, including manuscript
rounding. This closes the physical-well item in the earlier finalization
review's `OPEN_VERIFICATION.txt`. Upstream sequencing and other unrelated
submission requirements remain outside this analysis.
