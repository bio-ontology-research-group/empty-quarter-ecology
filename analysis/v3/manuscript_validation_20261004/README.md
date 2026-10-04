# Independent manuscript checks, 4 October 2026

`verify_recruitment.py --coverm /path/to/coverm_profiles.tar.gz --output analysis/v3/manuscript_validation_20261004` sums the original CoverM recruitment of all 991 catalogue genomes in each of 150 libraries. The mean is 28.518623%, confirming the manuscript's 28.5%. Using only the 975 genomes in the trait analysis would answer a different question. The output JSON pins the original archive hash; the TSV preserves every library value.

`verify_archaea.py --root /path/to/ecology` reconstructs the mean archaeal read fraction directly from the canonical feature table, SILVA domain labels and 1,227 core-profile identifiers. It checks all sample depths against the alpha-diversity cache. The mean share is 0.0239404%, confirming 0.02% at the manuscript's precision. Its JSON records all input hashes. These inputs are unchanged in data commit 5a17782.
