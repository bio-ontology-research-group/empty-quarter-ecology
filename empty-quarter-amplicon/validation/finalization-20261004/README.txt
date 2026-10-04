Review package: ecology manuscript, 4 October 2026

Baseline: Rund's 2 October Overleaf revision a5d5623.
Repository: https://git.overleaf.com/6a7422eef2624277cf100b3b

OPEN_VERIFICATION.txt identifies the physical well-to-well result still needing
Rund's source files before submission. VALIDATION.txt records the checks performed.

CHANGES.txt groups every substantive correction, including the two commits
already made after Rund's revision. source-changes.patch records exact changes
to manuscript, bibliography, table and generated LaTeX source and build settings.

main-redline.pdf and supplement-redline.pdf are flattened latexdiff comparisons:
blue underlined text is added; red struck-through text is removed. Existing
broken references and citation keys in removed text are shown as [old ref] or [old cite]; the source patch retains their keys.
The reference list is rendered from the current bibliography; all bibliography
source edits are preserved in source-changes.patch.

figure-comparison.pdf shows changed and added graphics. The three original
main composite figures regenerate byte for byte; their caption edits appear
in the text redlines. The core-genus overlap figure is new, rainfall reflects
the corrected Site 52 daily series, and the PMA plot has a reproducible vector
renderer with the same nine pairs and three group means.

Final manuscript PDFs and editable source are at the project root. The ecology
GitHub repository mirrors these files under empty-quarter-amplicon/ and contains
the analysis code, result tables and manifests. The LaTeX comparison generator
is scripts/release/build_rund_redlines.py in that repository. Export the baseline
with git archive a5d5623 and provide it via --baseline, the built current paper
via --paper, and a separate temporary directory via --output.
