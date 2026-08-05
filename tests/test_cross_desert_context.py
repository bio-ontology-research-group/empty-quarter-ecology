import hashlib
import json
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / "analysis/v3/cross_desert_context"
MAIN = (ROOT / "empty-quarter-amplicon/main.tex").read_text(encoding="utf-8")
SUPPLEMENT = (ROOT / "empty-quarter-amplicon/supplement.tex").read_text(
    encoding="utf-8"
)


def test_cross_desert_outputs_use_independent_biological_units() -> None:
    manifest = json.loads((RESULTS / "run_manifest.json").read_text())
    assert manifest["status"] == "contextual_comparisons_complete"
    assert manifest["cohorts"]["atacama_gradient"]["analysis_sites"] == 16
    assert manifest["cohorts"]["atacama_pit"]["profiles_retained"] == 62
    assert manifest["cohorts"]["atacama_pit"]["sampled_depths"] == 23
    assert "no feature-table merging" in manifest["scope"]

    table = pd.read_csv(RESULTS / "comparison_statistics.tsv", sep="\t")
    assert set(table["independent_unit"]) == {"site", "sampled depth"}
    assert set(table["n_units"]) == {16, 23}


def test_cross_desert_statistics_match_manuscript_claims() -> None:
    table = pd.read_csv(RESULTS / "comparison_statistics.tsv", sep="\t")
    gradient = table[
        (table["question"] == "soil relative humidity versus Shannon diversity")
        & (table["estimate_name"] == "partial_rho")
    ].iloc[0]
    pit = table[
        table["question"] == "community composition among three depth zones"
    ].iloc[0]
    assert round(float(gradient["estimate"]), 2) == 0.68
    assert round(float(gradient["p_value"]), 3) == 0.004
    assert round(float(pit["estimate"]), 2) == 1.46
    assert round(float(pit["p_value"]), 4) == 0.0014
    main = " ".join(MAIN.split())
    supplement = " ".join(SUPPLEMENT.split())
    # The main text explains the effects in plain language; the Supplement
    # preserves the exact adjusted statistic and model label.
    assert "rank correlation $0.68$, $p=0.004$" in main
    assert "depth zones ($p=0.0014$)" in main
    assert "partial $\\rho=0.6753$" in supplement
    assert "pseudo-$F=1.463$" in supplement


def test_cross_desert_checksums_are_current() -> None:
    for line in (RESULTS / "SHA256SUMS").read_text().splitlines():
        expected, name = line.split("  ", 1)
        path = RESULTS / name
        assert hashlib.sha256(path.read_bytes()).hexdigest() == expected


def test_soil_position_figure_adds_descriptive_taxon_context() -> None:
    script = (ROOT / "analysis/v3/make_submission_figures.py").read_text()
    figure = script.split("def make_soil_position_figure", 1)[1].split(
        "def make_function_control_figure", 1
    )[0]
    assert "add_gridspec(2, 3" in figure
    assert "paired_displacement_loadings" in script
    assert "Genera contributing most to position differences" in figure
    assert "do not test each genus separately" in " ".join(MAIN.split())
    manifest = pd.read_csv(
        ROOT / "empty-quarter-amplicon/figures/figure_manifest.tsv", sep="\t"
    )
    assert "paired_composition_loadings" in set(manifest["name"])
