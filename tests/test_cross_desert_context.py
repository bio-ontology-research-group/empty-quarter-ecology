import hashlib
import json
from pathlib import Path

import pandas as pd
from manuscript_paths import PAPER
import manuscript_text as mt


ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / "analysis/v3/cross_desert_context"
MAIN = (PAPER / "main.tex").read_text(encoding="utf-8")
SUPPLEMENT = (PAPER / "supplement.tex").read_text(
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
    # Corrected partial-rank approximation accounts for the nuisance-design rank.
    assert round(float(gradient["p_value"]), 6) == 0.011308
    assert round(float(pit["estimate"]), 2) == 1.46
    assert round(float(pit["p_value"]), 4) == 0.0014
    # The between-desert reanalysis is not part of the current manuscript.
    # Literature on the Atacama is cited, but no reanalysed Atacama statistic
    # or Empty Quarter-versus-Atacama comparison is presented as a result.
    combined = mt.combined()
    for value in ("0.6753", "0.0113", "1.463", "PRJEB39249", "PRJEB17617", "QIITA"):
        assert value not in combined, value
    results = mt.section("main", r"\section*Results", r"\section*Discussion")
    assert "Atacama" not in results


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
    # Panel titles were removed in the October 2026 figure revision; the
    # caption carries the panel description.
    # The Figure 2d caption keeps the descriptive boundary; the Results text now
    # reports the per-genus family (analysis/v3/biology_context, Sep 2026).
    caption = " ".join(MAIN.split()).split("\\caption{Compartment differences", 1)[1].split("\\label", 1)[0]
    assert "Genera with the largest contributions to the three composition contrasts" in caption
    assert "higher CLR abundance" in caption
    manifest = pd.read_csv(
        PAPER / "figures/figure_manifest.tsv", sep="\t"
    )
    # The figure manifest records the paired-displacement loadings that drive
    # panel d, with the checksum of the committed analysis output.
    inputs = manifest[manifest["role"] == "input"]
    loadings = inputs[inputs["file"].astype(str).str.endswith("paired_displacement_loadings.tsv")]
    assert len(loadings) == 1
    source = ROOT / "analysis/v3/compartment_composition/paired_displacement_loadings.tsv"
    assert loadings["sha256"].iloc[0] == hashlib.sha256(source.read_bytes()).hexdigest()
