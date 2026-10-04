"""Guard the approved cohort change and the scientific corrections around it."""
import importlib.util
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pytest

from manuscript_paths import PAPER
import manuscript_text as mt

ROOT = Path(__file__).resolve().parents[1]
ANALYSIS = ROOT / "analysis/v3"
sys.path.insert(0, str(ANALYSIS))
spec = importlib.util.spec_from_file_location("rund_renderer", ANALYSIS / "render_review_figures.py")
renderer = importlib.util.module_from_spec(spec)
spec.loader.exec_module(renderer)


def test_figure_uses_matched_common_profile_estimates(tmp_path, monkeypatch):
    source = pd.read_csv(ANALYSIS / "paired_alpha_sensitivity_20260909/paired_alpha_estimands.tsv", sep="\t")
    paired, normalized = renderer.campaign_matched_alpha(source)
    expected = source.query("cohort == 'common_profiles' and estimand == 'campaign_matched'")
    assert len(paired) == 9
    assert set(paired.n_sites) == {60}
    for metric, plotted_metric in (("shannon", "shannon"), ("normalized_shannon", "evenness_h_over_log_hurlbert")):
        reference = expected.query("metric == @metric").set_index("contrast")
        plotted = paired.query("metric == @plotted_metric").set_index("contrast")
        for column, target in (("mean", "mean_difference"), ("ci_low", "bootstrap_ci_low"), ("ci_high", "bootstrap_ci_high")):
            np.testing.assert_array_equal(reference[column], plotted.loc[reference.index, target])
    with monkeypatch.context() as patch:
        patch.setattr(renderer.figures.plt, "close", lambda _: None)
        renderer.figures.make_soil_position_figure(
            paired, normalized,
            pd.read_csv(ANALYSIS / "compartment_composition/compartment_location_results.tsv", sep="\t"),
            pd.read_csv(ANALYSIS / "compartment_composition/paired_displacement_loadings.tsv", sep="\t"),
            tmp_path / "matched.pdf",
        )
        figure = plt.gcf()
        assert figure.axes[1].containers[0].lines[0].get_xdata()[0] == pytest.approx(-0.269, abs=.0005)
        assert figure.axes[2].containers[0].lines[0].get_xdata()[0] == pytest.approx(-0.0425, abs=.00005)
    plt.close(figure)
    with pytest.raises(ValueError, match="Incomplete campaign-matched"):
        renderer.campaign_matched_alpha(source.loc[~((source.cohort == "common_profiles") & (source.contrast == "Deep-Surface"))])


def test_reconciled_claims_and_prose():
    main = mt.text("main")
    supplement = mt.text("supplement")
    source = mt.strip_comments((PAPER / "main.tex").read_text())
    for retired in ("compartment alone", "agnostic to sequencing depth", "all 25 tracked conclusions unchanged",
                    "group 1 [NiFe]", "p=0.00781", "p=0.0432"):
        assert retired not in main, retired
    # No rendered author notes remain (commented-out notes are not rendered).
    assert "\\todo" not in source
    # Current reconciled statements, tied to their sources in other tests.
    assert "campaign-matched" in main.lower()
    estimands = pd.read_csv(ANALYSIS / "paired_alpha_sensitivity_20260909/paired_alpha_estimands.tsv", sep="\t")
    root_deep = estimands.query(
        "cohort == 'common_profiles' and estimand == 'campaign_matched' and metric == 'shannon' and contrast == 'Rhizosphere-Deep'"
    ).iloc[0]
    assert f"{root_deep['mean']:.3f}" == "-0.269"
    assert "-0.269" in supplement and "[-0.492, -0.031]" in supplement
    assert "membrane-bound group 1 or soluble group 3d" in main
    assert "both membrane-bound group 1 and soluble group 3d" in supplement
    assert "permutations restricted within sites (999)" in main
    assert main.index("nearest-sequenced-taxon index") < main.index("Of the 462 predicted MetaCyc pathways")
