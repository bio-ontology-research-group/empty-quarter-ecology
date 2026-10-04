import importlib.util
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT = Path(__file__).parents[1]
spec = importlib.util.spec_from_file_location("review_figures", ROOT / "analysis/v3/make_submission_figures.py")
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)


def test_all_compartment_panels_align_contrast_rows(tmp_path, monkeypatch):
    def read(name):
        return pd.read_csv(ROOT / "analysis/v3" / name, sep="\t")
    with monkeypatch.context() as patch:
        patch.setattr(module.plt, "close", lambda _: None)
        module.make_soil_position_figure(
            read("results/paired_compartment_effects.tsv"),
            read("evenness_decomposition/paired_contrasts.tsv"),
            read("compartment_composition/compartment_location_results.tsv"),
            read("compartment_composition/paired_displacement_loadings.tsv"),
            tmp_path / "figure.pdf",
        )
        figure = plt.gcf()
        assert all(axis.yaxis_inverted() for axis in figure.axes[:3])
        shannon = figure.axes[1].containers[0].lines[0]
        assert np.isclose(shannon.get_xdata()[0], -.206, atol=.001)
        assert shannon.get_ydata()[0] == 0
        assert "Root-adjacent" in figure.axes[0].get_yticklabels()[0].get_text()
    plt.close(figure)
