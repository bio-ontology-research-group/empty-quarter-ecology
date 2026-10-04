"""Semantic checks for the current core-set graphic and figure replay scope."""
import importlib.util
import json
import re
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def test_core_intersections_are_exclusive_and_exhaustive():
    spec = importlib.util.spec_from_file_location("core_overlap", ROOT / "analysis/v3/render_core_overlap.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    rows = [{"compartment": name, "core_genera": taxa, "n_core_genera_ge_90pct_sites": "3"}
            for name, taxa in (("surface", "a,b,c"), ("shallow_subsurface", "b,c,d"), ("root_adjacent", "c,d,e"))]
    regions = module.intersections(rows)
    assert regions == {"100": ["a"], "010": [], "001": ["e"], "110": ["b"],
                       "101": [], "011": ["d"], "111": ["c"]}


def test_replay_covers_every_submitted_figure():
    from manuscript_paths import PAPER
    spec = importlib.util.spec_from_file_location("figure_replay", ROOT / "scripts/release/render_figures.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    sources = "\n".join((PAPER / name).read_text() for name in ("main.tex", "supplement.tex"))
    sources = re.sub(r"(?m)^\s*%.*$", "", sources)
    figures = re.findall(r"\\includegraphics(?:\[[^]]*\])?\{figures/([^}]+)\}", sources)
    assert set(figures) == set(module.FIGURES)
    data = json.loads((PAPER / "figures/core_overlap_manifest.json").read_text())
    assert len(data["regions"]["111"]) == 59
    assert sum(len(v) for k,v in data["regions"].items() if k[0] == "1") == 82
    assert sum(len(v) for k,v in data["regions"].items() if k[1] == "1") == 69
    assert sum(len(v) for k,v in data["regions"].items() if k[2] == "1") == 85


def test_repository_verifier_detects_inline_todos_but_ignores_comments(tmp_path):
    import csv
    spec = importlib.util.spec_from_file_location("release_verify", ROOT / "scripts/release/verify_repository.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    paper = tmp_path / "empty-quarter-amplicon"
    (paper / "figures").mkdir(parents=True)
    (paper / "supplement.tex").write_text("")
    for name in ("main.pdf", "supplement.pdf", "figures/fig1_landscape.pdf",
                 "figures/fig2_soil_position.pdf", "figures/fig3_function_controls.pdf",
                 "figures/rain_calendar_refit.pdf", "figures/fig_core_overlap.pdf", "figures/fig_pma_richness.pdf",
                 "figures/pma_richness_manifest.json",
                 "figures/figure_review_manifest.json", "figures/core_overlap_manifest.json",
                 "figures/figure_runtime.json"):
        (paper / name).write_text("")
    (paper / "figures/figure_manifest.tsv").write_text("role\tfile\tbytes\tsha256\n")
    (paper / "main.tex").write_text(r"\todo[inline]{unfinished figure}")
    failures = []
    module.verify_manuscript(tmp_path, failures)
    assert failures == ["active manuscript contains unresolved TODO markers"]
    (paper / "main.tex").write_text("% " + r"\todo[inline]{historical comment}")
    failures = []
    module.verify_manuscript(tmp_path, failures)
    assert failures == []


def test_catalogue_recruitment_uses_all_991_genomes():
    import pandas as pd
    source = ROOT / "analysis/v3/manuscript_validation_20261004"
    data = pd.read_csv(source / "library_recruitment.tsv", sep="\t")
    summary = json.loads((source / "recruitment_summary.json").read_text())
    assert len(data) == 150 and set(data.genomes) == {991}
    assert abs(data.recruited_pct.mean() - summary["mean_recruited_pct"]) < 1e-10
    assert round(summary["mean_recruited_pct"], 1) == 28.5
