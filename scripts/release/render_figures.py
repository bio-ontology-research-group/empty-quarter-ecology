#!/usr/bin/env python3
"""Recreate every figure cited by the current paper and compare its PDF bytes."""
from __future__ import annotations
import argparse
import filecmp
import subprocess
import sys
import tempfile
from pathlib import Path

FIGURES = ("fig1_landscape.pdf", "fig2_soil_position.pdf", "fig3_function_controls.pdf",
           "fig_core_overlap.pdf", "rain_calendar_refit.pdf", "fig_pma_richness.pdf")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("root", nargs="?", type=Path, default=Path.cwd())
    args = parser.parse_args()
    root = args.root.resolve()
    source = root / "analysis/v3"
    with tempfile.TemporaryDirectory(prefix="eq-current-figures-") as temporary:
        paper = Path(temporary)
        output = paper / "figures"
        for command in (
            [sys.executable, str(source / "render_review_figures.py"), "--landscape", "--output-dir", str(output)],
            [sys.executable, str(source / "render_rain_calendar_refit.py"), "--source", str(source / "rain_calendar_refit_20260909"), "--paper", str(paper)],
            [sys.executable, str(source / "render_core_overlap.py"), "--output-dir", str(output)],
            [sys.executable, str(source / "render_pma_richness.py"), "--output-dir", str(output)],
        ):
            subprocess.run(command, check=True)
        reviewed = root / "empty-quarter-amplicon/figures"
        differences = [name for name in FIGURES if not (reviewed / name).exists()
                       or not filecmp.cmp(output / name, reviewed / name, shallow=False)]
        if differences:
            print("FAIL: regenerated figure PDFs differ: " + ", ".join(differences), file=sys.stderr)
            return 1
        # Compare source and output hashes for the three composite figures and Venn.
        for manifest in ("figure_manifest.tsv", "figure_review_manifest.json", "core_overlap_manifest.json", "pma_richness_manifest.json"):
            if not filecmp.cmp(output / manifest, reviewed / manifest, shallow=False):
                print("FAIL: figure provenance differs: " + manifest, file=sys.stderr)
                return 1
    print("PASS: all six submitted figure PDFs and review manifests are byte-identical")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
