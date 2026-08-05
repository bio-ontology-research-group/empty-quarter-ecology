#!/usr/bin/env python3
"""Render the submitted figures and require byte-identical outputs."""

from __future__ import annotations

import argparse
import difflib
import filecmp
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path


FIGURES = (
    "fig1_landscape.pdf",
    "fig2_soil_position.pdf",
    "fig3_function_controls.pdf",
    "figS_campaign_rainfall.pdf",
    "figure_manifest.tsv",
    "figure_runtime.json",
)


def link(destination: Path, target: Path) -> None:
    destination.parent.mkdir(parents=True, exist_ok=True)
    destination.symlink_to(target.resolve(), target_is_directory=target.is_dir())


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("root", nargs="?", type=Path, default=Path.cwd())
    args = parser.parse_args()
    root = args.root.resolve()
    v3 = root / "analysis/v3"

    with tempfile.TemporaryDirectory(prefix="eq-ecology-figures-") as temporary:
        temporary_root = Path(temporary)
        core = temporary_root / "core"
        output = temporary_root / "figures"
        boundary = temporary_root / "boundary.kml"
        output.mkdir()
        link(core / "cache", root / "analysis/v2/review/cache")
        link(core / "spatial_turnover_rescue", v3 / "spatial_turnover_rescue")
        link(core / "pma_endpoints", v3 / "pma_endpoint_results")
        link(core / "claim_rescue", v3 / "results")
        link(core / "evenness_decomposition", v3 / "evenness_decomposition")
        link(core / "compartment_composition", v3 / "compartment_composition")
        link(core / "distance_decay_turnover", v3 / "distance_decay_turnover")
        shutil.copyfile(
            root / "metadata/geodata/empty_quarter_boundary.kml", boundary
        )

        command = [
            sys.executable,
            str(v3 / "make_submission_figures.py"),
            "--core-dir",
            str(core),
            "--environment-dir",
            str(v3 / "environment_associations"),
            "--picrust-dir",
            str(v3 / "picrust2_ecology"),
            "--rain-dir",
            str(v3 / "rain_pulse_response"),
            "--control-dir",
            str(v3 / "control_audit"),
            "--pma-dir",
            str(v3 / "pma_endpoint_results"),
            "--measured-function-dir",
            str(v3 / "measured_function_summary_results"),
            "--boundary-file",
            str(boundary),
            "--output-dir",
            str(output),
        ]
        subprocess.run(command, check=True)

        reviewed = root / "empty-quarter-amplicon/figures"
        differences = [
            name for name in FIGURES
            if not filecmp.cmp(output / name, reviewed / name, shallow=False)
        ]
        if differences:
            print("FAIL: regenerated figure artifacts differ", file=sys.stderr)
            for name in differences:
                print(f"- {name}", file=sys.stderr)
                if name.endswith(".tsv"):
                    expected_lines = (reviewed / name).read_text(
                        encoding="utf-8"
                    ).splitlines(keepends=True)
                    observed_lines = (output / name).read_text(
                        encoding="utf-8"
                    ).splitlines(keepends=True)
                    sys.stderr.writelines(
                        difflib.unified_diff(
                            expected_lines,
                            observed_lines,
                            fromfile=f"reviewed/{name}",
                            tofile=f"regenerated/{name}",
                        )
                    )
            return 1

    print("PASS: all four figures, manifest, and runtime are byte-identical")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
