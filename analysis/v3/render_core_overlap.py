#!/usr/bin/env python3
"""Draw the exact core-genus set intersections, with reproducible provenance."""
import argparse
import csv
import hashlib
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Circle


def intersections(rows):
    names = ("surface", "shallow_subsurface", "root_adjacent")
    sets = {row["compartment"]: set(row["core_genera"].split(",")) for row in rows}
    for row in rows:
        assert len(sets[row["compartment"]]) == int(row["n_core_genera_ge_90pct_sites"])
    result = {}
    for mask in ("100", "010", "001", "110", "101", "011", "111"):
        included = [sets[name] for bit, name in zip(mask, names) if bit == "1"]
        excluded = [sets[name] for bit, name in zip(mask, names) if bit == "0"]
        region = set.intersection(*included) - set.union(set(), *excluded)
        result[mask] = sorted(region)
    assert sum(map(len, result.values())) == len(set.union(*sets.values()))
    return result


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    root = Path(__file__).resolve().parents[2]
    source = root / "analysis/v3/biology_context_corrected_20260909/core_genera_by_compartment.tsv"
    with source.open() as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    regions = intersections(rows)
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 11, "pdf.fonttype": 42})
    fig, ax = plt.subplots(figsize=(6.4, 4.7))
    centres = ((-.58, .32), (.58, .32), (0, -.65))
    for centre, colour in zip(centres, ("#DDAA33", "#4477AA", "#228833")):
        ax.add_patch(Circle(centre, 1.1, facecolor=colour, edgecolor=colour, alpha=.25, lw=1.5))
    positions = {"100": (-1.08,.55), "010": (1.08,.55), "001": (0,-1.18),
                 "110": (0,.75), "101": (-.59,-.48), "011": (.59,-.48), "111": (0,-.05)}
    for mask, xy in positions.items():
        ax.text(*xy, str(len(regions[mask])), ha="center", va="center", fontsize=14)
    ax.text(-.85,1.58,"Surface (82)",ha="center")
    ax.text(.85,1.58,"Shallow-subsurface (69)",ha="center")
    ax.text(0,-1.95,"Root-adjacent (85)",ha="center")
    ax.set(xlim=(-1.85,1.85),ylim=(-2.13,1.82),aspect="equal")
    ax.axis("off")
    fig.tight_layout()
    args.output_dir.mkdir(parents=True,exist_ok=True)
    output = args.output_dir / "fig_core_overlap.pdf"
    fig.savefig(output, metadata={"CreationDate": None, "ModDate": None})
    plt.close(fig)
    digest = lambda p: hashlib.sha256(p.read_bytes()).hexdigest()
    manifest = {"definition": "Core in at least 90% of 60 sites; pooled site-compartment counts; 20 draws at 12865 reads; detection in at least half of draws. Circle areas are schematic.",
                "input": str(source.relative_to(root)), "input_sha256": digest(source),
                "generator": str(Path(__file__).resolve().relative_to(root)),
                "generator_sha256": digest(Path(__file__)), "regions": regions,
                "output_sha256": digest(output)}
    (args.output_dir / "core_overlap_manifest.json").write_text(json.dumps(manifest,indent=2)+"\n")


if __name__ == "__main__":
    main()
