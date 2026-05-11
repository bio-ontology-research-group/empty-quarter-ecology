#!/usr/bin/env python3
"""Render wind-Mantel sweep heatmaps:
  - sensitivity (window x angle, per dust threshold)
  - score-type comparison (window x score_type, per compartment)
  - distance stratification (window x stratum, per compartment)
"""
from __future__ import annotations

from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

REPO = Path(__file__).resolve().parents[1]
CACHE = REPO / "cache" / "wind_dispersal"
FIG = REPO / "figures"
FIG.mkdir(exist_ok=True)


def heatmap(ax, df, idx, col, val, title, vmin=-0.4, vmax=0.4):
    pv = df.pivot_table(index=idx, columns=col, values=val).round(3)
    im = ax.imshow(pv.values, aspect="auto", cmap="RdBu_r", vmin=vmin, vmax=vmax)
    ax.set_xticks(range(len(pv.columns)))
    ax.set_xticklabels(pv.columns)
    ax.set_yticks(range(len(pv.index)))
    ax.set_yticklabels(pv.index)
    ax.set_title(title, fontsize=10)
    ax.set_xlabel(col)
    ax.set_ylabel(idx)
    for i in range(pv.shape[0]):
        for j in range(pv.shape[1]):
            v = pv.values[i, j]
            if pd.notna(v):
                ax.text(j, i, f"{v:.2f}", ha="center", va="center",
                        fontsize=7, color=("white" if abs(v) > 0.25 else "black"))
    return im


def main():
    df = pd.read_csv(CACHE / "sweep_mantel_full.tsv", sep="\t")
    print(f"sweep rows: {len(df)}", flush=True)

    # --- Figure 1: sensitivity heatmaps (window x angle), one per threshold,
    #               for compartment=surface, score=max, stratum=all
    fig, axes = plt.subplots(1, 3, figsize=(13, 4), constrained_layout=True)
    for ax, thr in zip(axes, [5.0, 7.0, 10.0]):
        s = df[(df["score_type"] == "max") & (df["stratum"] == "all") &
               (df["compartment"] == "surface") & (df["dust_thr"] == thr)]
        # median across trips
        s_med = s.groupby(["window_days", "angle_tol"])["r_part"].median().reset_index()
        im = heatmap(ax, s_med, "window_days", "angle_tol", "r_part",
                     f"Surface | dust_thr={thr} m/s")
    fig.suptitle("Sensitivity sweep: median partial Mantel r(BC, wind | distance) "
                 "across trips\n(score=max, stratum=all sites)", fontsize=11)
    fig.colorbar(im, ax=axes[-1], label="r_part", shrink=0.8)
    fig.savefig(FIG / "fig_wind_sensitivity.pdf")
    plt.close(fig)

    # --- Figure 2: score-type comparison (window x score_type), per compartment
    fig, axes = plt.subplots(1, 3, figsize=(13, 4), constrained_layout=True)
    for ax, comp in zip(axes, ["rhizosphere", "surface", "deep"]):
        s = df[(df["dust_thr"] == 7.0) & (df["angle_tol"] == 30) &
               (df["stratum"] == "all") & (df["compartment"] == comp)]
        s_med = s.groupby(["window_days", "score_type"])["r_part"].median().reset_index()
        im = heatmap(ax, s_med, "window_days", "score_type", "r_part", comp)
    fig.suptitle("Score-type comparison: median r_part (thr=7 m/s, angle=30°, all sites)")
    fig.colorbar(im, ax=axes[-1], label="r_part", shrink=0.8)
    fig.savefig(FIG / "fig_wind_score_type.pdf")
    plt.close(fig)

    # --- Figure 3: distance stratification (window x stratum), per compartment
    fig, axes = plt.subplots(1, 3, figsize=(13, 4), constrained_layout=True)
    for ax, comp in zip(axes, ["rhizosphere", "surface", "deep"]):
        s = df[(df["dust_thr"] == 7.0) & (df["angle_tol"] == 30) &
               (df["score_type"] == "max") & (df["compartment"] == comp)]
        s_med = s.groupby(["window_days", "stratum"])["r_part"].median().reset_index()
        # Re-order strata
        s_med["stratum"] = pd.Categorical(s_med["stratum"],
            categories=["all", "lt100km", "100_500km", "gt500km"], ordered=True)
        s_med = s_med.sort_values(["window_days", "stratum"])
        im = heatmap(ax, s_med, "window_days", "stratum", "r_part", comp)
    fig.suptitle("Distance stratification: median r_part (score=max, thr=7, ang=30°)")
    fig.colorbar(im, ax=axes[-1], label="r_part", shrink=0.8)
    fig.savefig(FIG / "fig_wind_distance_strata.pdf")
    plt.close(fig)

    print(f"Wrote {FIG}/fig_wind_sensitivity.pdf, fig_wind_score_type.pdf, "
          f"fig_wind_distance_strata.pdf")


if __name__ == "__main__":
    main()
