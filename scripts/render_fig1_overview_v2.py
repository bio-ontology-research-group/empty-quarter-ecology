#!/usr/bin/env python3
"""Main Fig 1 (Overview) — bistable-narrative version.

Panels:
  (a) 60-site map across the Rub' al-Khali transect.
  (b) Per-(trip, compartment) sample count.
  (c) Shannon diversity by compartment (boxplot).
  (d) XRF sabkha-vs-sandy axis (PC1 of major elements; or Si/Ca ratio).

Reads:
  data/geodata/trip*_geodata.tsv
  cache/metadata.parquet
  cache/per_sample_shannon.tsv
  cache/xrf_summary_all_trips.tsv  (per-compartment summaries)
  cache/xrf_lithology_pca.tsv      (if present)
"""
from __future__ import annotations

from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec

mpl.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
    "font.size": 8,
    "axes.titlesize": 9,
    "axes.labelsize": 8,
    "xtick.labelsize": 7,
    "ytick.labelsize": 7,
    "legend.fontsize": 7,
    "axes.linewidth": 0.7,
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
})

REPO = Path(__file__).resolve().parents[1]
CACHE = REPO / "cache"
FIG = REPO.parent / "figures"
FIG.mkdir(exist_ok=True)

COMP_COLORS = {"surface": "#ee7733", "deep": "#0077bb", "rhizosphere": "#33aa55"}
COMP_ORDER = ["surface", "deep", "rhizosphere"]


def panel_a_map(ax):
    sites = []
    for t in (1,):  # T1 has full 60 sites
        g = pd.read_csv(REPO / "data" / "geodata" / f"trip{t}_geodata.tsv",
                          sep="\t")
        g["Site"] = pd.to_numeric(g["Site"], errors="coerce")
        g = g.dropna(subset=["Site"])
        sites.append(g[["Site", "Latitude", "Longitude"]])
    s = pd.concat(sites, ignore_index=True).drop_duplicates("Site")
    s = s[s["Site"] <= 60]
    ax.scatter(s["Longitude"], s["Latitude"], s=22, c="#cc3311",
                  edgecolors="black", linewidths=0.4, zorder=3)
    # rough Saudi/Empty Quarter outline
    ax.set_xlim(44, 56)
    ax.set_ylim(18, 22)
    ax.set_xlabel("Longitude (°E)")
    ax.set_ylabel("Latitude (°N)")
    ax.set_title("a  60 sites · Rub' al-Khali transect", loc="left",
                    fontweight="bold")
    ax.set_aspect("equal", adjustable="box")
    ax.grid(alpha=0.3, linewidth=0.4)
    ax.text(0.97, 0.05, f"n = {len(s)} sites", transform=ax.transAxes,
              ha="right", fontsize=7,
              bbox=dict(facecolor="white", edgecolor="none", alpha=0.8))


def panel_b_design(ax):
    md = pd.read_parquet(CACHE / "metadata.parquet")
    md = md[md["control"] != True]
    cnt = (md.groupby(["trip", "compartment"]).size()
              .unstack(fill_value=0)[COMP_ORDER])
    cnt.plot(kind="bar", stacked=True,
              color=[COMP_COLORS[c] for c in COMP_ORDER],
              ax=ax, width=0.7, edgecolor="white", linewidth=0.4)
    ax.set_xlabel("Sampling trip")
    ax.set_ylabel("Samples (n)")
    ax.set_title("b  Sampling design", loc="left", fontweight="bold")
    ax.set_xticklabels([f"T{i}" for i in cnt.index], rotation=0)
    ax.legend(title="Compartment", frameon=False, fontsize=6,
                title_fontsize=7, loc="upper left")
    ax.grid(axis="y", alpha=0.3, linewidth=0.4)
    total = cnt.sum().sum()
    ax.text(0.97, 0.95, f"N = {total:,} samples", transform=ax.transAxes,
              ha="right", va="top", fontsize=7,
              bbox=dict(facecolor="white", edgecolor="none", alpha=0.8))


def panel_c_shannon(ax):
    sh = pd.read_csv(CACHE / "per_sample_shannon.tsv", sep="\t")
    md = pd.read_parquet(CACHE / "metadata.parquet").reset_index()
    df = sh.merge(md[["sample", "compartment", "control"]], on="sample")
    df = df[df["control"] != True]
    data = [df[df["compartment"] == c]["shannon"].dropna().values
              for c in COMP_ORDER]
    bp = ax.boxplot(data, labels=COMP_ORDER, patch_artist=True,
                       widths=0.5, showfliers=False,
                       medianprops=dict(color="black", linewidth=0.9))
    for patch, c in zip(bp["boxes"], COMP_ORDER):
        patch.set_facecolor(COMP_COLORS[c])
        patch.set_alpha(0.7)
        patch.set_edgecolor("black")
        patch.set_linewidth(0.5)
    ax.set_ylabel("Shannon α-diversity (ASV)")
    ax.set_title("c  α-diversity by compartment", loc="left",
                    fontweight="bold")
    ax.grid(axis="y", alpha=0.3, linewidth=0.4)


def panel_d_xrf(ax):
    # Si:Ca proxy for sabkha-vs-sandy
    try:
        xrf = pd.read_csv(REPO / "data" / "geochemistry" /
                            "xrf_lab_table_all_trips.tsv", sep="\t")
        xrf["compartment"] = xrf["compartment"].str.lower()
        # Approximate sabkha index = sum of (SO3, Na, Cl, Ca, Mg) / sum of (Si, Al, Fe)
        sab_cols = [c for c in ["SO3", "Na", "Cl", "Ca", "Mg"]
                       if c in xrf.columns]
        san_cols = [c for c in ["Si", "Al", "Fe"] if c in xrf.columns]
        if not sab_cols or not san_cols:
            ax.text(0.5, 0.5, "XRF columns not available",
                       ha="center", va="center", transform=ax.transAxes)
            ax.set_title("d  Substrate axis (sabkha · sandy)",
                            loc="left", fontweight="bold")
            return
        xrf["sabkha_index"] = (np.log10(xrf[sab_cols].sum(axis=1) + 1) -
                                   np.log10(xrf[san_cols].sum(axis=1) + 1))
        xrf = xrf.dropna(subset=["sabkha_index", "compartment"])
        for c in COMP_ORDER:
            sub = xrf[xrf["compartment"] == c]
            if len(sub) == 0:
                continue
            ax.hist(sub["sabkha_index"], bins=25, alpha=0.55,
                      color=COMP_COLORS[c], label=c, edgecolor="white",
                      linewidth=0.3)
        ax.axvline(0, color="grey", linestyle="--", linewidth=0.6,
                      alpha=0.5)
        ax.set_xlabel("Sabkha index (log evaporite / silicate)")
        ax.set_ylabel("Samples")
        ax.set_title("d  Substrate axis (sandy ↔ sabkha)",
                        loc="left", fontweight="bold")
        ax.legend(frameon=False, fontsize=6, loc="upper left")
        ax.grid(axis="y", alpha=0.3, linewidth=0.4)
    except Exception as e:
        ax.text(0.5, 0.5, f"XRF panel: {e}", ha="center", va="center",
                   transform=ax.transAxes, fontsize=6)
        ax.set_title("d  Substrate axis", loc="left", fontweight="bold")


def main():
    fig = plt.figure(figsize=(7.2, 5.4))
    gs = GridSpec(2, 2, figure=fig, wspace=0.32, hspace=0.42,
                     left=0.08, right=0.97, top=0.95, bottom=0.10)
    ax_a = fig.add_subplot(gs[0, 0])
    ax_b = fig.add_subplot(gs[0, 1])
    ax_c = fig.add_subplot(gs[1, 0])
    ax_d = fig.add_subplot(gs[1, 1])
    panel_a_map(ax_a)
    panel_b_design(ax_b)
    panel_c_shannon(ax_c)
    panel_d_xrf(ax_d)
    out = FIG / "fig1_overview_v2.pdf"
    fig.savefig(out, dpi=300, bbox_inches="tight")
    print(f"wrote {out}", flush=True)


if __name__ == "__main__":
    main()
