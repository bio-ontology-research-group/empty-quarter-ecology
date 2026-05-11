#!/usr/bin/env python3
"""Main Fig 2 (Dispersal-driven + relic-DNA) — bistable-narrative version.

Panels:
  (a) iCAMP process attribution per compartment (stacked bar; the 67% point).
  (b) Wind partial Mantel r per (trip, compartment).
  (c) Relic-likelihood-score distribution with alive cutoff (0.3).
  (d) Shannon vs MAT — all / alive / relic stratified (THE relic artifact).

Reads:
  cache/icamp/process_summary_all.tsv
  cache/mantel_per_trip_compartment.tsv
  cache/relic_priors/relic_score_with_mag_prior.tsv
  cache/per_sample_shannon.tsv  + cache/feature_table_alive,relic.parquet
  cache/metadata.parquet  + data/geodata/trip1_geodata.tsv
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
})

REPO = Path(__file__).resolve().parents[1]
CACHE = REPO / "cache"
FIG = REPO.parent / "figures"

COMP_COLORS = {"surface": "#ee7733", "deep": "#0077bb", "rhizosphere": "#33aa55"}
COMP_ORDER = ["surface", "deep", "rhizosphere"]
PROCESS_COLORS = {
    "homogeneous_selection": "#882255",
    "variable_selection":    "#cc6677",
    "homogenizing_dispersal":"#117733",   # green = the headline 67%
    "dispersal_limitation":  "#88ccee",
    "drift_or_weak":         "#aaaaaa",
}
PROCESS_ORDER = ["homogeneous_selection", "variable_selection",
                    "homogenizing_dispersal", "dispersal_limitation",
                    "drift_or_weak"]


def panel_a_icamp(ax):
    d = pd.read_csv(CACHE / "icamp" / "process_summary_all.tsv", sep="\t")
    d = d.set_index("compartment").loc[COMP_ORDER]
    fracs = d[[f"frac_{p}" for p in PROCESS_ORDER]].values  # rows = comp
    bottom = np.zeros(len(COMP_ORDER))
    for i, p in enumerate(PROCESS_ORDER):
        ax.bar(np.arange(len(COMP_ORDER)), fracs[:, i] * 100,
                  bottom=bottom * 100, color=PROCESS_COLORS[p],
                  edgecolor="white", linewidth=0.5,
                  label=p.replace("_", " "))
        bottom = bottom + fracs[:, i]
    # Highlight the 67% homogenizing dispersal band
    for i, c in enumerate(COMP_ORDER):
        h = d.loc[c, "frac_homogenizing_dispersal"]
        ax.text(i, h * 100 * 0.5 + (1 - h) * 0  # midpoint approx
                  + (1 - d.loc[c, "frac_homogenizing_dispersal"] -
                      d.loc[c, "frac_dispersal_limitation"] -
                      d.loc[c, "frac_drift_or_weak"]) * 100,
                  f"{h*100:.0f}%",
                  ha="center", va="center", fontsize=9, fontweight="bold",
                  color="white")
    ax.set_xticks(np.arange(len(COMP_ORDER)))
    ax.set_xticklabels(COMP_ORDER)
    ax.set_ylabel("% of pairwise comparisons")
    ax.set_ylim(0, 100)
    ax.set_title("a  iCAMP assembly processes", loc="left",
                    fontweight="bold")
    ax.legend(frameon=False, fontsize=5.5, loc="upper right",
                bbox_to_anchor=(1.02, 1.02), ncol=1)
    ax.grid(axis="y", alpha=0.3, linewidth=0.4)


def panel_b_wind(ax):
    m = pd.read_csv(CACHE / "mantel_per_trip_compartment.tsv", sep="\t")
    # Plot wind partial-Mantel r per compartment, jittered by trip
    for c in COMP_ORDER:
        sub = m[m["compartment"] == c]
        if len(sub) == 0:
            continue
        x = np.array(sub["trip"]) + np.random.RandomState(42).uniform(
            -0.12, 0.12, size=len(sub))
        ax.scatter(x, sub["mantel_r"], s=30, c=COMP_COLORS[c],
                      edgecolors="black", linewidths=0.4, label=c,
                      zorder=3, alpha=0.85)
    ax.axhline(0, color="grey", linestyle="--", linewidth=0.6, alpha=0.6)
    ax.set_xlabel("Sampling trip")
    ax.set_ylabel("Partial Mantel r\n(BC ~ wind | distance)")
    ax.set_title("b  Wind-dispersal partial Mantel",
                    loc="left", fontweight="bold")
    ax.legend(frameon=False, fontsize=6, loc="upper left")
    ax.grid(axis="y", alpha=0.3, linewidth=0.4)
    ax.text(0.97, 0.97, "annual-scale\nwind memory",
              transform=ax.transAxes, ha="right", va="top", fontsize=6,
              color="#555555",
              bbox=dict(facecolor="white", edgecolor="#cccccc",
                          alpha=0.85, linewidth=0.4))


def panel_c_relic(ax):
    rs = pd.read_csv(CACHE / "relic_priors" / "relic_score_with_mag_prior.tsv",
                       sep="\t")
    s = rs["relic_score_with_mag"].dropna()
    ax.hist(s, bins=40, color="#bbbbbb", edgecolor="white", linewidth=0.3)
    ax.axvline(0.3, color="#cc3311", linestyle="--", linewidth=1.0,
                  label="alive cutoff 0.3")
    frac_relic = (s > 0.3).mean()
    ax.text(0.97, 0.85,
              f"{frac_relic*100:.0f}% of ASVs\nclassified relic",
              transform=ax.transAxes, ha="right", va="top", fontsize=7,
              color="#cc3311",
              bbox=dict(facecolor="white", edgecolor="#cc3311",
                          linewidth=0.5))
    ax.set_xlabel("Relic-likelihood score (with MAG prior)")
    ax.set_ylabel("ASVs")
    ax.set_title("c  Composite relic indicator (n = {:,} ASVs)".format(len(s)),
                    loc="left", fontweight="bold")
    ax.legend(frameon=False, fontsize=6, loc="upper left")


def _shannon_panel(df, col, ax, color, title, rho_label_xy=(0.05, 0.95)):
    sub = df.dropna(subset=[col, "AnnualMeanTemp"])
    if len(sub) == 0:
        ax.text(0.5, 0.5, "n=0", ha="center", va="center",
                  transform=ax.transAxes)
        return
    rho = sub[[col, "AnnualMeanTemp"]].corr(method="spearman").iloc[0, 1]
    ax.scatter(sub["AnnualMeanTemp"], sub[col], s=4, alpha=0.35,
                  c=color, edgecolors="none")
    # rolling mean line
    x = sub["AnnualMeanTemp"].values
    y = sub[col].values
    idx = np.argsort(x)
    xs = x[idx]; ys = y[idx]
    win = max(20, len(sub) // 25)
    if len(xs) > win:
        rm_y = pd.Series(ys).rolling(win, center=True).mean().values
        ax.plot(xs, rm_y, color=color, linewidth=1.2, alpha=0.9)
    ax.text(rho_label_xy[0], rho_label_xy[1],
              f"ρ = {rho:+.2f}", transform=ax.transAxes, va="top", ha="left",
              fontsize=8, fontweight="bold", color=color,
              bbox=dict(facecolor="white", edgecolor=color, linewidth=0.5))
    ax.set_title(title, fontsize=8, loc="left")
    ax.grid(alpha=0.3, linewidth=0.4)


def panel_d_shannon_mat(parent_ax):
    sh = pd.read_csv(CACHE / "per_sample_shannon.tsv", sep="\t")
    md = pd.read_parquet(CACHE / "metadata.parquet").reset_index()
    geo = pd.read_csv(REPO / "data" / "geodata" / "trip1_geodata.tsv",
                          sep="\t")
    geo["site"] = pd.to_numeric(geo["Site"], errors="coerce")
    geo = geo.dropna(subset=["site"])
    geo["site"] = geo["site"].astype(int)
    df = sh.merge(md[["sample", "site", "compartment", "control"]],
                     on="sample")
    df = df[df["control"] != True]
    df["site"] = df["site"].astype(int)
    df = df.merge(geo[["site", "AnnualMeanTemp"]], on="site", how="left")
    # alive + relic Shannon
    for label, ft_path in [("shannon_alive", "feature_table_alive.parquet"),
                                ("shannon_relic", "feature_table_relic.parquet")]:
        ft = pd.read_parquet(CACHE / ft_path)
        ra = ft.div(ft.sum(axis=0).replace(0, 1), axis=1)
        sh_x = (-(ra * np.log(ra.where(ra > 0, 1))).sum(axis=0))
        sh_x = sh_x.rename(label).reset_index().rename(
            columns={"index": "sample"})
        df = df.merge(sh_x, on="sample", how="left")
    # Build 3 sub-axes inside parent_ax position
    bbox = parent_ax.get_position()
    parent_ax.axis("off")
    parent_ax.set_title("d  Shannon ~ MAT is a relic artefact",
                            loc="left", fontweight="bold")
    n = 3
    w = (bbox.width - 0.02 * (n - 1)) / n
    fig = parent_ax.figure
    axes = [fig.add_axes([bbox.x0 + i * (w + 0.02), bbox.y0, w, bbox.height])
              for i in range(n)]
    _shannon_panel(df, "shannon", axes[0], "#222222",
                    "All ASVs")
    _shannon_panel(df, "shannon_alive", axes[1], "#117733",
                    "Alive only (score ≤ 0.3)")
    _shannon_panel(df, "shannon_relic", axes[2], "#cc3311",
                    "Relic only (score > 0.3)")
    for ax in axes:
        ax.set_xlabel("MAT (°C)")
    axes[0].set_ylabel("Shannon")


def main():
    fig = plt.figure(figsize=(7.2, 6.0))
    gs = GridSpec(2, 2, figure=fig, wspace=0.30, hspace=0.50,
                     left=0.08, right=0.97, top=0.95, bottom=0.10)
    ax_a = fig.add_subplot(gs[0, 0])
    ax_b = fig.add_subplot(gs[0, 1])
    ax_c = fig.add_subplot(gs[1, 0])
    ax_d = fig.add_subplot(gs[1, 1])
    panel_a_icamp(ax_a)
    panel_b_wind(ax_b)
    panel_c_relic(ax_c)
    panel_d_shannon_mat(ax_d)
    out = FIG / "fig2_dispersal_relic_v2.pdf"
    fig.savefig(out, dpi=300, bbox_inches="tight")
    print(f"wrote {out}", flush=True)


if __name__ == "__main__":
    main()
