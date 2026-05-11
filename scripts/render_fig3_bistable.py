#!/usr/bin/env python3
"""Main Fig 3 (Bistable architecture) — A vs B strategies.

Panels:
  (a) Bimodal log2(A/B) distribution across all samples.
  (b) Per-(site, comp, trip) dominance heatmap.
  (c) Stratified networks A-dominant vs B-dominant (summary stats).
  (d) Strategy A vs B functional fingerprint (Nibribacter MAGs corrected KEGG).

Reads:
  cache/two_strategy_temporal/per_sample_strategy_with_precip.tsv
  cache/network_A_vs_B/network_summary.tsv
  cache/network_A_vs_B/keystone_{A,B}_dominant.tsv
  cache/nibribacter_mags/corrected_function_summary.tsv
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

C_A = "#117733"  # Strategy A (DOM-cyclers, green)
C_B = "#882255"  # Strategy B (halotolerant, magenta)
C_NEUTRAL = "#888888"


def panel_a_bimodal(ax):
    df = pd.read_csv(CACHE / "two_strategy_temporal" /
                       "per_sample_strategy_with_precip.tsv", sep="\t")
    s = df["log2_A_over_B"].replace([np.inf, -np.inf], np.nan).dropna()
    # Cap at ±20 for display
    s_disp = s.clip(-20, 20)
    ax.hist(s_disp[s_disp > 0], bins=40, color=C_A, alpha=0.7,
              edgecolor="white", linewidth=0.3, label="A-dominant")
    ax.hist(s_disp[s_disp < 0], bins=40, color=C_B, alpha=0.7,
              edgecolor="white", linewidth=0.3, label="B-dominant")
    ax.axvline(0, color="grey", linestyle="--", linewidth=0.7, alpha=0.6)
    n_A = (df["dominant"] == "A").sum()
    n_B = (df["dominant"] == "B").sum()
    ax.text(0.97, 0.95, f"A-dominant: {n_A}\nB-dominant: {n_B}",
              transform=ax.transAxes, ha="right", va="top", fontsize=7,
              bbox=dict(facecolor="white", edgecolor="none", alpha=0.85))
    ax.set_xlabel("log₂(Σ Strategy A / Σ Strategy B)")
    ax.set_ylabel("Samples")
    ax.set_title("a  Bimodal dominance distribution",
                    loc="left", fontweight="bold")
    ax.legend(frameon=False, fontsize=6, loc="upper left")


def panel_b_heatmap(ax):
    df = pd.read_csv(CACHE / "two_strategy_temporal" /
                       "per_sample_strategy_with_precip.tsv", sep="\t")
    # Aggregate per (site, comp, trip)
    df["site_comp"] = df["site"].astype(str) + "/" + df["compartment"]
    agg = (df.groupby(["site", "compartment", "trip"])["log2_A_over_B"]
              .median().reset_index())
    # Pivot to (site_comp) × trip
    agg["site_comp"] = agg["site"].astype(str) + "/" + agg["compartment"]
    p = agg.pivot_table(index="site_comp", columns="trip",
                              values="log2_A_over_B", aggfunc="median")
    # Sort by mean dominance
    p["_m"] = p.mean(axis=1)
    p = p.sort_values("_m", ascending=False).drop(columns="_m")
    p_disp = p.clip(-10, 10)
    im = ax.imshow(p_disp.values, aspect="auto", cmap="RdYlGn",
                       vmin=-10, vmax=10, interpolation="nearest")
    ax.set_xticks(np.arange(p.shape[1]))
    ax.set_xticklabels([f"T{int(t)}" for t in p.columns], rotation=0)
    ax.set_yticks([])
    ax.set_xlabel("Trip")
    ax.set_ylabel(f"(Site, compartment) cells (n = {p.shape[0]})")
    ax.set_title("b  Per-cell dominance across trips",
                    loc="left", fontweight="bold")
    cax = ax.inset_axes([1.02, 0.0, 0.04, 1.0])
    cbar = plt.colorbar(im, cax=cax)
    cbar.set_label("log₂(A / B)", fontsize=7)
    cbar.ax.tick_params(labelsize=6)


def panel_c_networks(ax):
    """Network summary: nodes/edges/modules + top keystones."""
    summ = pd.read_csv(CACHE / "network_A_vs_B" / "network_summary.tsv",
                         sep="\t")
    ksA = pd.read_csv(CACHE / "network_A_vs_B" / "keystone_A_dominant.tsv",
                         sep="\t").head(6)
    ksB = pd.read_csv(CACHE / "network_A_vs_B" / "keystone_B_dominant.tsv",
                         sep="\t").head(6)
    ax.axis("off")
    ax.set_title("c  Stratified network architecture",
                    loc="left", fontweight="bold")
    # Two columns: A-dom and B-dom
    fig = ax.figure
    bbox = ax.get_position()
    w = bbox.width / 2 - 0.005
    h = bbox.height
    ax_left = fig.add_axes([bbox.x0, bbox.y0, w, h])
    ax_right = fig.add_axes([bbox.x0 + w + 0.01, bbox.y0, w, h])
    for axi, klass, kdf, color, header_color in [
            (ax_left, "A", ksA, C_A, "#117733"),
            (ax_right, "B", ksB, C_B, "#882255")]:
        axi.axis("off")
        sub = summ[summ["dominance_class"] == klass].iloc[0]
        n_nodes = int(sub["n_genera"])
        n_edges = int(sub["n_edges"])
        n_mod = int(sub["n_modules"])
        n_samp = int(sub["n_samples"])
        axi.text(0.5, 0.95,
                    f"{klass}-dominant network",
                    transform=axi.transAxes, ha="center", va="top",
                    fontsize=8.5, fontweight="bold", color=header_color)
        axi.text(0.5, 0.85,
                    f"n = {n_samp:,} samples",
                    transform=axi.transAxes, ha="center", va="top",
                    fontsize=7, color="#444444")
        axi.text(0.5, 0.77,
                    f"{n_nodes} genera · {n_edges} edges · {n_mod} modules",
                    transform=axi.transAxes, ha="center", va="top",
                    fontsize=7, color="#444444")
        axi.text(0.5, 0.68, "Top keystones:",
                    transform=axi.transAxes, ha="center", va="top",
                    fontsize=7, fontweight="bold")
        for i, (_, row) in enumerate(kdf.iterrows()):
            axi.text(0.5, 0.60 - i * 0.08,
                        f"{i+1}. {row['node']}",
                        transform=axi.transAxes, ha="center", va="top",
                        fontsize=7,
                        color=header_color if i == 0 else "#333333")
    # Annotation box: density comparison
    summ_A = summ[summ["dominance_class"] == "A"].iloc[0]
    summ_B = summ[summ["dominance_class"] == "B"].iloc[0]
    density_ratio = (summ_B["n_edges"] / summ_B["n_genera"]) / \
                       (summ_A["n_edges"] / summ_A["n_genera"])
    fig.text(bbox.x0 + bbox.width / 2, bbox.y0 - 0.04,
                f"B-network is {density_ratio:.1f}× denser, fewer modules",
                ha="center", va="top", fontsize=7, fontstyle="italic",
                color="#555555")


def panel_d_function(ax):
    """Functional fingerprint: Strategy A categories from Nibribacter
    corrected KEGG analysis."""
    d = pd.read_csv(CACHE / "nibribacter_mags" /
                       "corrected_function_summary.tsv", sep="\t")
    # Aggregate across MAGs → unique KOs sum
    agg = (d.groupby("category")
              .agg(n_kos_in_list=("n_kos_in_list", "first"),
                    n_unique_kos_found_sum=("n_unique_kos_found", "sum"))
              .reset_index())
    agg["pct_of_target"] = (
        agg["n_unique_kos_found_sum"] / agg["n_kos_in_list"] * 100)
    # Highlight A-signature and B-signature categories
    A_signature = ["DNA_repair_specific", "Heat_shock_specific",
                       "Trehalose_biosynth_strict",
                       "OxStress_specific",
                       "GH_glycoside_hydrolase", "TonB_SusCD"]
    B_signature = ["Ectoine_biosynth_strict", "Betaine_uptake_specific",
                       "Sporulation_initiation_negctrl"]
    NEG_CTRL = ["Photosynthesis_negctrl", "Sporulation_initiation_negctrl"]
    # Display order: A-sig then B-sig (which should mostly be absent in
    # Nibribacter MAGs)
    keep = A_signature + ["Betaine_biosynth_strict",
                              "Ectoine_biosynth_strict",
                              "Sporulation_initiation_negctrl",
                              "Photosynthesis_negctrl"]
    sub = agg[agg["category"].isin(keep)].copy()
    sub["order"] = sub["category"].apply(lambda c: keep.index(c) if c in keep else 99)
    sub = sub.sort_values("order")
    # Color: A-sig green; B-sig magenta; neg-ctrl grey
    def colour(c):
        if c in A_signature: return C_A
        if c in NEG_CTRL: return C_NEUTRAL
        return C_B
    ax.barh(np.arange(len(sub)), sub["pct_of_target"].values,
              color=[colour(c) for c in sub["category"]],
              edgecolor="white", linewidth=0.4)
    labels = [c.replace("_strict", "").replace("_specific", "")
                  .replace("_negctrl", " (ctrl)").replace("_", " ")
                  for c in sub["category"]]
    ax.set_yticks(np.arange(len(sub)))
    ax.set_yticklabels(labels, fontsize=7)
    ax.invert_yaxis()
    ax.axvline(100, color="grey", linestyle="--", linewidth=0.6, alpha=0.6)
    ax.set_xlabel("% of target KO list detected")
    ax.set_title("d  Strategy A MAGs: functional fingerprint\n"
                    "(green = A signature; magenta = B signature; grey = control)",
                    loc="left", fontweight="bold", fontsize=8.5)
    ax.grid(axis="x", alpha=0.3, linewidth=0.4)


def main():
    fig = plt.figure(figsize=(7.2, 7.5))
    gs = GridSpec(3, 2, figure=fig, wspace=0.45, hspace=0.55,
                     left=0.10, right=0.94, top=0.95, bottom=0.07,
                     height_ratios=[1, 1.2, 1.5])
    ax_a = fig.add_subplot(gs[0, 0])
    ax_b = fig.add_subplot(gs[0, 1])
    ax_c = fig.add_subplot(gs[1, :])
    ax_d = fig.add_subplot(gs[2, :])
    panel_a_bimodal(ax_a)
    panel_b_heatmap(ax_b)
    panel_c_networks(ax_c)
    panel_d_function(ax_d)
    out = FIG / "fig3_bistable_v2.pdf"
    fig.savefig(out, dpi=300, bbox_inches="tight")
    print(f"wrote {out}", flush=True)


if __name__ == "__main__":
    main()
