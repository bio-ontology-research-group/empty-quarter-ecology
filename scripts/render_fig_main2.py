#!/usr/bin/env python3
"""Render Main Figure 2 (Trans-biome generalisation) — Nature-quality.

Panels:
  (a) Cross-desert prevalence of CSP1-2 at 85% V4 identity (5 deserts).
  (b) CSP1-2 relative abundance vs Shannon — Empty Quarter + Namib.
  (c) Salinity proxy vs Shannon — Empty Quarter (XRF S) + Atacama (EC).
  (d) Compositional co-occurrence ego networks centred on CSP1-2.

Reads:
  cache/crossdesert/per_sample.tsv         (panels a/b/c)
  cache/crossdesert/comparison_summary.tsv (panel a)
  cache/causal_frame_tier1.parquet         (panels b/c — Empty Quarter)
  cache/network_edges_*.tsv                (panel d)
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import networkx as nx

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
    "xtick.major.width": 0.7,
    "ytick.major.width": 0.7,
    "lines.linewidth": 1.0,
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
})

REPO = Path(__file__).resolve().parents[1]
CACHE = REPO / "cache"
FIG = REPO.parent / "figures"

# Colour palette per desert (colorblind-safe)
DESERT_COLORS = {
    "Namib":         "#117733",   # green
    "Gurbantunggut": "#88ccee",   # sky blue
    "EmptyQuarter":  "#000000",   # black — focal
    "McMurdo":       "#7e57c2",   # purple
    "Atacama":       "#cc3311",   # red
}
DESERT_LABEL = {
    "Namib": "Namib",
    "Gurbantunggut": "Gurbantunggut",
    "EmptyQuarter": "Empty Quarter",
    "McMurdo": "McMurdo (cold)",
    "Atacama": "Atacama",
}
C_KEYSTONE = "#cc3311"
C_NEUTRAL = "#bbbbbb"

# ----------------------------------------------------------------------------
# Load
# ----------------------------------------------------------------------------
xs = pd.read_csv(CACHE / "crossdesert" / "per_sample.tsv", sep="\t")
summ = pd.read_csv(CACHE / "crossdesert" / "comparison_summary.tsv", sep="\t")

# Empty Quarter local frame for panels b and c
frame = pd.read_parquet(CACHE / "causal_frame_tier1.parquet")

# ----------------------------------------------------------------------------
# Figure
# ----------------------------------------------------------------------------
fig = plt.figure(figsize=(7.2, 9.2))
gs = GridSpec(
    3, 6, figure=fig,
    width_ratios=[1, 1, 1, 1, 1, 1],
    height_ratios=[1.0, 1.0, 0.9],
    hspace=0.65, wspace=1.6,
    left=0.09, right=0.97, top=0.95, bottom=0.05,
)

# ============================================================================
# (a) Cross-desert prevalence
# ============================================================================
ax_a = fig.add_subplot(gs[0, :3])
desert_order = ["Namib", "Gurbantunggut", "EmptyQuarter", "McMurdo", "Atacama"]

# Gurbantunggut isn't in the summary table — read its per-sample TSV.
gurb = pd.read_csv(CACHE / "crossdesert" / "gurbantunggut_per_sample.tsv",
                   sep="\t")
gurb_n = len(gurb)
gurb_prev = (gurb["csp_rel_85"] > 0).mean()

prev = []
for d in desert_order:
    if d == "Gurbantunggut":
        prev.append((d, float(gurb_prev), gurb_n))
        continue
    row = summ[summ["desert"].str.replace(" (baseline)", "", regex=False) == d]
    if len(row):
        f = row.iloc[0]["frac_samples_with_CSP85"]
        n = row.iloc[0]["n_samples"]
        prev.append((d, float(f), int(n)))
    else:
        row = summ[summ["desert"].str.contains(d, case=False, regex=False)]
        if len(row):
            f = row.iloc[0]["frac_samples_with_CSP85"]
            n = row.iloc[0]["n_samples"]
            prev.append((d, float(f), int(n)))

x = np.arange(len(prev))
heights = [p[1] for p in prev]
labels = [DESERT_LABEL[p[0]] for p in prev]
ns = [p[2] for p in prev]
colors = [DESERT_COLORS[p[0]] for p in prev]

bars = ax_a.bar(x, heights, color=colors, edgecolor="black",
                linewidth=0.5, width=0.7)
for xi, h, n in zip(x, heights, ns):
    ax_a.text(xi, h + 0.025, f"{h*100:.0f}%",
              ha="center", va="bottom", fontsize=7.5, fontweight="bold")

ax_a.set_xticks(x)
xtick_labels = [f"{labels[i]}\n(n={ns[i]})" for i in range(len(labels))]
ax_a.set_xticklabels(xtick_labels, fontsize=6.5, rotation=25, ha="right",
                     rotation_mode="anchor")
ax_a.set_ylabel("Fraction of samples with CSP1-2\n(85% V4 identity)", fontsize=7.5)
ax_a.set_ylim(0, 1.05)
ax_a.set_yticks([0, 0.25, 0.5, 0.75, 1.0])
ax_a.set_yticklabels(["0%", "25%", "50%", "75%", "100%"])
ax_a.tick_params(direction="out", length=2.5)
for sp in ("top", "right"):
    ax_a.spines[sp].set_visible(False)
ax_a.grid(True, axis="y", alpha=0.25, lw=0.4)
ax_a.set_axisbelow(True)

ax_a.text(-0.13, 1.08, "a", transform=ax_a.transAxes,
          fontsize=12, fontweight="bold", va="top", ha="left")
ax_a.set_title("CSP1-2 prevalence across global drylands",
               loc="left", pad=4, fontsize=8.5)

# ============================================================================
# (b) Abundance vs Shannon — EQ + Namib
# ============================================================================
ax_b = fig.add_subplot(gs[0, 3:])
# Empty Quarter: from local frame (csp_relab vs shannon)
eq = frame[(frame["csp_relab"] > 0)
           & (frame["shannon"].notna())].copy()
eq["x"] = eq["csp_relab"] * 100  # %
ax_b.scatter(eq["x"], eq["shannon"], s=10, alpha=0.45,
             c=DESERT_COLORS["EmptyQuarter"],
             edgecolors="none", zorder=2,
             label=f"Empty Quarter (n={len(eq)})")

# Namib subset
namib = xs[(xs["desert"] == "Namib") & (xs["csp_rel_85"] > 0)].copy()
namib["x"] = namib["csp_rel_85"] * 100
ax_b.scatter(namib["x"], namib["shannon_otu97"], s=18, alpha=0.85,
             c=DESERT_COLORS["Namib"],
             edgecolors="black", linewidths=0.3, zorder=3,
             label=f"Namib (n={len(namib)})")

# Log-x regression lines
def loglinfit(x, y):
    mask = (x > 0) & np.isfinite(x) & np.isfinite(y)
    lx = np.log10(x[mask])
    yv = y[mask].astype(float)
    if len(lx) < 5:
        return None
    p = np.polyfit(lx, yv, 1)
    xx = np.linspace(lx.min(), lx.max(), 100)
    yy = np.polyval(p, xx)
    rho = np.corrcoef(lx, yv)[0, 1]
    return 10 ** xx, yy, rho

for grp_x, grp_y, c, lbl in [
    (eq["x"].values, eq["shannon"].values, DESERT_COLORS["EmptyQuarter"], "EQ"),
    (namib["x"].values, namib["shannon_otu97"].values,
     DESERT_COLORS["Namib"], "Namib"),
]:
    res = loglinfit(grp_x, grp_y)
    if res:
        xx, yy, rho = res
        ax_b.plot(xx, yy, color=c, lw=1.4, ls="--", alpha=0.85, zorder=4)
        # Stack ρ annotations top-left to avoid data overlap
        x_pos = 0.04
        y_pos = 0.96 if lbl == "EQ" else 0.88
        ax_b.text(x_pos, y_pos, f"{lbl} ρ = {rho:+.2f}",
                  transform=ax_b.transAxes, ha="left", va="top",
                  fontsize=7, color=c, fontweight="bold")

ax_b.set_xscale("log")
ax_b.set_xlabel("CSP1-2 relative abundance (%)", fontsize=7.5)
ax_b.set_ylabel("Shannon diversity", fontsize=7.5)
ax_b.tick_params(direction="out", length=2.5)
for sp in ("top", "right"):
    ax_b.spines[sp].set_visible(False)
ax_b.legend(loc="lower right", frameon=False, fontsize=6.5,
            handletextpad=0.4, borderpad=0.3, scatterpoints=1)
ax_b.grid(True, alpha=0.2, lw=0.4)
ax_b.set_axisbelow(True)

ax_b.text(-0.18, 1.08, "b", transform=ax_b.transAxes,
          fontsize=12, fontweight="bold", va="top", ha="left")
ax_b.set_title("Diversity tracks CSP1-2 abundance — two deserts",
               loc="left", pad=4, fontsize=8.5)

# ============================================================================
# (c) Salinity vs Shannon — EQ XRF S + Atacama EC
# ============================================================================
ax_c = fig.add_subplot(gs[1, :])

# EQ — use XRF S as % dry mass
eq_sal = frame[(frame["S"].notna()) & (frame["S"] > 0)
               & (frame["shannon"].notna())].copy()
eq_z = (np.log10(eq_sal["S"]) - np.log10(eq_sal["S"]).mean()) / \
        np.log10(eq_sal["S"]).std()
ax_c.scatter(eq_z, eq_sal["shannon"], s=10, alpha=0.4,
             c=DESERT_COLORS["EmptyQuarter"],
             edgecolors="none", zorder=2,
             label=f"Empty Quarter — XRF S (n={len(eq_sal)})")

# Atacama — EC proxy
atac = xs[(xs["desert"] == "Atacama")
          & xs["electrical_conductivity"].notna()
          & (xs["electrical_conductivity"].astype(str).str.replace(".", "")
             .str.replace("-", "").str.isdigit())].copy()
atac["EC"] = pd.to_numeric(atac["electrical_conductivity"], errors="coerce")
atac = atac.dropna(subset=["EC"])
atac = atac[atac["EC"] > 0]
atac_z = (np.log10(atac["EC"]) - np.log10(atac["EC"]).mean()) / \
          np.log10(atac["EC"]).std()
ax_c.scatter(atac_z, atac["shannon_otu97"], s=18, alpha=0.9,
             c=DESERT_COLORS["Atacama"],
             edgecolors="black", linewidths=0.3, zorder=3,
             label=f"Atacama — EC (n={len(atac)})")

# Linear fits in z-space
def linfit(x, y):
    mask = np.isfinite(x) & np.isfinite(y)
    if mask.sum() < 5:
        return None
    p = np.polyfit(x[mask], y[mask], 1)
    xx = np.linspace(np.min(x[mask]), np.max(x[mask]), 50)
    return xx, np.polyval(p, xx), float(np.corrcoef(x[mask], y[mask])[0, 1])

for x_arr, y_arr, c, lbl, x_pos, y_pos, ha in [
    (eq_z.values, eq_sal["shannon"].values,
     DESERT_COLORS["EmptyQuarter"], "EQ", 0.04, 0.10, "left"),
    (atac_z.values, atac["shannon_otu97"].values,
     DESERT_COLORS["Atacama"], "Atacama", 0.96, 0.10, "right"),
]:
    res = linfit(x_arr, y_arr)
    if res:
        xx, yy, rho = res
        ax_c.plot(xx, yy, color=c, lw=1.4, ls="--", alpha=0.85, zorder=4)
        ax_c.text(x_pos, y_pos, f"{lbl} ρ = {rho:+.2f}",
                  transform=ax_c.transAxes, ha=ha, va="bottom",
                  fontsize=7, color=c, fontweight="bold")

ax_c.set_xlabel("Salinity proxy (per-desert log z-score)", fontsize=7.5)
ax_c.set_ylabel("Shannon diversity", fontsize=7.5)
ax_c.tick_params(direction="out", length=2.5)
for sp in ("top", "right"):
    ax_c.spines[sp].set_visible(False)
ax_c.legend(loc="upper right", frameon=False, fontsize=6.5,
            handletextpad=0.4, borderpad=0.3, scatterpoints=1)
ax_c.grid(True, alpha=0.2, lw=0.4)
ax_c.set_axisbelow(True)

ax_c.text(-0.13, 1.08, "c", transform=ax_c.transAxes,
          fontsize=12, fontweight="bold", va="top", ha="left")
ax_c.set_title("Salinity suppresses Shannon — two deserts",
               loc="left", pad=4, fontsize=8.5)

# ============================================================================
# (d) Ego networks centred on CSP1-2 per compartment
# ============================================================================
def load_net(comp):
    nodes = pd.read_csv(CACHE / f"network_nodes_{comp}.tsv", sep="\t")
    edges = pd.read_csv(CACHE / f"network_edges_{comp}.tsv", sep="\t")
    return nodes, edges

def draw_ego(ax, comp, label, top_n=12):
    nodes, edges = load_net(comp)
    # Filter to CSP1-2 + its top-N strongest neighbors
    csp_edges = edges[(edges["source"] == "CSP1-2") |
                      (edges["target"] == "CSP1-2")].copy()
    if csp_edges.empty:
        ax.text(0.5, 0.5, f"{label}\n(CSP1-2 not connected)",
                ha="center", va="center", transform=ax.transAxes,
                fontsize=7)
        ax.set_xticks([]); ax.set_yticks([])
        for sp in ("top", "right", "left", "bottom"):
            ax.spines[sp].set_visible(False)
        return
    csp_edges["nbr"] = np.where(csp_edges["source"] == "CSP1-2",
                                csp_edges["target"], csp_edges["source"])
    csp_edges["abs_w"] = csp_edges["weight"].abs()
    csp_edges = csp_edges.sort_values("abs_w", ascending=False).head(top_n)

    nbrs = list(csp_edges["nbr"])
    G = nx.Graph()
    G.add_node("CSP1-2")
    for _, r in csp_edges.iterrows():
        G.add_edge("CSP1-2", r["nbr"], weight=float(r["weight"]),
                   sign=int(r["sign"]))

    # Radial layout
    n = len(nbrs)
    pos = {"CSP1-2": (0.0, 0.0)}
    for i, nm in enumerate(nbrs):
        ang = 2 * np.pi * i / n - np.pi / 2
        pos[nm] = (np.cos(ang), np.sin(ang))

    # Edges
    for u, v, d in G.edges(data=True):
        x0, y0 = pos[u]; x1, y1 = pos[v]
        c = "#0077bb" if d["sign"] > 0 else "#cc3311"
        lw = 0.4 + 2.2 * abs(d["weight"])
        ax.plot([x0, x1], [y0, y1], color=c, lw=lw, alpha=0.75,
                solid_capstyle="round", zorder=2)

    # Hub node — moderate dot, label sits just to the upper-left
    ax.scatter(0, 0, s=200, c=C_KEYSTONE, edgecolors="black",
               linewidths=0.8, zorder=4)
    ax.text(0.0, -0.22, "CSP1-2", ha="center", va="top",
            fontsize=7.5, fontweight="bold", color=C_KEYSTONE, zorder=5,
            bbox=dict(boxstyle="round,pad=0.18", fc="white",
                      ec=C_KEYSTONE, lw=0.5))

    # Neighbor nodes
    for nm in nbrs:
        x, y_ = pos[nm]
        ax.scatter(x, y_, s=42, c="#888", edgecolors="black",
                   linewidths=0.4, zorder=4)
        # Label outside the dot
        r = 1.18
        lx, ly = r * x, r * y_
        ha = "left" if x > 0.05 else ("right" if x < -0.05 else "center")
        va = "bottom" if y_ > 0.05 else ("top" if y_ < -0.05 else "center")
        ax.text(lx, ly, nm.replace("_", " "), ha=ha, va=va,
                fontsize=5.6, color="#222")

    ax.set_xlim(-1.7, 1.7); ax.set_ylim(-1.55, 1.55)
    ax.set_xticks([]); ax.set_yticks([])
    for sp in ("top", "right", "left", "bottom"):
        ax.spines[sp].set_visible(False)
    ax.set_title(label, loc="center", fontsize=8, pad=2, fontweight="bold")
    ax.set_aspect("equal")

ax_d1 = fig.add_subplot(gs[2, 0:2])
ax_d2 = fig.add_subplot(gs[2, 2:4])
ax_d3 = fig.add_subplot(gs[2, 4:6])

draw_ego(ax_d1, "surface", "Surface")
draw_ego(ax_d2, "deep", "Deep")
draw_ego(ax_d3, "rhizosphere", "Rhizosphere")

# Row label for panel d
fig.text(0.05, 0.30, "d",
         fontsize=12, fontweight="bold", va="bottom", ha="left")
fig.text(0.10, 0.302,
         "CSP1-2 ego networks per compartment (top-12 |ρ| neighbours; "
         "blue = positive, red = negative co-occurrence)",
         fontsize=8.5, va="bottom", ha="left")

# ============================================================================
out_pdf = FIG / "fig3_csp12.pdf"
out_png = FIG / "fig3_csp12.png"
fig.savefig(out_pdf, bbox_inches="tight")
fig.savefig(out_png, bbox_inches="tight", dpi=300)
print(f"wrote {out_pdf} ({out_pdf.stat().st_size} bytes)")
print(f"wrote {out_png} ({out_png.stat().st_size} bytes)")
