#!/usr/bin/env python3
"""Atacama within-desert figure: positive variables explaining CSP1-2 presence.

Panels:
  (a) Box-violin: elevation in CSP+ vs CSP-
  (b) Box-violin: soil RH in CSP+ vs CSP-
  (c) Box-violin: EC in CSP+ vs CSP-
  (d) Map / scatter: elevation vs RH, coloured by CSP+/-
"""
from __future__ import annotations

from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt

mpl.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial","Helvetica","DejaVu Sans"],
    "font.size": 8, "axes.titlesize":9, "axes.labelsize":8,
    "xtick.labelsize":7, "ytick.labelsize":7, "legend.fontsize":7,
    "axes.linewidth":0.7, "pdf.fonttype":42, "ps.fonttype":42,
})

REPO = Path(__file__).resolve().parents[1]
CACHE = REPO / "cache"
FIG = REPO.parent / "figures"

df = pd.read_csv(CACHE / "crossdesert" / "per_sample.tsv", sep="\t")
ata = df[df["desert"]=="Atacama"].copy()
for c in ["elevation_m","avg_soil_rh","avg_soil_temp",
          "electrical_conductivity","ph","soil_organic_carbon"]:
    ata[c] = pd.to_numeric(ata[c], errors="coerce")
ata["has_csp"] = (ata["csp_rel_85"]>0).astype(int)

CPOS = "#117733"
CNEG = "#bb5566"

def violin_pair(ax, ata, col, ylab, ylog=False):
    pos = ata[ata["has_csp"]==1][col].dropna().values
    neg = ata[ata["has_csp"]==0][col].dropna().values
    parts = ax.violinplot([neg, pos], positions=[0,1], widths=0.7,
                           showextrema=False, showmedians=False)
    for i,pc in enumerate(parts["bodies"]):
        pc.set_facecolor([CNEG, CPOS][i])
        pc.set_alpha(0.4)
        pc.set_edgecolor("black")
        pc.set_linewidth(0.7)
    # boxplot overlay
    bp = ax.boxplot([neg,pos], positions=[0,1], widths=0.18,
                     patch_artist=True, showfliers=False, zorder=3)
    for i,box in enumerate(bp["boxes"]):
        box.set(facecolor=[CNEG,CPOS][i], edgecolor="black", linewidth=0.7)
    for med in bp["medians"]:
        med.set(color="white", linewidth=1.6)
    # jitter points
    rng = np.random.default_rng(7)
    for i,(arr,col2) in enumerate(zip([neg,pos],[CNEG,CPOS])):
        x = i + rng.uniform(-0.05,0.05,len(arr))
        ax.scatter(x, arr, s=10, color=col2, edgecolor="black",
                    linewidth=0.3, alpha=0.7, zorder=4)
    ax.set_xticks([0,1])
    ax.set_xticklabels([f"CSP-\n(n={len(neg)})", f"CSP+\n(n={len(pos)})"])
    ax.set_ylabel(ylab)
    if ylog: ax.set_yscale("log")
    for sp in ("top","right"):
        ax.spines[sp].set_visible(False)
    ax.tick_params(direction="out", length=2.5)
    ax.grid(True, alpha=0.25, lw=0.4, axis="y")
    ax.set_axisbelow(True)
    # MW p
    from scipy.stats import mannwhitneyu
    if len(pos)>=3 and len(neg)>=3:
        u,p = mannwhitneyu(pos, neg, alternative="two-sided")
        ax.text(0.5, 0.97, f"p = {p:.2g}", transform=ax.transAxes,
                ha="center", va="top", fontsize=7, fontweight="bold")

fig, axes = plt.subplots(1, 4, figsize=(8.5, 2.6))

violin_pair(axes[0], ata, "elevation_m", "Elevation (m)")
axes[0].axhline(3000, color="gray", lw=0.5, ls=":")
axes[0].text(-0.15, 1.10, "a", transform=axes[0].transAxes,
              fontsize=12, fontweight="bold")
axes[0].set_title("Elevation", loc="left", pad=4)

violin_pair(axes[1], ata, "avg_soil_rh", "Soil RH (%)")
axes[1].text(-0.15, 1.10, "b", transform=axes[1].transAxes,
              fontsize=12, fontweight="bold")
axes[1].set_title("Soil RH", loc="left", pad=4)

violin_pair(axes[2], ata, "electrical_conductivity",
             "Electrical conductivity (dS/m)", ylog=True)
axes[2].text(-0.15, 1.10, "c", transform=axes[2].transAxes,
              fontsize=12, fontweight="bold")
axes[2].set_title("EC", loc="left", pad=4)

# Panel (d): elevation vs RH scatter coloured by has_csp
ax = axes[3]
for label, sub, col in [("CSP-", ata[ata["has_csp"]==0], CNEG),
                          ("CSP+", ata[ata["has_csp"]==1], CPOS)]:
    ax.scatter(sub["elevation_m"], sub["avg_soil_rh"],
                c=col, edgecolor="black", linewidth=0.4,
                s=42, alpha=0.85, label=label)
ax.set_xlabel("Elevation (m)")
ax.set_ylabel("Soil RH (%)")
ax.legend(loc="lower right", frameon=False)
for sp in ("top","right"):
    ax.spines[sp].set_visible(False)
ax.tick_params(direction="out", length=2.5)
ax.grid(True, alpha=0.25, lw=0.4)
ax.set_axisbelow(True)
ax.text(-0.15, 1.10, "d", transform=ax.transAxes,
         fontsize=12, fontweight="bold")
ax.set_title("Niche space", loc="left", pad=4)
# annotate the Altiplano cluster
ax.annotate("Altiplano\n(CSP+)", xy=(3700, 100), xytext=(2800, 105),
             fontsize=7, ha="center",
             arrowprops=dict(arrowstyle="->", lw=0.6, color="black"))
ax.annotate("Hyperarid core\n(CSP-)", xy=(1200, 17), xytext=(1800, 35),
             fontsize=7, ha="center",
             arrowprops=dict(arrowstyle="->", lw=0.6, color="black"))

plt.tight_layout()
out_pdf = FIG / "extended_fig_atacama_within_desert.pdf"
out_png = FIG / "extended_fig_atacama_within_desert.png"
fig.savefig(out_pdf, bbox_inches="tight")
fig.savefig(out_png, bbox_inches="tight", dpi=300)
print(f"wrote {out_pdf}")
print(f"wrote {out_png}")
