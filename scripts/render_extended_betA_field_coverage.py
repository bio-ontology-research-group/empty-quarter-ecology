#!/usr/bin/env python3
"""Extended Data figure: betA gene coverage in metagenomes vs salinity.

Reads:
  cache/betA_field_coverage.tsv  (produced via tblastn on unimatrix01)

Panels:
  (a) betA hits (e<1e-30) vs salinity (XRF S, log scale).
  (b) betA hits vs CSP1-2 16S relative abundance.
  (c) betA hits vs Shannon.
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt

mpl.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
    "font.size": 8,
    "axes.titlesize": 9, "axes.labelsize": 8,
    "xtick.labelsize": 7, "ytick.labelsize": 7,
    "legend.fontsize": 7, "axes.linewidth": 0.7,
    "pdf.fonttype": 42, "ps.fonttype": 42,
})

REPO = Path(__file__).resolve().parents[1]
CACHE = REPO / "cache"
FIG = REPO.parent / "figures"

cov = pd.read_csv(CACHE / "betA_field_coverage.tsv", sep="\t")
print(f"Samples with betA coverage: {len(cov)}")

cov["S_pos"] = cov["S"].clip(lower=1e-3)

fig, axes = plt.subplots(1, 3, figsize=(7.2, 2.4))

# Panel (a) S vs hits
ax = axes[0]
ax.scatter(cov["S_pos"], cov["n_hits_e30"], s=42, c="#cc3311",
           edgecolor="black", linewidth=0.5)
ax.set_xscale("log")
ax.set_xlabel("Salinity (XRF S, % dry mass)")
ax.set_ylabel("betA tBLASTn hits (e<1e-30)")
for sp in ("top", "right"):
    ax.spines[sp].set_visible(False)
ax.tick_params(direction="out", length=2.5)
ax.grid(True, alpha=0.25, lw=0.4)
ax.set_axisbelow(True)
ax.text(-0.20, 1.10, "a", transform=ax.transAxes,
        fontsize=12, fontweight="bold", va="top", ha="left")
ax.set_title("Salinity vs betA gene hits", loc="left", pad=4, fontsize=8.5)
# spearman
from scipy.stats import spearmanr
rho_a, p_a = spearmanr(cov["S_pos"], cov["n_hits_e30"])
ax.text(0.05, 0.95, f"ρ = {rho_a:+.2f}\np = {p_a:.2g}",
        transform=ax.transAxes, ha="left", va="top",
        fontsize=7.5, fontweight="bold")

# Panel (b) CSP1-2 abundance vs hits
ax = axes[1]
ax.scatter(cov["csp_relab"].clip(lower=1e-6) * 100, cov["n_hits_e30"],
           s=42, c="#0077bb", edgecolor="black", linewidth=0.5)
ax.set_xscale("log")
ax.set_xlabel("CSP1-2 16S relative abundance (%)")
ax.set_ylabel("betA tBLASTn hits (e<1e-30)")
for sp in ("top", "right"):
    ax.spines[sp].set_visible(False)
ax.tick_params(direction="out", length=2.5)
ax.grid(True, alpha=0.25, lw=0.4)
ax.set_axisbelow(True)
ax.text(-0.20, 1.10, "b", transform=ax.transAxes,
        fontsize=12, fontweight="bold", va="top", ha="left")
ax.set_title("CSP1-2 16S vs betA hits", loc="left", pad=4, fontsize=8.5)
rho_b, p_b = spearmanr(cov["csp_relab"].clip(lower=1e-6),
                        cov["n_hits_e30"])
ax.text(0.05, 0.95, f"ρ = {rho_b:+.2f}\np = {p_b:.2g}",
        transform=ax.transAxes, ha="left", va="top",
        fontsize=7.5, fontweight="bold")

# Panel (c) Shannon vs hits
ax = axes[2]
ax.scatter(cov["shannon"], cov["n_hits_e30"], s=42, c="#117733",
           edgecolor="black", linewidth=0.5)
ax.set_xlabel("Community Shannon")
ax.set_ylabel("betA tBLASTn hits (e<1e-30)")
for sp in ("top", "right"):
    ax.spines[sp].set_visible(False)
ax.tick_params(direction="out", length=2.5)
ax.grid(True, alpha=0.25, lw=0.4)
ax.set_axisbelow(True)
ax.text(-0.20, 1.10, "c", transform=ax.transAxes,
        fontsize=12, fontweight="bold", va="top", ha="left")
ax.set_title("Community Shannon vs betA hits", loc="left", pad=4, fontsize=8.5)
rho_c, p_c = spearmanr(cov["shannon"], cov["n_hits_e30"])
ax.text(0.05, 0.95, f"ρ = {rho_c:+.2f}\np = {p_c:.2g}",
        transform=ax.transAxes, ha="left", va="top",
        fontsize=7.5, fontweight="bold")

plt.tight_layout()
out_pdf = FIG / "extended_fig_betA_field_coverage.pdf"
out_png = FIG / "extended_fig_betA_field_coverage.png"
fig.savefig(out_pdf, bbox_inches="tight")
fig.savefig(out_png, bbox_inches="tight", dpi=300)
print(f"wrote {out_pdf}")
print(f"wrote {out_png}")
