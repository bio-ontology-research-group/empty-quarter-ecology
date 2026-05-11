#!/usr/bin/env python3
"""Hill cross-cohort prediction figure.

Shows: pre-registered EQ Hill curve (solid) vs free Atacama Hill (dashed)
vs Atacama linear (dotted) vs Atacama observations.
"""
from __future__ import annotations
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt

mpl.rcParams.update({
    "font.family":"sans-serif",
    "font.sans-serif":["Arial","Helvetica","DejaVu Sans"],
    "font.size":8, "axes.titlesize":9, "axes.labelsize":8,
    "xtick.labelsize":7, "ytick.labelsize":7, "legend.fontsize":7,
    "axes.linewidth":0.7, "pdf.fonttype":42, "ps.fonttype":42,
})

REPO = Path(__file__).resolve().parents[1]
CACHE = REPO / "cache"
FIG = REPO.parent / "figures"

# Atacama data
df = pd.read_csv(CACHE/"crossdesert"/"per_sample.tsv", sep="\t")
ata = df[df["desert"]=="Atacama"].copy()
ata["EC"] = pd.to_numeric(ata["electrical_conductivity"], errors="coerce")
m = ata.dropna(subset=["EC","shannon_otu97"])

# Fits from cache (re-run quickly)
fits = pd.read_csv(CACHE/"hill_cross_cohort_atacama.tsv", sep="\t")

EQ_K, EQ_n = 0.02, 9.4
EQ_K_pool, EQ_n_pool = 0.09, 4.6
median_EC = m["EC"].clip(lower=1e-3).median()
EQ_med_S = 0.10  # rough EQ within-site median

# Plot
fig, axes = plt.subplots(1, 2, figsize=(7.0, 3.0))

# Panel a: EQ pre-registered curve in z(salinity) units, vs Atacama curve
ax = axes[0]
S_grid = np.logspace(-3, 1, 200)

# EQ within-site Hill curve: Vmax * K^n / (K^n + S^n)
def hill(S, Vmax, K, n):
    return Vmax * (K**n) / (K**n + S**n)

eq_curve = hill(S_grid, 6, EQ_K, EQ_n)
ax.plot(S_grid, eq_curve, color="#0077bb", lw=1.6,
         label=f"EQ within-site Hill\nK={EQ_K}%, n={EQ_n}")
eq_curve_pool = hill(S_grid, 6, EQ_K_pool, EQ_n_pool)
ax.plot(S_grid, eq_curve_pool, color="#33aabb", lw=1.2, ls="--",
         label=f"EQ pooled Hill\nK={EQ_K_pool}%, n={EQ_n_pool}")
ax.set_xscale("log")
ax.set_xlabel("Salinity (% S, EQ XRF — pre-registered scale)")
ax.set_ylabel("Predicted Shannon")
ax.legend(loc="lower left", frameon=False, fontsize=6.5)
for sp in ("top","right"):
    ax.spines[sp].set_visible(False)
ax.tick_params(direction="out", length=2.5)
ax.grid(True, alpha=0.25, lw=0.4)
ax.set_axisbelow(True)
ax.text(-0.15, 1.08, "a", transform=ax.transAxes, fontweight="bold", fontsize=12)
ax.set_title("Pre-registered EQ Hill curve", loc="left", pad=4)

# Panel b: Atacama free Hill vs linear vs data
ax = axes[1]
EC_grid = np.logspace(-3, 1, 200)
ax.scatter(m["EC"].clip(lower=1e-3), m["shannon_otu97"], s=42,
            c="#bb5566", edgecolor="black", linewidth=0.4, alpha=0.85,
            label=f"Atacama (n={len(m)})", zorder=3)

# Free Hill from fits table
h1 = fits[fits["hypothesis"]=="H1_free_Hill"].iloc[0]
free_curve = hill(EC_grid, h1["Vmax"], h1["K"], h1["n_Hill"])
ax.plot(EC_grid, free_curve, color="black", lw=1.4,
         label=f"Atacama free Hill\nK={h1['K']:.2f} dS/m, n={h1['n_Hill']:.2f}")

# Linear fit
slope, inter = np.polyfit(m["EC"].clip(lower=1e-3), m["shannon_otu97"], 1)
lin_curve = inter + slope * EC_grid
ax.plot(EC_grid, lin_curve, color="gray", lw=1.0, ls=":",
         label=f"Atacama linear\nb={slope:+.2f}")

ax.set_xscale("log")
ax.set_xlim(1e-3, 5)
ax.set_ylim(0, 8.5)
ax.set_xlabel("Electrical conductivity (dS/m)")
ax.set_ylabel("Shannon (OTU97)")
ax.legend(loc="lower left", frameon=False, fontsize=6.5)
for sp in ("top","right"):
    ax.spines[sp].set_visible(False)
ax.tick_params(direction="out", length=2.5)
ax.grid(True, alpha=0.25, lw=0.4)
ax.set_axisbelow(True)
ax.text(-0.15, 1.08, "b", transform=ax.transAxes, fontweight="bold", fontsize=12)
ax.set_title("Atacama: pre-registered shape does NOT transfer", loc="left", pad=4)

# Annotate failure
h0_aic = fits[fits["hypothesis"]=="H0_linear"].iloc[0]["AIC"]
h1_aic = fits[fits["hypothesis"]=="H1_free_Hill"].iloc[0]["AIC"]
delta = h0_aic - h1_aic
ax.text(0.05, 0.20, f"ΔAIC(linear-free Hill) = {delta:+.1f}\nfree n_Atacama = {h1['n_Hill']:.2f} (vs EQ n=9.4)\nLinear is preferred",
         transform=ax.transAxes, ha="left", va="top", fontsize=6.5,
         bbox=dict(boxstyle="round,pad=0.3", facecolor="#fff3e6",
                    edgecolor="#cc6611", lw=0.6))

plt.tight_layout()
out_pdf = FIG / "extended_fig_hill_cross_cohort.pdf"
out_png = FIG / "extended_fig_hill_cross_cohort.png"
fig.savefig(out_pdf, bbox_inches="tight")
fig.savefig(out_png, bbox_inches="tight", dpi=300)
print(f"wrote {out_pdf}")
print(f"wrote {out_png}")
