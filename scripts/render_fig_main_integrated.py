#!/usr/bin/env python3
"""Render integrated forecast figure — combines Fig 4 (restoration)
and Fig 5 (climate damage) into one composite.

Panels:
  (a) Salinity → Shannon Hill curve (mechanism, top-left).
  (b) Per-compartment intervention posteriors: do(S −1 SD), do(P +1 SD),
      do(reclam) overlaid with CMIP6 climate-amplified salinity damage.
  (c) Time horizon comparison — baseline → 2050 climate damage → 2100
      damage → reclamation gain — per compartment.

Reads:
  cache/causal_tier3_interventions_expanded.tsv
  cache/cmip6_interventions.tsv
  cache/causal_frame_tier1.parquet
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
    "xtick.major.width": 0.7,
    "ytick.major.width": 0.7,
    "lines.linewidth": 1.0,
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
})

REPO = Path(__file__).resolve().parents[1]
CACHE = REPO / "cache"
FIG = REPO.parent / "figures"

COMP_COLOR = {"surface": "#ee7733", "deep": "#0077bb",
              "rhizosphere": "#117733"}

frame = pd.read_parquet(CACHE / "causal_frame_tier1.parquet")
inter = pd.read_csv(CACHE / "causal_tier3_interventions_expanded.tsv", sep="\t")
cmip6 = pd.read_csv(CACHE / "cmip6_interventions.tsv", sep="\t")

# ============================================================================
fig = plt.figure(figsize=(7.2, 7.6))
gs = GridSpec(2, 2, figure=fig,
              height_ratios=[1.0, 1.05],
              hspace=0.50, wspace=0.40,
              left=0.10, right=0.97, top=0.93, bottom=0.07)

# ============================================================================
# (a) Hill curve recap
# ============================================================================
ax_a = fig.add_subplot(gs[0, 0])
cells = (frame.groupby(["trip", "site", "compartment"])
         .agg(S=("S", "mean"), csp=("csp_relab", "mean"))
         .reset_index())
cells = cells.dropna(subset=["S", "csp"]).query("S > 0").copy()

for comp, c in COMP_COLOR.items():
    sub = cells[cells["compartment"] == comp]
    ax_a.scatter(sub["S"], sub["csp"] * 100, s=12, alpha=0.55,
                 c=c, edgecolors="none", label=comp.title(), zorder=2)

M0, V, K, n_hill = 0.0046, 0.0035, 0.093, 4.6
xx = np.logspace(-3, 1.3, 400)
yy = (M0 - V * xx ** n_hill / (K ** n_hill + xx ** n_hill)) * 100
ax_a.plot(xx, yy, color="black", lw=1.6, zorder=4,
          label=f"Hill (n={n_hill:.1f}, K={K:.2f}%)")
ax_a.axvline(K, color="#666", lw=0.6, ls=(0, (2, 1.5)), alpha=0.7)
ax_a.text(K * 1.15, 0.05,
          f"threshold S ≈ {K:.2f}% dry mass",
          fontsize=6.8, color="#444", style="italic")

ax_a.set_xscale("log")
ax_a.set_xlabel("Salinity (XRF S, % dry mass)", fontsize=7.5)
ax_a.set_ylabel("CSP1-2 relative abundance (%)", fontsize=7.5)
ax_a.set_xlim(0.001, 20)
ax_a.set_ylim(0, 1.0)
ax_a.tick_params(direction="out", length=2.5)
for sp in ("top", "right"):
    ax_a.spines[sp].set_visible(False)
ax_a.legend(loc="upper right", frameon=False, fontsize=6.3,
            handletextpad=0.3, borderpad=0.3, scatterpoints=1,
            labelspacing=0.3)

ax_a.text(-0.16, 1.07, "a", transform=ax_a.transAxes,
          fontsize=12, fontweight="bold", va="top", ha="left")
ax_a.set_title("Mechanism: salinity threshold collapses CSP1-2",
               loc="left", pad=4, fontsize=8.5)

# ============================================================================
# (b) Intervention forest plot — per compartment, with climate row
# ============================================================================
ax_b = fig.add_subplot(gs[0, 1])

interventions = [
    ("Climate (SSP3-7.0, 2100)\n+ salinity amplification",
     None, None, "climate", "#888"),
    ("Reclamation\ndo(S −2 SD, P +1 SD)",
     "do(S -2 SD, P +1 SD) [reclam]", None, "intervention", "#117733"),
    ("Desalination only\ndo(S −1 SD)",
     "do(S -1 SD)", None, "intervention", "#0077bb"),
    ("Phosphorus only\ndo(P +1 SD)",
     "do(P +1 SD)", None, "intervention", "#cc8855"),
]
CLIMATE_SCEN = "SSP3-7.0"  # high-end available in cache

comp_order = ["surface", "deep", "rhizosphere"]
y_offsets = {"surface": -0.22, "deep": 0.0, "rhizosphere": 0.22}

for i, (lbl, intkey, _, kind, c) in enumerate(interventions):
    if kind == "intervention":
        for comp in comp_order:
            row = inter[(inter["intervention"] == intkey)
                        & (inter["compartment"] == comp)]
            if row.empty:
                continue
            m = float(row.iloc[0]["median_z"])
            lo = float(row.iloc[0]["ci_lo"])
            hi = float(row.iloc[0]["ci_hi"])
            cred = bool(row.iloc[0]["credible_95"])
            y = i + y_offsets[comp]
            cc = COMP_COLOR[comp]
            ax_b.errorbar(m, y, xerr=[[m - lo], [hi - m]],
                          fmt="o", color=cc, ecolor=cc, ms=6,
                          mec="black", mew=0.4, lw=1.0, capsize=2,
                          mfc=cc if cred else "white")
    elif kind == "climate":
        # Use highest-available SSP at 2100 with_S_amp pathway per compartment
        ssp = cmip6[(cmip6["scenario"] == CLIMATE_SCEN)
                    & (cmip6["horizon"].astype(str) == "2100")
                    & (cmip6["pathway"] == "with_S_amp")]
        for comp in comp_order:
            row = ssp[ssp["compartment"] == comp]
            if row.empty:
                continue
            m = float(row.iloc[0]["delta_shannon_z_median"])
            lo = float(row.iloc[0]["ci_lo"])
            hi = float(row.iloc[0]["ci_hi"])
            cred = bool(row.iloc[0]["credible_95"])
            y = i + y_offsets[comp]
            cc = COMP_COLOR[comp]
            ax_b.errorbar(m, y, xerr=[[m - lo], [hi - m]],
                          fmt="s", color=cc, ecolor=cc, ms=5,
                          mec="black", mew=0.4, lw=1.0, capsize=2,
                          mfc=cc if cred else "white")

ax_b.axvline(0, color="black", lw=0.6, alpha=0.5, ls=(0, (1, 1.2)))
ax_b.set_yticks(range(len(interventions)))
ax_b.set_yticklabels([t[0] for t in interventions], fontsize=7)
ax_b.invert_yaxis()
ax_b.set_xlabel("Posterior median ΔShannon (z, 95% CI)", fontsize=7.5)
ax_b.set_xlim(-0.7, 1.0)
ax_b.tick_params(direction="out", length=2.5)
for sp in ("top", "right"):
    ax_b.spines[sp].set_visible(False)
ax_b.grid(True, axis="x", alpha=0.25, lw=0.4)
ax_b.set_axisbelow(True)

# Mini compartment legend
for comp, c in COMP_COLOR.items():
    ax_b.scatter([], [], c=c, s=36, edgecolors="black", linewidths=0.4,
                 label=comp.title())
ax_b.legend(loc="upper right", frameon=False, fontsize=6.3,
            handletextpad=0.3, borderpad=0.3, scatterpoints=1,
            labelspacing=0.3)

ax_b.text(-0.32, 1.07, "b", transform=ax_b.transAxes,
          fontsize=12, fontweight="bold", va="top", ha="left")
ax_b.set_title("Climate damage vs reclamation intervention",
               loc="left", pad=4, fontsize=8.5)

# ============================================================================
# (c) Time-horizon trajectory per compartment
# ============================================================================
ax_c = fig.add_subplot(gs[1, :])

# Compute trajectories per compartment per scenario
horizons = [2025, 2050, 2100]
scenarios = ["SSP1-2.6", "SSP2-4.5", "SSP3-7.0"]
sc_color = {"SSP1-2.6": "#0077bb", "SSP2-4.5": "#cc8855",
            "SSP3-7.0": "#cc3311"}

# Get reclamation effect: stays constant (intervention applied at t=2025)
reclam_per_comp = {}
for comp in comp_order:
    row = inter[(inter["intervention"] == "do(S -2 SD, P +1 SD) [reclam]")
                & (inter["compartment"] == comp)]
    reclam_per_comp[comp] = (float(row.iloc[0]["median_z"]),
                             float(row.iloc[0]["ci_lo"]),
                             float(row.iloc[0]["ci_hi"]))

# Plot one panel per compartment side-by-side
n_comp = len(comp_order)
inner = gs[1, :].subgridspec(1, n_comp, wspace=0.35)
for ci, comp in enumerate(comp_order):
    ax_inner = fig.add_subplot(inner[0, ci])

    # Baseline at 0 (z-anchored)
    ax_inner.scatter(2025, 0, s=80, color="black", marker="o", zorder=5,
                     label="Baseline 2025")
    ax_inner.errorbar(2025, 0, yerr=[[0], [0]], fmt="none", color="black")

    # CMIP6 trajectories with_S_amp
    for sc in scenarios:
        ys, ylos, yhis = [0], [0], [0]
        for h in (2050, 2100):
            row = cmip6[(cmip6["scenario"] == sc)
                        & (cmip6["horizon"].astype(str) == str(h))
                        & (cmip6["compartment"] == comp)
                        & (cmip6["pathway"] == "with_S_amp")]
            if not row.empty:
                ys.append(float(row.iloc[0]["delta_shannon_z_median"]))
                ylos.append(float(row.iloc[0]["ci_lo"]))
                yhis.append(float(row.iloc[0]["ci_hi"]))
            else:
                ys.append(np.nan); ylos.append(np.nan); yhis.append(np.nan)
        ax_inner.plot([2025, 2050, 2100], ys, "-o", lw=1.3,
                      ms=4, color=sc_color[sc],
                      mec="black", mew=0.4,
                      label=sc if ci == 0 else None)
        ax_inner.fill_between([2025, 2050, 2100], ylos, yhis,
                              alpha=0.15, color=sc_color[sc])

    # Reclamation reference line (applied as offset to baseline)
    rm, rlo, rhi = reclam_per_comp[comp]
    ax_inner.axhline(rm, color=COMP_COLOR[comp], lw=1.3, ls=(0, (4, 2)),
                     alpha=0.85, label="Reclamation gain" if ci == 0 else None)
    ax_inner.fill_between([2020, 2105], rlo, rhi,
                          alpha=0.10, color=COMP_COLOR[comp])
    ax_inner.text(2103, rm, f"{rm:+.2f}σ",
                  ha="right", va="bottom", fontsize=7,
                  color=COMP_COLOR[comp], fontweight="bold")

    ax_inner.axhline(0, color="black", lw=0.5, alpha=0.5)
    ax_inner.set_xlim(2020, 2105)
    ax_inner.set_ylim(-0.85, 0.85)
    ax_inner.set_xticks([2025, 2050, 2100])
    ax_inner.set_xlabel("Year", fontsize=7.5)
    if ci == 0:
        ax_inner.set_ylabel("ΔShannon (z, vs 2025 baseline)", fontsize=7.5)
    ax_inner.set_title(comp.title(), loc="left", pad=4, fontsize=8.5,
                       fontweight="bold")
    ax_inner.tick_params(direction="out", length=2.5)
    for sp in ("top", "right"):
        ax_inner.spines[sp].set_visible(False)
    if ci == 0:
        ax_inner.legend(loc="lower left", frameon=False, fontsize=6.0,
                        handletextpad=0.5, borderpad=0.3, labelspacing=0.3)

# Panel-c label
fig.text(0.05, 0.46, "c",
         fontsize=12, fontweight="bold", va="top", ha="left")
fig.text(0.10, 0.46,
         "Climate damage vs reclamation gain — per compartment, century horizon",
         fontsize=8.5, va="top", ha="left")

# Save
out_pdf = FIG / "fig5_integrated_forecast.pdf"
out_png = FIG / "fig5_integrated_forecast.png"
fig.savefig(out_pdf, bbox_inches="tight")
fig.savefig(out_png, bbox_inches="tight", dpi=300)
print(f"wrote {out_pdf} ({out_pdf.stat().st_size} bytes)")
print(f"wrote {out_png} ({out_png.stat().st_size} bytes)")
