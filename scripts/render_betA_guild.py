#!/usr/bin/env python3
"""Figure: EQ betA-encoder guild map across 312 metagenomes.

Panels:
  (a) Per-compartment producer-bin abundance: bins encoding K00108
       per sample (boxplot + jitter, by compartment).
  (b) Per-sample K00108 density (per million CDS) by compartment.
  (c) Spatial: lat/lon scatter coloured by guild density (if metadata).
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

samp = pd.read_csv(CACHE / "betA_per_sample_summary.tsv", sep="\t")

# Compartment from sample id
import re
def comp_of(sid):
    s = sid.split("_")[-1] if "_" in sid else sid
    m = re.match(r"^[0-9]+([A-Z]+)r[0-9]+$", s)
    if not m: return "?"
    return {"PR":"rhizosphere","S":"surface","D":"deep"}.get(m.group(1), m.group(1))
samp["compartment"] = samp["sample_id"].apply(comp_of)

# COLORS
CMAP = {"surface":"#cc6677","rhizosphere":"#117733","deep":"#332288"}
ORDER = ["surface","rhizosphere","deep"]

fig, axes = plt.subplots(1, 3, figsize=(8.6, 2.8))

# Panel (a): bins-with-K00108 per sample by compartment
ax = axes[0]
data_a = [samp[samp["compartment"]==c]["n_bins_with_K00108"].values for c in ORDER]
bp = ax.boxplot(data_a, positions=range(len(ORDER)), widths=0.45,
                 patch_artist=True, showfliers=False)
for i, box in enumerate(bp["boxes"]):
    box.set(facecolor=CMAP[ORDER[i]], alpha=0.45, edgecolor="black", linewidth=0.7)
for med in bp["medians"]:
    med.set(color="white", linewidth=1.6)
rng = np.random.default_rng(7)
for i, arr in enumerate(data_a):
    x = i + rng.uniform(-0.12, 0.12, len(arr))
    ax.scatter(x, arr, s=10, color=CMAP[ORDER[i]], edgecolor="black",
                linewidth=0.3, alpha=0.7, zorder=3)
ax.set_xticks(range(len(ORDER)))
ax.set_xticklabels(ORDER, rotation=15)
ax.set_ylabel("# bins encoding K00108 per sample")
for sp in ("top","right"):
    ax.spines[sp].set_visible(False)
ax.tick_params(direction="out", length=2.5)
ax.grid(True, alpha=0.25, lw=0.4, axis="y")
ax.set_axisbelow(True)
ax.text(-0.15, 1.10, "a", transform=ax.transAxes, fontweight="bold", fontsize=12)
ax.set_title("Producer-bin count", loc="left", pad=4)

# Panel (b): K00108 density per million CDS
ax = axes[1]
data_b = [samp[samp["compartment"]==c]["K00108_density_per_M_CDS"].values for c in ORDER]
bp = ax.boxplot(data_b, positions=range(len(ORDER)), widths=0.45,
                 patch_artist=True, showfliers=False)
for i, box in enumerate(bp["boxes"]):
    box.set(facecolor=CMAP[ORDER[i]], alpha=0.45, edgecolor="black", linewidth=0.7)
for med in bp["medians"]:
    med.set(color="white", linewidth=1.6)
for i, arr in enumerate(data_b):
    x = i + rng.uniform(-0.12, 0.12, len(arr))
    ax.scatter(x, arr, s=10, color=CMAP[ORDER[i]], edgecolor="black",
                linewidth=0.3, alpha=0.7, zorder=3)
ax.set_xticks(range(len(ORDER)))
ax.set_xticklabels(ORDER, rotation=15)
ax.set_ylabel("K00108 hits per million CDS")
for sp in ("top","right"):
    ax.spines[sp].set_visible(False)
ax.tick_params(direction="out", length=2.5)
ax.grid(True, alpha=0.25, lw=0.4, axis="y")
ax.set_axisbelow(True)
ax.text(-0.15, 1.10, "b", transform=ax.transAxes, fontweight="bold", fontsize=12)
ax.set_title("Functional density", loc="left", pad=4)

# Panel (c): producer prevalence per compartment (% samples with at least 1)
ax = axes[2]
prev_data = []
for c in ORDER:
    sub = samp[samp["compartment"]==c]
    prev = (sub["n_K00108_strict"]>0).mean() * 100
    n = len(sub)
    prev_data.append({"compartment":c, "prevalence_%":prev, "n":n})
df_prev = pd.DataFrame(prev_data)
ax.bar(range(len(ORDER)), df_prev["prevalence_%"],
        color=[CMAP[c] for c in ORDER], edgecolor="black", linewidth=0.6)
ax.set_xticks(range(len(ORDER)))
ax.set_xticklabels(ORDER, rotation=15)
ax.set_ylabel("Samples with $\\geq$1 K00108 producer (%)")
for i, row in df_prev.iterrows():
    ax.text(i, row["prevalence_%"]+1, f"{row['prevalence_%']:.0f}%\n(n={row['n']})",
            ha="center", va="bottom", fontsize=7)
ax.set_ylim(0, 105)
for sp in ("top","right"):
    ax.spines[sp].set_visible(False)
ax.tick_params(direction="out", length=2.5)
ax.grid(True, alpha=0.25, lw=0.4, axis="y")
ax.set_axisbelow(True)
ax.text(-0.15, 1.10, "c", transform=ax.transAxes, fontweight="bold", fontsize=12)
ax.set_title("Producer prevalence", loc="left", pad=4)

plt.tight_layout()
out_pdf = FIG / "extended_fig_betA_guild.pdf"
out_png = FIG / "extended_fig_betA_guild.png"
fig.savefig(out_pdf, bbox_inches="tight")
fig.savefig(out_png, bbox_inches="tight", dpi=300)
print(f"wrote {out_pdf}\nwrote {out_png}")
