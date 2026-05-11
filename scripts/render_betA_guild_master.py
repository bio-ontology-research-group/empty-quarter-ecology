#!/usr/bin/env python3
"""Master figure: EQ betA guild map + leak asymmetry + functional Hill.

Composite Extended Data figure: 4 panels.
  (a) Per-compartment: # bins encoding K00108 per sample
  (b) Per-compartment: K00108 density (per million CDS)
  (c) Leak asymmetry: betA+ vs uptake+ contingency by compartment
  (d) Functional Hill: K00108 density vs Shannon (saturating)
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
import re
def comp_of(sid):
    s = sid.split("_")[-1] if "_" in sid else sid
    m = re.match(r"^[0-9]+([A-Z]+)r[0-9]+$", s)
    if not m: return "?"
    return {"PR":"rhizosphere","S":"surface","D":"deep"}.get(m.group(1), m.group(1))
samp["compartment"] = samp["sample_id"].apply(comp_of)

leak = pd.read_csv(CACHE / "leak_asymmetry_per_bin.tsv", sep="\t")
hillf = pd.read_csv(CACHE / "functional_hill_fit.tsv", sep="\t")

CMAP = {"surface":"#cc6677","rhizosphere":"#117733","deep":"#332288"}
ORDER = ["surface","rhizosphere","deep"]

fig = plt.figure(figsize=(8.6, 5.8))
gs = fig.add_gridspec(2, 2, hspace=0.42, wspace=0.32)
axes = [fig.add_subplot(gs[i,j]) for i in (0,1) for j in (0,1)]

# Panel (a): bins-with-K00108 per sample
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
ax.set_ylabel("# bins with K00108 per sample")
for sp in ("top","right"): ax.spines[sp].set_visible(False)
ax.tick_params(direction="out", length=2.5)
ax.grid(True, alpha=0.25, lw=0.4, axis="y")
ax.set_axisbelow(True)
ax.text(-0.15, 1.07, "a", transform=ax.transAxes, fontweight="bold", fontsize=12)
ax.set_title("Producer-bin count", loc="left", pad=4)

# Panel (b): K00108 density
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
for sp in ("top","right"): ax.spines[sp].set_visible(False)
ax.tick_params(direction="out", length=2.5)
ax.grid(True, alpha=0.25, lw=0.4, axis="y")
ax.set_axisbelow(True)
ax.text(-0.15, 1.07, "b", transform=ax.transAxes, fontweight="bold", fontsize=12)
ax.set_title("Functional density", loc="left", pad=4)

# Panel (c): leak asymmetry per compartment
ax = axes[2]
labels = []; betA_uptake_neg = []; betA_uptake_pos = []
nonbetA_uptake_neg = []; nonbetA_uptake_pos = []
for comp in ORDER:
    sub = leak[leak["compartment"]==comp]
    if len(sub) < 5: continue
    labels.append(comp)
    bp_n = sub[sub["has_betA"]==1]
    bn_n = sub[sub["has_betA"]==0]
    if len(bp_n)>0:
        betA_uptake_pos.append((bp_n["has_uptake_strict"]==1).mean()*100)
        betA_uptake_neg.append((bp_n["has_uptake_strict"]==0).mean()*100)
    else:
        betA_uptake_pos.append(0); betA_uptake_neg.append(0)
    nonbetA_uptake_pos.append((bn_n["has_uptake_strict"]==1).mean()*100)
    nonbetA_uptake_neg.append((bn_n["has_uptake_strict"]==0).mean()*100)
x = np.arange(len(labels))
w = 0.35
ax.bar(x - w/2, betA_uptake_pos, w, color="#117733",
        edgecolor="black", linewidth=0.6, label="betA+ : uptake+")
ax.bar(x + w/2, nonbetA_uptake_pos, w, color="#882255",
        edgecolor="black", linewidth=0.6, label="betA- : uptake+")
ax.set_xticks(x); ax.set_xticklabels(labels, rotation=15)
ax.set_ylabel("% bins with uptake genes")
ax.legend(loc="upper right", frameon=False, fontsize=7)
for sp in ("top","right"): ax.spines[sp].set_visible(False)
ax.tick_params(direction="out", length=2.5)
ax.grid(True, alpha=0.25, lw=0.4, axis="y")
ax.set_axisbelow(True)
ax.text(-0.15, 1.07, "c", transform=ax.transAxes, fontweight="bold", fontsize=12)
ax.set_title("Leak asymmetry", loc="left", pad=4)

# Panel (d): functional Hill — Shannon vs K00108 density
ax = axes[3]
# Recompute merged for plotting
xrf = pd.read_csv(CACHE/"xrf_site_compartment_panel.tsv", sep="\t")
shan = pd.read_csv(CACHE/"per_sample_shannon.tsv", sep="\t")
shan["mg_token"] = shan["sample"].str.split("_").str[1]
samp_x = samp.copy()
import re as _re
def parse_sid(sid):
    s = sid.split("_")[-1] if "_" in sid else sid
    m = _re.match(r"^([0-9]+)([A-Z]+)r([0-9]+)$", s)
    if not m: return None, None, None
    return int(m.group(1)), {"PR":"rhizosphere","S":"surface","D":"deep"}[m.group(2)], int(m.group(3))
samp_x[["site","compartment","replicate"]] = samp_x["sample_id"].apply(
    lambda x: pd.Series(parse_sid(x)))
samp_x = samp_x.merge(shan[["mg_token","shannon"]], left_on="sample_id",
                       right_on="mg_token", how="left")
samp_x = samp_x.merge(xrf[["site","compartment","S"]],
                      on=["site","compartment"], how="left")
m = samp_x.dropna(subset=["S","shannon","K00108_density_per_M_CDS"])

for comp in ORDER:
    sub = m[m["compartment"]==comp]
    ax.scatter(sub["K00108_density_per_M_CDS"].clip(lower=0.01),
               sub["shannon"], s=18, color=CMAP[comp],
               edgecolor="black", linewidth=0.3, alpha=0.75, label=comp)
# Overlay free Hill if we have it
try:
    h3 = hillf[hillf["specification"].str.contains("density", na=False)].iloc[0]
    if not np.isnan(h3.get("K", np.nan)):
        Vmax, K, n = h3["Vmax"], h3["K"], h3["n_Hill"]
        D = np.logspace(-2, 4, 200)
        y = Vmax*(D**n)/(K**n + D**n)
        ax.plot(D, y, color="black", lw=1.4,
                 label=f"Hill: K={K:.1f}, n={n:.1f}")
except Exception:
    pass
ax.set_xscale("log")
ax.set_xlabel("K00108 density (hits per million CDS)")
ax.set_ylabel("Community Shannon")
ax.legend(loc="lower right", frameon=False, fontsize=7)
for sp in ("top","right"): ax.spines[sp].set_visible(False)
ax.tick_params(direction="out", length=2.5)
ax.grid(True, alpha=0.25, lw=0.4)
ax.set_axisbelow(True)
ax.text(-0.15, 1.07, "d", transform=ax.transAxes, fontweight="bold", fontsize=12)
ax.set_title("Functional Hill (B)", loc="left", pad=4)

out_pdf = FIG / "extended_fig_betA_guild_master.pdf"
out_png = FIG / "extended_fig_betA_guild_master.png"
fig.savefig(out_pdf, bbox_inches="tight")
fig.savefig(out_png, bbox_inches="tight", dpi=300)
print(f"wrote {out_pdf}\nwrote {out_png}")
