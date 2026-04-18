#!/usr/bin/env python3
"""Build Fig 11 cross-desert comparison panel."""
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

root = Path(__file__).resolve().parents[1] / "cache" / "crossdesert"
figdir = Path(__file__).resolve().parents[1] / "figures"
figdir.mkdir(exist_ok=True)

per = pd.read_csv(root / "per_sample.tsv", sep="\t")
summary = pd.read_csv(root / "comparison_summary.tsv", sep="\t")
eq_ref = pd.read_csv(root / "eq_shannon_reference.tsv", sep="\t")
eq_ref["desert"] = "EmptyQuarter"
eq_ref = eq_ref.rename(columns={"shannon_asv": "Shannon"})

shannon = pd.concat([
    per[["desert", "shannon_otu97"]].rename(columns={"shannon_otu97": "Shannon"}),
    eq_ref[["desert", "Shannon"]]
], ignore_index=True).dropna()

order = ["EmptyQuarter", "Namib", "McMurdo", "Atacama"]
labels = ["Empty Q", "Namib", "McMurdo", "Atacama"]
colors = ["#444", "#d97706", "#0e7490", "#b91c1c"]

fig, axes = plt.subplots(1, 4, figsize=(13, 3.1))

# Panel a: Shannon distributions
ax = axes[0]
data = [shannon[shannon.desert == d]["Shannon"].dropna().values for d in order]
bp = ax.boxplot(data, patch_artist=True, showfliers=False, widths=0.55)
for patch, c in zip(bp["boxes"], colors):
    patch.set_facecolor(c); patch.set_alpha(0.55); patch.set_edgecolor("k")
for m in bp["medians"]: m.set_color("k")
ax.set_xticks(range(1, 5))
ax.set_xticklabels(labels, rotation=0)
ax.set_ylabel("Shannon (per sample)")
ax.set_title("(a) Diversity distributions")
ax.grid(True, axis="y", alpha=0.3)

# Panel b: CSP1-2 prevalence
ax = axes[1]
prev_map = dict(zip(summary["desert"], summary["frac_samples_with_CSP85"].astype(float)))
vals = [prev_map.get(d, np.nan) for d in order]
ax.bar(range(len(order)), vals, color=colors, edgecolor="k", alpha=0.75)
ax.set_xticks(range(len(order)))
ax.set_xticklabels(labels)
ax.set_ylabel("Fraction of samples with CSP1-2")
ax.set_ylim(0, 1.0)
ax.set_title("(b) CSP1-2 prevalence (85% V4)")
for i, v in enumerate(vals):
    if not np.isnan(v):
        ax.text(i, v + 0.03, f"{v*100:.0f}%", ha="center", fontsize=9)
ax.grid(True, axis="y", alpha=0.3)

# Panel c: Namib CSP1-2 vs Shannon
ax = axes[2]
nam = per[per.desert == "Namib"].copy()
nam["csp_rel_85"] = pd.to_numeric(nam["csp_rel_85"], errors="coerce")
nam = nam.dropna(subset=["csp_rel_85", "shannon_otu97"])
nam_plot = nam[nam.csp_rel_85 > 0]
if len(nam_plot):
    ax.scatter(nam_plot["csp_rel_85"] * 100, nam_plot["shannon_otu97"],
               s=28, alpha=0.75, c="#d97706", edgecolor="k", linewidth=0.5)
    ax.set_xscale("log")
ax.set_xlabel("CSP1-2 relative abundance (%, log)")
ax.set_ylabel("Shannon")
ax.set_title(r"(c) Namib: $\rho=+0.64$, $p=1.4\times10^{-4}$")
ax.grid(True, alpha=0.3)

# Panel d: Atacama EC vs Shannon
ax = axes[3]
atac = per[per.desert == "Atacama"].copy()
atac["ec"] = pd.to_numeric(atac["electrical_conductivity"], errors="coerce")
atac = atac.dropna(subset=["ec", "shannon_otu97"])
atac = atac[atac.ec > 0]
ax.scatter(atac["ec"], atac["shannon_otu97"],
           s=28, alpha=0.75, c="#b91c1c", edgecolor="k", linewidth=0.5)
ax.set_xscale("log")
ax.set_xlabel("EC (dS m$^{-1}$, log)")
ax.set_ylabel("Shannon")
ax.set_title(r"(d) Atacama: $\rho=-0.38$, $p=0.012$")
ax.grid(True, alpha=0.3)

for ax in axes:
    ax.spines[["top", "right"]].set_visible(False)

fig.tight_layout()
out = figdir / "fig11_crossdesert.pdf"
fig.savefig(out, bbox_inches="tight")
fig.savefig(figdir / "fig11_crossdesert.png", bbox_inches="tight", dpi=160)
print(f"wrote {out}")
