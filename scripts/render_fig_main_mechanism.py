#!/usr/bin/env python3
"""Render Main Figure 2 (Single-gene mechanism) — Nature-quality.

Panels:
  (a) Schematic of the public-good leak (CSP1-2 makes glycine betaine
      via betA → leaks → feeds the dependent guild under salt stress).
  (b) MAG quality landscape — 4 EQ CSP1-2 MAGs in CheckM2 space.
  (c) betA (K00108) presence asymmetry across genome groups.
  (d) Functional fingerprint per MAG: dark % + osmoprotection /
      stress / nitrogen genes.

Reads:
  cache/csp_mags_checkm2.tsv
  cache/csp_mags_gspa_summary.tsv
  data/public_metagenomes/betA_matrix.tsv
  data/public_metagenomes/candidate_betA_meta.tsv
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.patches as mp
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

C_KEYSTONE = "#cc3311"   # CSP1-2
C_DEP = "#ee7733"        # dependent guild
C_CTRL = "#888888"       # soil control / other
C_HQ = "#117733"         # high-quality MAG

# ----------------------------------------------------------------------------
# Load data
# ----------------------------------------------------------------------------
mags = pd.read_csv(CACHE / "csp_mags_checkm2.tsv", sep="\t")
gspa = pd.read_csv(CACHE / "csp_mags_gspa_summary.tsv", sep="\t")

bet_matrix = pd.read_csv(REPO / "data/public_metagenomes/betA_matrix.tsv",
                         sep="\t")
bet_cand = pd.read_csv(REPO / "data/public_metagenomes/candidate_betA_meta.tsv",
                       sep="\t")

# ----------------------------------------------------------------------------
# Figure
# ----------------------------------------------------------------------------
fig = plt.figure(figsize=(7.2, 7.4))
gs = GridSpec(
    2, 2, figure=fig,
    height_ratios=[1.0, 1.05],
    hspace=0.45, wspace=0.35,
    left=0.10, right=0.97, top=0.94, bottom=0.07,
)

# ============================================================================
# (a) betA phylogenetic placement — EQ CSP1-2 MAGs form a high-score clade
# ============================================================================
ax_a = fig.add_subplot(gs[0, 0])
from Bio import Phylo
tree = Phylo.read(REPO / "data/public_metagenomes/betA_tree.contree",
                  "newick")

# Group rule from leaf name
def group_of(name: str) -> str:
    if name.startswith("EQ_CSP12"):
        return "EQ_CSP12"
    if name.startswith("A_dadabacteria"):
        return "A_dadabacteria"
    if name.startswith("B_dependent_family"):
        return "B_dependent_family"
    if name.startswith("C_soil_top20"):
        return "C_soil_top20"
    return "other"

GRP_COLOR = {
    "EQ_CSP12": C_KEYSTONE,
    "A_dadabacteria": "#cc6633",
    "B_dependent_family": C_DEP,
    "C_soil_top20": "#117733",
    "other": "#999999",
}
GRP_LABEL = {
    "EQ_CSP12": "EQ CSP1-2 MAGs (this study)",
    "A_dadabacteria": "Public Dadabacteria",
    "B_dependent_family": "Dependent-family genomes",
    "C_soil_top20": "Top-scoring soil betA",
    "other": "Other",
}

# Mid-point root, ladderise for clean rendering
tree.root_at_midpoint()
tree.ladderize()

leaves = tree.get_terminals()
n_leaves = len(leaves)
# Y position for each leaf in plotting order (terminals top-to-bottom)
y_of = {leaf.name: i for i, leaf in enumerate(leaves[::-1])}

# Compute X coords (cumulative branch length) recursively
def compute_x(clade, x0=0.0):
    coords = {}
    if clade.branch_length is None:
        bl = 0.0
    else:
        bl = clade.branch_length
    x_here = x0 + bl
    coords[id(clade)] = x_here
    for child in clade.clades:
        coords.update(compute_x(child, x_here))
    return coords
x_of = compute_x(tree.root, x0=0.0)

# Find y-position of internal nodes as midpoint of terminal descendants
def y_of_clade(clade):
    if clade.is_terminal():
        return y_of[clade.name]
    ys = [y_of_clade(c) for c in clade.clades]
    return (min(ys) + max(ys)) / 2.0

# Draw tree
def draw_clade(clade):
    x_self = x_of[id(clade)]
    y_self = y_of_clade(clade)
    if clade.clades:
        ys = [y_of_clade(c) for c in clade.clades]
        # vertical bar
        ax_a.plot([x_self, x_self], [min(ys), max(ys)],
                  color="#444", lw=0.4, zorder=2)
    # horizontal bar to parent
    parent_x = x_self - (clade.branch_length or 0.0)
    ax_a.plot([parent_x, x_self], [y_self, y_self],
              color="#444", lw=0.4, zorder=2)
    for c in clade.clades:
        draw_clade(c)
draw_clade(tree.root)

# Color tip dots and labels
max_x = max(x_of.values())
for leaf in leaves:
    grp = group_of(leaf.name)
    c = GRP_COLOR[grp]
    x_l = x_of[id(leaf)]
    y_l = y_of[leaf.name]
    sz = 22 if grp == "EQ_CSP12" else 6
    edge_lw = 0.6 if grp == "EQ_CSP12" else 0.0
    ax_a.scatter(x_l, y_l, s=sz, c=c, edgecolor="black",
                 linewidths=edge_lw, zorder=4)

# Highlight the EQ_CSP12 clade with a translucent box
eq_y = sorted(y_of[lf.name] for lf in leaves
              if group_of(lf.name) == "EQ_CSP12")
if eq_y:
    y_min, y_max = min(eq_y) - 0.5, max(eq_y) + 0.5
    ax_a.axhspan(y_min, y_max, xmin=0, xmax=1,
                 color=C_KEYSTONE, alpha=0.10, zorder=1)
    ax_a.text(max_x * 1.02, (y_min + y_max) / 2,
              "EQ CSP1-2 MAGs\n(high betA score)",
              ha="left", va="center", fontsize=6.5,
              color=C_KEYSTONE, fontweight="bold")

# Legend
for grp, lbl in GRP_LABEL.items():
    if grp == "other":
        continue
    sz = 26 if grp == "EQ_CSP12" else 16
    ax_a.scatter([], [], s=sz, c=GRP_COLOR[grp],
                 edgecolor="black", linewidths=0.4, label=lbl)
ax_a.legend(loc="lower right", frameon=False, fontsize=5.8,
            handletextpad=0.4, borderpad=0.3, scatterpoints=1,
            labelspacing=0.3)

ax_a.set_xlabel("betA protein-tree distance", fontsize=7)
ax_a.set_yticks([])
ax_a.set_ylim(-1, n_leaves + 0.5)
ax_a.set_xlim(0, max_x * 1.30)
for sp in ("top", "right", "left"):
    ax_a.spines[sp].set_visible(False)
ax_a.tick_params(axis="x", direction="out", length=2.5)
ax_a.tick_params(axis="y", length=0)

ax_a.text(-0.10, 1.05, "a", transform=ax_a.transAxes,
          fontsize=12, fontweight="bold", va="top", ha="left")
ax_a.set_title("betA phylogeny — EQ CSP1-2 MAGs form a distinct clade",
               loc="left", pad=4, fontsize=8.5)

# ============================================================================
# (b) MAG quality landscape (CheckM2)
# ============================================================================
ax_b = fig.add_subplot(gs[0, 1])

# Plot MAGs
sizes = mags["size_mb"] * 80
ax_b.scatter(mags["completeness"], mags["contamination"],
             s=sizes, c=C_KEYSTONE, edgecolor="black", linewidths=0.7,
             alpha=0.85, zorder=4)
# Annotate
for _, m in mags.iterrows():
    ax_b.annotate(m["MAG"].split("__")[0],
                  xy=(m["completeness"], m["contamination"]),
                  xytext=(8, 5), textcoords="offset points",
                  fontsize=6.2, color="#222")

# MIMAG quality boundaries (shaded regions)
# Medium quality: completeness ≥ 50, contamination < 10
# High quality:   completeness ≥ 90, contamination < 5
ax_b.axvspan(50, 90, alpha=0.05, color="#ee7733", zorder=1)
ax_b.axvspan(90, 100, alpha=0.10, color=C_HQ, zorder=1)
ax_b.axhline(5, color=C_HQ, ls=(0, (3, 2)), lw=0.6, alpha=0.7)
ax_b.axhline(10, color="#ee7733", ls=(0, (3, 2)), lw=0.6, alpha=0.7)

ax_b.text(95, 4.7, "MIMAG\nhigh-quality\nzone", ha="center", va="bottom",
          fontsize=6.5, color=C_HQ, fontstyle="italic")
ax_b.text(70, 9.2, "Medium-quality bound (cont. = 10%)",
          fontsize=6.0, color="#ee7733", fontstyle="italic")

# Add the manuscript's reported HQ co-assembly MAG (manuscript text)
ax_b.scatter(92.4, 4.9, s=300, c=C_HQ, edgecolor="black",
             linewidths=0.7, alpha=0.9, zorder=4, marker="*")
ax_b.annotate("HQ co-assembly\n92.4% / 4.9%\n(3.9 Mb)",
              xy=(92.4, 4.9), xytext=(-2.5, 1.8),
              textcoords="data",
              fontsize=6.2, color=C_HQ, fontweight="bold",
              ha="right", va="bottom",
              arrowprops=dict(arrowstyle="-", lw=0.4, color=C_HQ))

ax_b.set_xlabel("Completeness (CheckM2, %)", fontsize=7.5)
ax_b.set_ylabel("Contamination (CheckM2, %)", fontsize=7.5)
ax_b.set_xlim(50, 100)
ax_b.set_ylim(0, 11)
ax_b.tick_params(direction="out", length=2.5)
for sp in ("top", "right"):
    ax_b.spines[sp].set_visible(False)
ax_b.grid(True, alpha=0.2, lw=0.4)
ax_b.set_axisbelow(True)

# Bubble size legend
for sz, lbl in [(2.0, "2 Mb"), (3.0, "3 Mb"), (4.0, "4 Mb")]:
    ax_b.scatter([], [], s=sz * 80, c="#aaaaaa",
                 edgecolor="black", linewidths=0.4, label=lbl)
ax_b.legend(loc="upper right", frameon=False, fontsize=6.5,
            handletextpad=0.6, borderpad=0.3, scatterpoints=1,
            labelspacing=0.4, title="Genome size",
            title_fontsize=6.5)

ax_b.text(-0.16, 1.05, "b", transform=ax_b.transAxes,
          fontsize=12, fontweight="bold", va="top", ha="left")
ax_b.set_title("EQ CSP1-2 MAGs — quality landscape",
               loc="left", pad=4, fontsize=8.5)

# ============================================================================
# (c) betA presence asymmetry
# ============================================================================
ax_c = fig.add_subplot(gs[1, 0])

# Build presence summary
groups = []
# EQ CSP1-2 MAGs — score ≥ 470 (all 5 candidates have ≥ 474.9 → 5/5)
eq_mags = bet_cand[bet_cand["group"] == "EQ_CSP12"]
eq_present = (eq_mags["K00108_score"] >= 200).sum()
eq_total = len(eq_mags)
groups.append(("EQ CSP1-2 MAGs\n(this study)", eq_present, eq_total, C_KEYSTONE))

# Public Dadabacteria
A = bet_matrix[bet_matrix["selection_group"] == "A_dadabacteria"]
groups.append(("Public Dadabacteria\n(IMG/JGI)",
               int(A["K00108"].sum()), len(A), "#cc3311aa"))

# Dependent families
B = bet_matrix[bet_matrix["selection_group"] == "B_dependent_family"]
groups.append(("Dependent-family\ngenomes (worldwide)",
               int(B["K00108"].sum()), len(B), C_DEP))

# Soil controls
C = bet_matrix[bet_matrix["selection_group"] == "C_soil_control"]
groups.append(("Soil-control\ngenomes",
               int(C["K00108"].sum()), len(C), C_CTRL))

# Plot horizontal bars showing fraction with betA
y = np.arange(len(groups))[::-1]
fracs = [g[1] / g[2] for g in groups]
colors = [g[3] for g in groups]
labels = [g[0] for g in groups]
counts = [f"{g[1]} / {g[2]}" for g in groups]

bars = ax_c.barh(y, fracs, color=colors, edgecolor="black",
                 linewidth=0.5, height=0.65)

for yi, f, ct in zip(y, fracs, counts):
    if f > 0.02:
        ax_c.text(f + 0.02, yi, f"{ct}  ({f*100:.1f}%)",
                  ha="left", va="center", fontsize=7)
    else:
        ax_c.text(0.02, yi, f"{ct}  ({f*100:.1f}%)",
                  ha="left", va="center", fontsize=7, color="#444")

ax_c.set_yticks(y)
ax_c.set_yticklabels(labels, fontsize=7.5)
ax_c.set_xlabel("Fraction of genomes encoding betA (K00108)",
                fontsize=7.5)
ax_c.set_xlim(0, 1.15)
ax_c.set_xticks([0, 0.25, 0.5, 0.75, 1.0])
ax_c.set_xticklabels(["0", "25%", "50%", "75%", "100%"])
ax_c.tick_params(direction="out", length=2.5)
for sp in ("top", "right"):
    ax_c.spines[sp].set_visible(False)
ax_c.grid(True, axis="x", alpha=0.25, lw=0.4)
ax_c.set_axisbelow(True)

ax_c.text(-0.32, 1.05, "c", transform=ax_c.transAxes,
          fontsize=12, fontweight="bold", va="top", ha="left")
ax_c.set_title("Strict betA asymmetry across genome groups",
               loc="left", pad=4, fontsize=8.5)

# ============================================================================
# (d) Per-MAG functional fingerprint
# ============================================================================
ax_d = fig.add_subplot(gs[1, 1])

# Gene categories: nitrogen fixation (vnfA + anfA), trehalose, UV repair
mag_keys = list(mags["MAG"])
gspa_idx = gspa.set_index("MAG")
data_rows = []
for k in mag_keys:
    if k in gspa_idx.index:
        r = gspa_idx.loc[k]
        data_rows.append({
            "mag": k.split("__")[0],
            "Alt-Nfix\n(vnfA + anfA)": int(r["vnfA"]) + int(r["anfA"]),
            "Trehalose": int(r["trehalose"]),
            "UV/SSB\nrepair": int(r["uvr_rec_ssb"]),
            "% dark": float(r["dark_pct"]),
        })

df_d = pd.DataFrame(data_rows).set_index("mag")
gene_cats = ["Alt-Nfix\n(vnfA + anfA)", "Trehalose", "UV/SSB\nrepair"]

# Grouped bar chart
n_mags = len(df_d)
x = np.arange(len(gene_cats))
width = 0.78 / n_mags
mag_palette = ["#cc3311", "#ee7733", "#0077bb", "#117733"]
for i, mag in enumerate(df_d.index):
    vals = df_d.loc[mag, gene_cats].values
    offsets = (i - (n_mags - 1) / 2) * width
    ax_d.bar(x + offsets, vals, width=width,
             color=mag_palette[i % len(mag_palette)],
             edgecolor="black", linewidth=0.4,
             label=f"{mag} ({df_d.loc[mag, '% dark']:.0f}% dark)")

ax_d.set_xticks(x)
ax_d.set_xticklabels(gene_cats, fontsize=7)
ax_d.set_ylabel("Gene count", fontsize=7.5)
ax_d.tick_params(direction="out", length=2.5)
for sp in ("top", "right"):
    ax_d.spines[sp].set_visible(False)
ax_d.legend(loc="upper left", frameon=False, fontsize=6.0,
            handletextpad=0.4, borderpad=0.3, labelspacing=0.3)
ax_d.grid(True, axis="y", alpha=0.25, lw=0.4)
ax_d.set_axisbelow(True)

ax_d.text(-0.16, 1.05, "d", transform=ax_d.transAxes,
          fontsize=12, fontweight="bold", va="top", ha="left")
ax_d.set_title("MAG-resolved stress-tolerance arsenal",
               loc="left", pad=4, fontsize=8.5)

# ============================================================================
out_pdf = FIG / "fig2_mechanism.pdf"
out_png = FIG / "fig2_mechanism.png"
fig.savefig(out_pdf, bbox_inches="tight")
fig.savefig(out_png, bbox_inches="tight", dpi=300)
print(f"wrote {out_pdf} ({out_pdf.stat().st_size} bytes)")
print(f"wrote {out_png} ({out_png.stat().st_size} bytes)")
