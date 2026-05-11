#!/usr/bin/env python3
"""Extended Data Figure: betA biochemistry validation.

Panels:
  (a) Catalytic His473 conservation across genome groups.
  (b) Pfam GMC-oxidoreductase fold presence vs strict K00108 — the
      ``fold-but-no-function'' control showing absence is betA-specific.
  (c) Three-way evidence consistency (KOfam, BLAST-vs-EQ, BLAST-vs-Ec).
  (d) Identity to SwissProt betA references for each EQ CSP1-2 MAG.

Reads:
  data/public_metagenomes/active_site_conservation.tsv
  data/public_metagenomes/pfam_summary.tsv
  data/public_metagenomes/three_way_consistency.tsv
  data/public_metagenomes/eq_vs_swissprot.tsv
"""
from __future__ import annotations

from pathlib import Path
import json

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
DATA = REPO / "data" / "public_metagenomes"
FIG = REPO.parent / "figures"

C_KEYSTONE = "#cc3311"
C_DAD = "#ee7733"
C_DEP = "#117733"
C_SOIL = "#0077bb"
C_REF = "#888888"

GRP_COLOR = {
    "EQ_CSP12":           C_KEYSTONE,
    "A_dadabacteria":     C_DAD,
    "B_dependent_family": C_DEP,
    "C_soil_top20":       C_SOIL,
    "C_soil_control":     C_SOIL,
    "REF_betA":           C_REF,
}
GRP_LABEL = {
    "EQ_CSP12":           "EQ CSP1-2 MAGs",
    "A_dadabacteria":     "Public Dadabacteria",
    "B_dependent_family": "Dependent families",
    "C_soil_top20":       "Top-scoring soil betA",
    "C_soil_control":     "Soil-control genomes",
    "REF_betA":           "Reference SwissProt betA",
}

active = pd.read_csv(DATA / "active_site_conservation.tsv", sep="\t")
pfam = pd.read_csv(DATA / "pfam_summary.tsv", sep="\t")
threeway = pd.read_csv(DATA / "three_way_consistency.tsv", sep="\t")
eqswiss = pd.read_csv(DATA / "eq_vs_swissprot.tsv", sep="\t", header=None,
                       names=["query", "subject", "pid", "alen",
                              "evalue", "bits", "qcov", "scov"])
candmeta = pd.read_csv(DATA / "candidate_betA_meta.tsv", sep="\t")

# ============================================================================
fig = plt.figure(figsize=(7.2, 6.4))
gs = GridSpec(2, 2, figure=fig,
              hspace=0.55, wspace=0.40,
              left=0.10, right=0.97, top=0.93, bottom=0.10)

# ----------------------------------------------------------------------------
# (a) His473 conservation by group
# ----------------------------------------------------------------------------
ax_a = fig.add_subplot(gs[0, 0])
his = (active[active["residue"] == "His473_actsite"]
       .sort_values("group", key=lambda s: s.map(
           {"EQ_CSP12": 0, "A_dadabacteria": 1,
            "B_dependent_family": 2, "C_soil_top20": 3, "REF_betA": 4})))

groups_a = his["group"].tolist()
fracs = (his["n"].values - his["conserved"].values + his["conserved"].values) / his["n"].values
# = 1.0 always; we want fraction conserved
fracs = his["conserved"].values / his["n"].values
ns = his["n"].values
labels_a = [GRP_LABEL[g] for g in groups_a]
colors_a = [GRP_COLOR[g] for g in groups_a]

y = np.arange(len(groups_a))[::-1]
ax_a.barh(y, fracs, color=colors_a, edgecolor="black", linewidth=0.5,
          height=0.65)
for yi, frac, n in zip(y, fracs, ns):
    ax_a.text(frac + 0.015, yi, f"{int(frac*n)} / {n}  ({frac*100:.0f}%)",
              ha="left", va="center", fontsize=6.8)
ax_a.set_yticks(y)
ax_a.set_yticklabels(labels_a, fontsize=7)
ax_a.set_xlabel("Fraction with H at catalytic position 473", fontsize=7.5)
ax_a.set_xlim(0, 1.30)
ax_a.set_xticks([0, 0.25, 0.5, 0.75, 1.0])
ax_a.set_xticklabels(["0", "25%", "50%", "75%", "100%"])
ax_a.tick_params(direction="out", length=2.5)
for sp in ("top", "right"):
    ax_a.spines[sp].set_visible(False)
ax_a.grid(True, axis="x", alpha=0.25, lw=0.4)
ax_a.set_axisbelow(True)

ax_a.text(-0.34, 1.07, "a", transform=ax_a.transAxes,
          fontsize=12, fontweight="bold", va="top", ha="left")
ax_a.set_title("Catalytic His473 conservation",
               loc="left", pad=4, fontsize=8.5)

# ----------------------------------------------------------------------------
# (b) Pfam GMC-fold presence vs strict K00108 — fold-but-no-function
# ----------------------------------------------------------------------------
ax_b = fig.add_subplot(gs[0, 1])
# Combine pfam_summary with K00108 strict counts
strict_counts = {
    "A_dadabacteria":     (0, 16),
    "B_dependent_family": (0, 188),
    "C_soil_control":     (216, 1954),
}
groups_b = ["A_dadabacteria", "B_dependent_family", "C_soil_control"]
x = np.arange(len(groups_b))
width = 0.36
fold_pct = pfam.set_index("selection_group")["pct_any"].reindex(groups_b).values
strict_pct = np.array([100 * strict_counts[g][0] / strict_counts[g][1]
                        for g in groups_b])
labels_b = [GRP_LABEL.get(g, g) for g in groups_b]

ax_b.bar(x - width/2, fold_pct, width=width,
         color="#bbbbbb", edgecolor="black", linewidth=0.5,
         label="Pfam GMC-oxidoreductase fold")
ax_b.bar(x + width/2, strict_pct, width=width,
         color=C_KEYSTONE, edgecolor="black", linewidth=0.5,
         label="Strict K00108 (betA function)")
for xi, fp, sp in zip(x, fold_pct, strict_pct):
    ax_b.text(xi - width/2, fp + 1.5, f"{fp:.1f}%",
              ha="center", va="bottom", fontsize=6.5)
    ax_b.text(xi + width/2, sp + 1.5, f"{sp:.1f}%",
              ha="center", va="bottom", fontsize=6.5,
              color=C_KEYSTONE, fontweight="bold")
ax_b.set_xticks(x)
ax_b.set_xticklabels(labels_b, fontsize=7, rotation=20, ha="right")
ax_b.set_ylabel("% of genomes with detection", fontsize=7.5)
ax_b.set_ylim(0, 80)
ax_b.tick_params(direction="out", length=2.5)
for sp in ("top", "right"):
    ax_b.spines[sp].set_visible(False)
ax_b.legend(loc="upper right", frameon=False, fontsize=6.5,
            handletextpad=0.4, borderpad=0.3, labelspacing=0.3)
ax_b.grid(True, axis="y", alpha=0.25, lw=0.4)
ax_b.set_axisbelow(True)

ax_b.text(-0.18, 1.07, "b", transform=ax_b.transAxes,
          fontsize=12, fontweight="bold", va="top", ha="left")
ax_b.set_title("Fold-but-no-function: GMC fold ≠ betA",
               loc="left", pad=4, fontsize=8.5)

# ----------------------------------------------------------------------------
# (c) Three-way evidence consistency
# ----------------------------------------------------------------------------
ax_c = fig.add_subplot(gs[1, 0])
# For each group, count: KOfam-only, BLAST-only, both, neither
groups_c = ["A_dadabacteria", "B_dependent_family", "C_soil_control"]
cat_data = []
for g in groups_c:
    sub = threeway[threeway["selection_group"] == g]
    if sub.empty:
        continue
    n = len(sub)
    kofam = sub["kofam_trusted"].astype(bool).values
    blast = sub["blast_eq_csp12"].astype(bool).values
    both = (kofam & blast).sum()
    kofam_only = (kofam & ~blast).sum()
    blast_only = (~kofam & blast).sum()
    neither = (~kofam & ~blast).sum()
    cat_data.append({
        "group": g, "n": n,
        "both": both, "kofam_only": kofam_only,
        "blast_only": blast_only, "neither": neither,
    })

cat_df = pd.DataFrame(cat_data)
x = np.arange(len(cat_df))
y_n = cat_df["n"].values
neither_pct = 100 * cat_df["neither"] / y_n
blast_only_pct = 100 * cat_df["blast_only"] / y_n
kofam_only_pct = 100 * cat_df["kofam_only"] / y_n
both_pct = 100 * cat_df["both"] / y_n

ax_c.bar(x, neither_pct, color="#dddddd", edgecolor="black", lw=0.4,
         label="neither")
ax_c.bar(x, blast_only_pct, bottom=neither_pct, color="#88aacc",
         edgecolor="black", lw=0.4, label="BLAST only")
ax_c.bar(x, kofam_only_pct, bottom=neither_pct + blast_only_pct,
         color="#cc8855", edgecolor="black", lw=0.4, label="KOfam only")
ax_c.bar(x, both_pct,
         bottom=neither_pct + blast_only_pct + kofam_only_pct,
         color=C_KEYSTONE, edgecolor="black", lw=0.4, label="both")

for xi, n in zip(x, y_n):
    ax_c.text(xi, 102, f"n={n}", ha="center", va="bottom",
              fontsize=6.5, color="#444")

ax_c.set_xticks(x)
ax_c.set_xticklabels([GRP_LABEL.get(g, g) for g in cat_df["group"]],
                     fontsize=7, rotation=20, ha="right")
ax_c.set_ylabel("% of genomes", fontsize=7.5)
ax_c.set_ylim(0, 110)
ax_c.set_yticks([0, 25, 50, 75, 100])
ax_c.tick_params(direction="out", length=2.5)
for sp in ("top", "right"):
    ax_c.spines[sp].set_visible(False)
ax_c.legend(loc="lower center", bbox_to_anchor=(0.5, -0.45),
            frameon=False, fontsize=6.5,
            handletextpad=0.4, borderpad=0.3, labelspacing=0.3, ncol=4)

ax_c.text(-0.18, 1.07, "c", transform=ax_c.transAxes,
          fontsize=12, fontweight="bold", va="top", ha="left")
ax_c.set_title("KOfam vs BLAST evidence agreement",
               loc="left", pad=4, fontsize=8.5)

# ----------------------------------------------------------------------------
# (d) EQ MAG identity to SwissProt betA references
# ----------------------------------------------------------------------------
ax_d = fig.add_subplot(gs[1, 1])
# For each MAG, take the top 10 swissprot hits and box their %ID
def _short(q: str) -> str:
    parts = q.split("|")
    if len(parts) < 2:
        return q
    short = parts[1].replace("CSP12_", "")
    return short.split("__")[0]

eqswiss["mag_short"] = eqswiss["query"].map(_short)

mags_d = sorted(eqswiss["mag_short"].dropna().unique().tolist())
data_pid = [eqswiss[eqswiss["mag_short"] == m]["pid"].values for m in mags_d]
bp = ax_d.boxplot(data_pid, positions=range(len(mags_d)), widths=0.5,
                  patch_artist=True,
                  boxprops=dict(facecolor=C_KEYSTONE, alpha=0.55,
                                edgecolor="black", linewidth=0.5),
                  medianprops=dict(color="black", lw=1.0),
                  flierprops=dict(marker="o", markersize=2.5,
                                   markerfacecolor="#666",
                                   markeredgecolor="none"))
# Top-1 reference annotated
for i, m in enumerate(mags_d):
    sub = eqswiss[eqswiss["mag_short"] == m].sort_values("pid",
                                                          ascending=False)
    if not sub.empty:
        top = sub.iloc[0]
        ref = top["subject"].split("|")[2].split("_")[1]  # e.g. BRUME
        ax_d.text(i, top["pid"] + 0.5,
                  f"{ref}", ha="center", va="bottom",
                  fontsize=5.8, color="#444", style="italic")

ax_d.axhline(40, color="#888", lw=0.5, ls=(0, (3, 2)))
ax_d.text(len(mags_d) - 0.5, 40.5, "%ID = 40 (typical homology cutoff)",
          ha="right", va="bottom", fontsize=6, color="#888",
          fontstyle="italic")

ax_d.set_xticks(range(len(mags_d)))
ax_d.set_xticklabels(mags_d, fontsize=7, rotation=18, ha="right")
ax_d.set_ylabel("% identity to SwissProt betA", fontsize=7.5)
ax_d.set_ylim(35, 55)
ax_d.tick_params(direction="out", length=2.5)
for sp in ("top", "right"):
    ax_d.spines[sp].set_visible(False)
ax_d.grid(True, axis="y", alpha=0.25, lw=0.4)
ax_d.set_axisbelow(True)

ax_d.text(-0.20, 1.07, "d", transform=ax_d.transAxes,
          fontsize=12, fontweight="bold", va="top", ha="left")
ax_d.set_title("EQ MAG betA vs SwissProt references",
               loc="left", pad=4, fontsize=8.5)

# Save
out_pdf = FIG / "extended_fig_betA_biochemistry.pdf"
out_png = FIG / "extended_fig_betA_biochemistry.png"
fig.savefig(out_pdf, bbox_inches="tight")
fig.savefig(out_png, bbox_inches="tight", dpi=300)
print(f"wrote {out_pdf} ({out_pdf.stat().st_size} bytes)")
print(f"wrote {out_png} ({out_png.stat().st_size} bytes)")
