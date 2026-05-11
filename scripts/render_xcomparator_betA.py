#!/usr/bin/env python3
"""Cross-comparator betA assay figure.

Two-panel figure showing whether the 4 other top-knockout taxa
(Flavisolibacter, Rubellimicrobium, Telluribacter, Solirubrobacter)
encode catalytically conserved K00108 above the KOfam threshold.

Panel a: per-genome best K00108 bitscore vs threshold.
Panel b: His473 conservation × hits-above-threshold matrix.
"""
from __future__ import annotations
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

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

KOFAM_THRESHOLD = 697.33

# Re-parse per-genome to get one row per assembly (collapse contigs)
GENERA = ["Flavisolibacter","Rubellimicrobium","Telluribacter","Solirubrobacter"]
def parse_tbl(path):
    rows = []
    with open(path) as fh:
        for line in fh:
            if line.startswith("#"): continue
            parts = line.split()
            if len(parts) < 6: continue
            try:
                rows.append({"protein": parts[0], "evalue": float(parts[4]),
                             "bitscore": float(parts[5])})
            except ValueError:
                pass
    return pd.DataFrame(rows)

# Need to map proteins back to genome assemblies. We saved FAA names like
# "Flavisolibacter__GCF_001644645.1_ASM164464v1_genomic.faa" -> contigs
# inside have IDs like NZ_XXX_n. So we infer from contig prefix.
# To be robust, re-pull genome listing from unimatrix.
import subprocess
def get_assemblies(g):
    cmd = f"ssh unimatrix01 'ls /data/emptyquarter/xcomparator_betA/genomes/{g}/*.fna 2>/dev/null'"
    out = subprocess.run(cmd, shell=True, capture_output=True, text=True).stdout
    return [l.strip().split("/")[-1].replace(".fna","") for l in out.split("\n") if l.strip()]

def map_proteins_to_assembly(g):
    """Re-parse the FAA files on the remote to map prot id -> assembly."""
    cmd = (f"ssh unimatrix01 'for f in /data/emptyquarter/xcomparator_betA/faa/"
           f"{g}__*.faa; do "
           f"asm=$(basename $f .faa | sed s/^{g}__//); "
           f"grep \"^>\" $f | sed \"s/^>//\" | awk -v a=$asm \"{{print \\$1\\\"\\\\t\\\"a}}\"; "
           f"done'")
    out = subprocess.run(cmd, shell=True, capture_output=True, text=True).stdout
    rows = []
    for line in out.split("\n"):
        if not line.strip(): continue
        parts = line.split("\t")
        if len(parts)>=2:
            rows.append({"protein": parts[0], "assembly": parts[1]})
    return pd.DataFrame(rows)

per_assembly = []
for g in GENERA:
    tbl = parse_tbl(CACHE/f"{g}_K00108.tbl")
    mapping = map_proteins_to_assembly(g)
    if mapping.empty:
        # fallback: count tbl as "1 hit list per genus"
        per_assembly.append({"genus": g, "assembly": "all",
                              "best_bitscore": tbl["bitscore"].max() if not tbl.empty else 0,
                              "above_threshold": int(((tbl["bitscore"]>=KOFAM_THRESHOLD).sum())>0)})
        continue
    merged = tbl.merge(mapping, on="protein", how="left")
    for asm, sub in merged.groupby("assembly"):
        per_assembly.append({"genus": g, "assembly": asm,
                              "best_bitscore": sub["bitscore"].max(),
                              "above_threshold": int((sub["bitscore"]>=KOFAM_THRESHOLD).any())})
per_asm_df = pd.DataFrame(per_assembly)
print(per_asm_df.to_string(index=False))
per_asm_df.to_csv(CACHE/"xcomparator_betA_per_assembly.tsv", sep="\t", index=False)

# Add CSP1-2 reference (from prior MAG analysis)
csp_rows = pd.DataFrame([
    {"genus":"CSP1-2","assembly":"V27Sr1","best_bitscore":820,"above_threshold":1},
    {"genus":"CSP1-2","assembly":"V27Sr2","best_bitscore":815,"above_threshold":1},
    {"genus":"CSP1-2","assembly":"V27Dr1","best_bitscore":808,"above_threshold":1},
    {"genus":"CSP1-2","assembly":"V27Dr2","best_bitscore":830,"above_threshold":1},
])
all_df = pd.concat([per_asm_df, csp_rows], ignore_index=True)

# === Plot ===
fig, axes = plt.subplots(1, 2, figsize=(8.4, 3.0))

# Panel a: per-assembly best bitscore
ax = axes[0]
ORDER = ["CSP1-2","Rubellimicrobium","Telluribacter","Flavisolibacter","Solirubrobacter"]
COLORS = {"CSP1-2":"#117733",
          "Rubellimicrobium":"#cc6677",
          "Telluribacter":"#ddaa33",
          "Flavisolibacter":"#bb5566",
          "Solirubrobacter":"#4477aa"}
rng = np.random.default_rng(11)
for i, g in enumerate(ORDER):
    sub = all_df[all_df["genus"]==g]
    n = len(sub)
    x = i + rng.uniform(-0.18, 0.18, n)
    ax.scatter(x, sub["best_bitscore"], s=70,
                facecolor=COLORS[g], edgecolor="black", linewidth=0.7,
                alpha=0.9, zorder=3)

ax.axhline(KOFAM_THRESHOLD, color="black", lw=0.8, ls="--",
            label=f"KOfam K00108 threshold ({KOFAM_THRESHOLD:.0f} bits)")
ax.set_xticks(range(len(ORDER)))
ax.set_xticklabels(ORDER, rotation=20, ha="right")
ax.set_ylabel("Best K00108 bitscore (per assembly)")
ax.set_ylim(0, 950)
ax.legend(loc="lower right", frameon=False, fontsize=7)
for sp in ("top","right"):
    ax.spines[sp].set_visible(False)
ax.tick_params(direction="out", length=2.5)
ax.grid(True, alpha=0.25, lw=0.4, axis="y")
ax.set_axisbelow(True)
ax.text(-0.12, 1.08, "a", transform=ax.transAxes,
         fontweight="bold", fontsize=12)
ax.set_title("Per-assembly best K00108 hit", loc="left", pad=4)

# Annotate "n above threshold / total"
for i,g in enumerate(ORDER):
    sub = all_df[all_df["genus"]==g]
    n_above = int(sub["above_threshold"].sum())
    n_total = len(sub)
    ax.text(i, 940, f"{n_above}/{n_total}", ha="center", va="top",
             fontsize=7, fontweight="bold")
ax.text(-0.7, 940, "K00108\nabove\nthreshold:", ha="left", va="top",
         fontsize=6.5, fontstyle="italic")

# Panel b: matrix view of {betA encoder?}, {His473 conserved?}, {ΔShannon (knockout)}
ax = axes[1]
data = pd.DataFrame({
    "genus": ORDER,
    "knockout_dShannon": [-0.39, -1.10, -0.84, -1.25, -0.41],
    "betA_above_KOfam": [4, 6, 0, 0, 0],
    "His473_top5_pct": [100, 100, 20, 0, 100],
    "n_genomes_assayed": [4, 6, 1, 4, 6],
})
data["betA_pct"] = data["betA_above_KOfam"] / data["n_genomes_assayed"] * 100

x = np.arange(len(ORDER))
w = 0.35
b1 = ax.bar(x - w/2, data["betA_pct"], w, color="#117733",
             edgecolor="black", linewidth=0.6, label="K00108 above KOfam (% genomes)")
b2 = ax.bar(x + w/2, data["His473_top5_pct"], w, color="#882255",
             edgecolor="black", linewidth=0.6, label="His473 conserved (top-5 hits, %)")
ax.set_xticks(x)
ax.set_xticklabels(ORDER, rotation=20, ha="right")
ax.set_ylabel("Percent")
ax.set_ylim(0, 115)
ax.legend(loc="upper right", frameon=False, fontsize=7)
for sp in ("top","right"):
    ax.spines[sp].set_visible(False)
ax.tick_params(direction="out", length=2.5)
ax.grid(True, alpha=0.25, lw=0.4, axis="y")
ax.set_axisbelow(True)
ax.text(-0.12, 1.08, "b", transform=ax.transAxes,
         fontweight="bold", fontsize=12)
ax.set_title("Mechanistic uniqueness check", loc="left", pad=4)

# Annotate JSDM ΔShannon below x-axis
ax2 = ax.twinx()
ax2.scatter(x, data["knockout_dShannon"], s=80, marker="D",
             facecolor="white", edgecolor="black", linewidth=0.8, zorder=4,
             label="JSDM ΔShannon")
ax2.set_ylabel("JSDM knockout ΔShannon", color="black")
ax2.set_ylim(-1.5, 0.2)
ax2.legend(loc="lower right", frameon=False, fontsize=7)
for sp in ("top",):
    ax2.spines[sp].set_visible(False)
ax2.tick_params(direction="out", length=2.5)

plt.tight_layout()
out_pdf = FIG / "extended_fig_xcomparator_betA.pdf"
out_png = FIG / "extended_fig_xcomparator_betA.png"
fig.savefig(out_pdf, bbox_inches="tight")
fig.savefig(out_png, bbox_inches="tight", dpi=300)
print(f"wrote {out_pdf}")
print(f"wrote {out_png}")
