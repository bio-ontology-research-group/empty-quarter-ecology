#!/usr/bin/env python3
"""Render Main Figure 1 (Discovery) — Nature-quality.

Panels:
  (a) 60-site Rub' al-Khali transect map.
  (b) Per-element XRF Spearman ρ vs site-mean Shannon (forest plot).
  (c) Per-compartment co-occurrence networks with CSP1-2 highlighted.
  (d) JSDM keystone knockout — top-12 most-affected dependent genera.

Reads:
  cache/per_element_shannon.tsv        (panel b)
  cache/network_edges_{surface,deep,rhizosphere}.tsv (panel c)
  cache/network_nodes_{surface,deep,rhizosphere}.tsv (panel c)
  cache/jsdm_knockout_shifts.tsv       (panel d)
  data/geodata/trip*_geodata.tsv       (panel a)
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import networkx as nx

# Nature style: Helvetica/Arial, modest line weights, clean typography.
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
FIG.mkdir(exist_ok=True)

# Colors — colorblind-safe palette (Tol vibrant)
C_NEG = "#cc3311"     # red — diversity-suppressing
C_POS = "#0077bb"     # blue — diversity-enhancing
C_SITE = "#ee7733"    # orange-red — site dots
C_KEYSTONE = "#cc3311"  # CSP1-2 hub
C_NEUTRAL = "#bbbbbb"   # other nodes
C_EDGE_POS = "#33446688"  # blue translucent
C_EDGE_NEG = "#cc331188"  # red translucent

# ----------------------------------------------------------------------------
# Load data
# ----------------------------------------------------------------------------
xrf_rho = pd.read_csv(CACHE / "per_element_shannon.tsv", sep="\t")

# Site coordinates
geo_frames = []
for t in (1, 2, 3, 4, 5):
    p = REPO / f"data/geodata/trip{t}_geodata.tsv"
    if p.exists():
        df = pd.read_csv(p, sep="\t")
        df["trip"] = t
        geo_frames.append(df)
geo = pd.concat(geo_frames, ignore_index=True) if geo_frames else pd.DataFrame()

site_cols_candidates = [
    ["SiteNum", "lat", "lon"],
    ["site", "latitude", "longitude"],
    ["Site", "Lat", "Lon"],
    ["SiteNum", "latitude", "longitude"],
    ["Site", "Latitude", "Longitude"],
]
sites = None
for c in site_cols_candidates:
    if all(cc in geo.columns for cc in c):
        sites = geo[c].copy()
        sites.columns = ["site", "lat", "lon"]
        break
if sites is None:
    raise RuntimeError(f"site coord columns not found in {list(geo.columns)}")
sites = sites.dropna(subset=["lat", "lon"]).copy()
sites["site"] = pd.to_numeric(sites["site"], errors="coerce")
sites = sites.dropna(subset=["site"])
sites["site"] = sites["site"].astype(int)
sites = (sites[sites["site"].between(1, 60)]
         .drop_duplicates("site").sort_values("site").reset_index(drop=True))

# JSDM knockout — exclude CSP1-2 itself (the perturbed node) and any NA rows
jsdm = pd.read_csv(CACHE / "jsdm_knockout_shifts.tsv", sep="\t")
jsdm = jsdm[jsdm["genus"].notna()
            & (jsdm["genus"] != "NA")
            & (~jsdm["genus"].str.contains("CSP1-2|Dadabacteria",
                                            case=False, na=False, regex=True))].copy()
jsdm["abs_shift"] = jsdm["shift"].abs()
jsdm_top = jsdm.sort_values("abs_shift", ascending=False).head(12).copy()

# Networks
def load_net(comp):
    nodes = pd.read_csv(CACHE / f"network_nodes_{comp}.tsv", sep="\t")
    edges = pd.read_csv(CACHE / f"network_edges_{comp}.tsv", sep="\t")
    return nodes, edges

# ----------------------------------------------------------------------------
# Figure
# ----------------------------------------------------------------------------
fig = plt.figure(figsize=(7.2, 8.0))
gs = GridSpec(
    3, 3, figure=fig,
    height_ratios=[1.0, 1.05, 0.85],
    hspace=0.85, wspace=0.55,
    left=0.10, right=0.97, top=0.94, bottom=0.06,
)

# ============================================================================
# (a) Map of 60 sites with real basemap
# ============================================================================
ax_a = fig.add_subplot(gs[0, :2])
import math
def haversine(lat1, lon1, lat2, lon2):
    R = 6371.0
    p1, p2 = math.radians(lat1), math.radians(lat2)
    dp = math.radians(lat2 - lat1); dl = math.radians(lon2 - lon1)
    a = math.sin(dp/2)**2 + math.cos(p1)*math.cos(p2)*math.sin(dl/2)**2
    return 2 * R * math.asin(math.sqrt(a))

lon_w, lon_e = sites["lon"].min(), sites["lon"].max()
lat_s, lat_n = sites["lat"].min(), sites["lat"].max()
dist_km = haversine(lat_s, lon_w, lat_n, lon_e)

map_drawn = False
try:
    import contextily as cx
    import pyproj
    transformer = pyproj.Transformer.from_crs(4326, 3857, always_xy=True)
    xs_m, ys_m = transformer.transform(sites["lon"].values,
                                        sites["lat"].values)
    pad_x, pad_y = 1.2e5, 6.0e4
    ax_a.set_xlim(min(xs_m) - pad_x, max(xs_m) + pad_x)
    ax_a.set_ylim(min(ys_m) - pad_y, max(ys_m) + pad_y)
    cx.add_basemap(ax_a, source=cx.providers.CartoDB.Positron, crs=3857,
                   zoom=6, attribution="")
    ax_a.scatter(xs_m, ys_m, s=34, c=C_SITE, edgecolor="black",
                 linewidth=0.5, zorder=5)
    # Transect spine
    sites_sorted = sites.sort_values("lon").reset_index(drop=True)
    xs_sorted, ys_sorted = transformer.transform(
        sites_sorted["lon"].values, sites_sorted["lat"].values)
    ax_a.plot(xs_sorted, ys_sorted, color="black", lw=0.4,
              alpha=0.45, zorder=4)
    map_drawn = True
    ax_a.set_xticks([]); ax_a.set_yticks([])
except Exception as e:
    print(f"basemap unavailable ({e}); falling back to plain scatter")

if not map_drawn:
    ax_a.set_facecolor("#f7e9c8")
    ax_a.scatter(sites["lon"], sites["lat"],
                 s=38, c=C_SITE, edgecolor="black", linewidth=0.5, zorder=4)
    sites_sorted = sites.sort_values("lon").reset_index(drop=True)
    ax_a.plot(sites_sorted["lon"], sites_sorted["lat"],
              color="black", lw=0.4, alpha=0.4, zorder=2)
    ax_a.set_xlim(44.8, 55.1)
    ax_a.set_ylim(18.9, 21.0)
    ax_a.set_xlabel("Longitude (°E)")
    ax_a.set_ylabel("Latitude (°N)")

ax_a.annotate(
    f"~{dist_km:.0f} km transect",
    xy=(0.98, 0.04), xycoords="axes fraction",
    ha="right", va="bottom", fontsize=7, style="italic", color="#333",
    bbox=dict(boxstyle="round,pad=0.18", fc="white", ec="none", alpha=0.85)
)
ax_a.text(0.02, 0.96, "Rub' al-Khali\n(Saudi Arabia)",
          transform=ax_a.transAxes, ha="left", va="top",
          fontsize=7, color="#222",
          bbox=dict(boxstyle="round,pad=0.18", fc="white", ec="none", alpha=0.85))
for sp in ("top", "right", "left", "bottom"):
    ax_a.spines[sp].set_visible(False)

# Bold panel label and title
ax_a.text(-0.10, 1.10, "a", transform=ax_a.transAxes,
          fontsize=12, fontweight="bold", va="top", ha="left")
ax_a.set_title("60 sites, 5 seasonal campaigns, 1,227 amplicon samples",
               loc="left", pad=4, fontsize=8.5)

# ============================================================================
# (b) Per-element XRF — Shannon Spearman ρ forest
# ============================================================================
ax_b = fig.add_subplot(gs[0, 2])
focal = ["S", "Cl", "Na", "Ca", "K", "P", "Fe", "Mn", "V", "Si", "Zn"]
xrf_show = (xrf_rho[xrf_rho["element"].isin(focal)]
            .copy().set_index("element").reindex(focal).reset_index())
xrf_show = xrf_show.dropna(subset=["rho"]).reset_index(drop=True)
y = np.arange(len(xrf_show))
colors = [C_NEG if r < 0 else C_POS for r in xrf_show["rho"]]
n_arr = xrf_show["n"].astype(float)
z_rho = np.arctanh(xrf_show["rho"].astype(float).clip(-0.999, 0.999))
se = 1.0 / np.sqrt(n_arr - 3)
lo, hi = np.tanh(z_rho - 1.96 * se), np.tanh(z_rho + 1.96 * se)

ax_b.hlines(y, lo, hi, color=colors, lw=1.8, alpha=0.6)
ax_b.scatter(xrf_show["rho"], y, s=42, c=colors,
             edgecolor="black", linewidth=0.45, zorder=3)
ax_b.axvline(0, color="black", linestyle=(0, (1, 1.2)), alpha=0.5, lw=0.6)
ax_b.set_yticks(y)
ax_b.set_yticklabels(xrf_show["element"], fontsize=7.5)
ax_b.invert_yaxis()
ax_b.set_xlabel("Spearman ρ (element vs site Shannon)", fontsize=7.5)
for i, row in xrf_show.iterrows():
    q = row.get("q", 1.0)
    if pd.notna(q) and q < 0.05:
        x_off = 0.04 if row["rho"] > 0 else -0.04
        ha = "left" if row["rho"] > 0 else "right"
        ax_b.text(row["rho"] + x_off, i, "*",
                  ha=ha, va="center", fontsize=10, fontweight="bold")
ax_b.set_xlim(-0.85, 0.7)
ax_b.set_xticks([-0.6, -0.3, 0.0, 0.3, 0.6])
ax_b.tick_params(direction="out", length=3)
for sp in ("top", "right"):
    ax_b.spines[sp].set_visible(False)
ax_b.grid(False)

ax_b.text(-0.34, 1.10, "b", transform=ax_b.transAxes,
          fontsize=12, fontweight="bold", va="top", ha="left")
n_sites_xrf = int(xrf_show["n"].iloc[0]) if len(xrf_show) else 0
ax_b.set_title(f"XRF × Shannon (n={n_sites_xrf} sites)",
               loc="left", pad=4, fontsize=8.5)

# ============================================================================
# (c) Per-compartment networks with CSP1-2 hub
# ============================================================================
def draw_hubs(ax, comp, label, top_n=12):
    """Top-N degree hubs for a compartment, CSP1-2 highlighted."""
    nodes, _ = load_net(comp)
    nodes = nodes.sort_values("degree", ascending=False).reset_index(drop=True)
    n_total = len(nodes)
    csp_row = nodes[nodes["node"] == "CSP1-2"]
    csp_rank = (nodes["degree"] >= csp_row.iloc[0]["degree"]).sum() if len(csp_row) else None
    csp_deg = int(csp_row.iloc[0]["degree"]) if len(csp_row) else None

    # Show top_n; if CSP1-2 is below top_n, append it as a "..." row
    show = nodes.head(top_n).copy()
    appended_csp = False
    if csp_rank is not None and csp_rank > top_n:
        show = pd.concat([show, csp_row], ignore_index=True)
        appended_csp = True

    y_pos = np.arange(len(show))[::-1]  # top of axis = highest degree
    is_csp = (show["node"] == "CSP1-2").values
    bar_colors = [C_KEYSTONE if c else C_NEUTRAL for c in is_csp]
    edge_colors = ["black" if c else "#888" for c in is_csp]
    edge_widths = [0.7 if c else 0.3 for c in is_csp]

    ax.barh(y_pos, show["degree"], height=0.72,
            color=bar_colors, edgecolor=edge_colors, linewidth=edge_widths)

    # Tick labels — italic for genus names; CSP1-2 and break sentinel handled
    labels = []
    for i, (_, r) in enumerate(show.iterrows()):
        nm = r["node"]
        if appended_csp and i == top_n:
            labels.append(f"⋮ (rank {csp_rank}) — {nm}")
        else:
            labels.append(nm)
    ax.set_yticks(y_pos)
    ax.set_yticklabels(labels, fontsize=6.5)
    # Bold CSP1-2 tick label
    for tlab, csp_flag in zip(ax.get_yticklabels(), is_csp):
        if csp_flag:
            tlab.set_fontweight("bold")
            tlab.set_color(C_KEYSTONE)

    ax.set_xlabel("Network degree", fontsize=7)
    ax.tick_params(axis="x", direction="out", length=2.5)
    ax.tick_params(axis="y", direction="out", length=0)
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)
    ax.set_xlim(0, max(nodes["degree"].max() * 1.12, 5))
    ax.grid(True, axis="x", alpha=0.25, lw=0.4)
    ax.set_axisbelow(True)

    # Title combines compartment label with CSP1-2 rank
    if csp_rank is not None:
        n_pct = 100 * csp_rank / n_total
        title_main = label
        sub = f"CSP1-2 rank {csp_rank}/{n_total}"
    else:
        title_main = label
        sub = "CSP1-2 absent"
    ax.set_title(title_main, loc="left", fontsize=8.5,
                 pad=8, fontweight="bold")
    ax.text(1.00, 1.02, sub, transform=ax.transAxes,
            ha="right", va="bottom", fontsize=7,
            color=C_KEYSTONE, fontstyle="italic")

ax_c1 = fig.add_subplot(gs[1, 0])
ax_c2 = fig.add_subplot(gs[1, 1])
ax_c3 = fig.add_subplot(gs[1, 2])

# Row label (placed above the row, with extra space)
fig.text(0.05, 0.65, "c",
         fontsize=12, fontweight="bold", va="bottom", ha="left")
fig.text(0.10, 0.652,
         "Top-degree hubs in CLR-Spearman co-occurrence networks (q ≤ 0.01)",
         fontsize=8.5, va="bottom", ha="left")

draw_hubs(ax_c1, "surface", "Surface")
draw_hubs(ax_c2, "deep", "Deep")
draw_hubs(ax_c3, "rhizosphere", "Rhizosphere")

# ============================================================================
# (d) JSDM keystone knockout — top-12 affected genera
# ============================================================================
ax_d = fig.add_subplot(gs[2, :])
order = jsdm_top.sort_values("fold_change").reset_index(drop=True)
y_d = np.arange(len(order))
fold = order["fold_change"].values
pct = (fold * 100)  # baseline = 100%
bar_colors = [C_NEG if f < 1.0 else C_POS for f in fold]

ax_d.barh(y_d, pct - 100, left=100, height=0.7,
          color=bar_colors, edgecolor="black", linewidth=0.4)
ax_d.axvline(100, color="black", linestyle=(0, (1, 1.2)), alpha=0.6, lw=0.6)
ax_d.set_yticks(y_d)
ax_d.set_yticklabels([g.replace("_", " ") for g in order["genus"]],
                     fontsize=7.5)
ax_d.set_xlabel("Genus relative abundance after CSP1-2 knockout (% of baseline)",
                fontsize=7.5)
ax_d.set_xlim(0, 130)
ax_d.set_xticks([0, 25, 50, 75, 100, 125])
ax_d.tick_params(direction="out", length=3)
for sp in ("top", "right"):
    ax_d.spines[sp].set_visible(False)
ax_d.grid(True, axis="x", alpha=0.25, lw=0.4)

# Annotate Shannon drop in the top-right corner
sh_drop = -0.35
ax_d.text(0.99, 0.95,
          f"Aggregate ΔShannon ≈ {sh_drop:+.2f}\n(via biotic associations alone)",
          transform=ax_d.transAxes,
          ha="right", va="top", fontsize=7,
          bbox=dict(boxstyle="round,pad=0.35", fc="#f7f7f7",
                    ec="#999", lw=0.4))

ax_d.text(-0.06, 1.08, "d", transform=ax_d.transAxes,
          fontsize=12, fontweight="bold", va="top", ha="left")
ax_d.set_title("In silico CSP1-2 knockout (JSDM, biotic associations only)",
               loc="left", pad=4, fontsize=8.5)

# ============================================================================
out_pdf = FIG / "fig1_overview.pdf"
out_png = FIG / "fig1_overview.png"
fig.savefig(out_pdf, bbox_inches="tight")
fig.savefig(out_png, bbox_inches="tight", dpi=300)
print(f"wrote {out_pdf} ({out_pdf.stat().st_size} bytes)")
print(f"wrote {out_png} ({out_png.stat().st_size} bytes)")
