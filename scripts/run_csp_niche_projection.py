"""
Climate-niche projection of the CSP1-2 keystone signature.

Two-layer model:
  Layer 1 (presence): logit P(CSP1-2 detected at >=85% V4 ID) ~ MAT + MAP
                      fit on the pooled cross-desert + EQ sample-level
                      cohort (n=1,474 across 5 deserts).
  Layer 2 (operational threshold): keystone signature engages only above
                      MAT >= 10 degC.  This encodes the McMurdo-cold
                      drop-out reported in the paper (CSP1-2 still
                      present at 48% but signature ρ -> 0).

The product P_oper = P_present(MAT, MAP) * I(MAT >= 10 degC) is the
"operational keystone" probability mapped onto a global 10-arcmin
grid for current and CMIP6 SSP3-7.0 / SSP2-4.5 (UKESM1-0-LL,
2081–2100) climate using WorldClim 2.1 baseline + WorldClim CMIP6
future bioclim variables.

Outputs:
  cache/niche_model_coeffs.tsv
  cache/niche_grid_summary.tsv      (per-region area changes)
  figures/fig_niche_projection.pdf
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import rasterio
import statsmodels.api as sm

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "src"))
from eq import CACHE_DIR, FIGURES_DIR  # noqa: E402

WC_DIR = REPO / "data" / "worldclim"
# Thermal threshold for signature engagement: calibrated against the
# 5-cohort signature-ρ data — McMurdo (MAT=-21°C, ρ=+0.12 p=0.43, null)
# is the only non-engaging cohort, while Atacama (9°C, ρ=+0.27 p=0.047),
# Gurbantunggut (8°C, ρ=+0.38), Namib (21°C, ρ=+0.64) and EQ (29°C, ρ=+0.37)
# all show positive ρ.  0 °C cleanly separates the engaging from the null
# cohort and is biologically interpretable as the liquid-water boundary.
THERMAL_THRESH_C = 0.0
ARIDITY_MAX_MM_YR = 500.0  # MAP upper bound for the dryland operational envelope
                           # (5 training cohorts span 0–200 mm; 500 is a permissive cap)

# Cohort centroids for samples with no per-sample coordinates
COHORT_CENTROIDS = {
    "Gurbantunggut": (45.0, 87.0),  # Junggar Basin, NW China
}

# ---------------------------------------------------------------------
#  1. Load cross-desert per-sample data + EQ causal-frame
# ---------------------------------------------------------------------
print("loading per-sample cohorts ...")
xd = pd.read_csv(CACHE_DIR / "crossdesert" / "per_sample.tsv", sep="\t")
xd["csp_present"] = (xd["csp_rel_85"] > 0).astype(int)
xd_keep = xd[["desert", "lat", "lon", "csp_present"]].dropna()
print(f"  cross-desert (lat/lon): {len(xd_keep)} rows")

# Gurbantunggut: no lat/lon in per_sample.tsv -> assign cohort centroid
gurb = pd.read_csv(CACHE_DIR / "crossdesert" / "gurbantunggut_per_sample.tsv", sep="\t")
gurb["csp_present"] = (gurb["csp_rel_85"] > 0).astype(int)
gurb["lat"] = COHORT_CENTROIDS["Gurbantunggut"][0]
gurb["lon"] = COHORT_CENTROIDS["Gurbantunggut"][1]
gurb_keep = gurb[["desert", "lat", "lon", "csp_present"]]
print(f"  Gurbantunggut: {len(gurb_keep)} rows (centroid {COHORT_CENTROIDS['Gurbantunggut']})")

# EQ: load sample-level CSP from causal frame, cross-link to site coords
eq_frame = pd.read_parquet(CACHE_DIR / "causal_frame_tier1.parquet")
eq_frame["csp_present"] = (eq_frame["csp_relab"] > 0).astype(int)
GEODATA = REPO / "data" / "geodata"
sites = []
for t in (1, 2, 3, 4, 5):
    df = pd.read_csv(GEODATA / f"trip{t}_geodata.tsv", sep="\t")
    df["SiteNum"] = pd.to_numeric(df["Site"], errors="coerce")
    df = df.dropna(subset=["SiteNum"]).drop_duplicates("SiteNum")
    df = df[(df.SiteNum >= 1) & (df.SiteNum <= 60)]
    sites.append(df[["SiteNum", "Latitude", "Longitude"]])
site_coords = (
    pd.concat(sites)
    .groupby("SiteNum")[["Latitude", "Longitude"]]
    .mean()
    .reset_index()
)
site_coords.columns = ["site", "lat", "lon"]
eq_with_coords = eq_frame[["site", "csp_present"]].merge(site_coords, on="site", how="left")
eq_with_coords["desert"] = "Empty Quarter"
eq_keep = eq_with_coords[["desert", "lat", "lon", "csp_present"]].dropna()
print(f"  Empty Quarter: {len(eq_keep)} rows ({eq_keep.lat.nunique()} site centroids)")

pooled = pd.concat([xd_keep, gurb_keep, eq_keep], ignore_index=True)
print(f"\npooled: {len(pooled)} samples across {pooled.desert.nunique()} deserts")
print(pooled.groupby("desert").agg(n=("csp_present", "size"),
                                   p_present=("csp_present", "mean")).round(3))

# ---------------------------------------------------------------------
#  2. Extract WorldClim MAT/MAP at each sample location
# ---------------------------------------------------------------------
print("\nextracting WorldClim MAT (bio1) and MAP (bio12) ...")
hist_bio1 = WC_DIR / "wc2.1_10m_bio_1.tif"
hist_bio12 = WC_DIR / "wc2.1_10m_bio_12.tif"
assert hist_bio1.exists() and hist_bio12.exists()


def extract_at_points(tif: Path, lats: np.ndarray, lons: np.ndarray) -> np.ndarray:
    with rasterio.open(tif) as src:
        coords = list(zip(lons, lats))
        return np.array([v[0] for v in src.sample(coords)])


pooled["MAT"] = extract_at_points(hist_bio1, pooled.lat.values, pooled.lon.values)
pooled["MAP"] = extract_at_points(hist_bio12, pooled.lat.values, pooled.lon.values)
print(pooled.groupby("desert").agg(MAT_mean=("MAT", "mean"),
                                   MAP_mean=("MAP", "mean")).round(2))

# Drop rows with NoData (pixel falls in nodata e.g. ocean)
pooled = pooled.replace({-3.4e+38: np.nan}).dropna(subset=["MAT", "MAP"]).reset_index(drop=True)
print(f"after NoData drop: {len(pooled)}")

# ---------------------------------------------------------------------
#  3. Fit logistic GLM: logit P(CSP) ~ MAT + MAP
# ---------------------------------------------------------------------
print("\nfitting logistic GLM ...")
X = pd.DataFrame({
    "MAT": pooled.MAT.values,
    "MAP": np.log1p(pooled.MAP.values),  # log-scale rainfall
})
X = sm.add_constant(X)
y = pooled.csp_present.values
glm = sm.Logit(y, X).fit(disp=False)
print(glm.summary())

coeffs = pd.DataFrame({
    "term": ["intercept", "MAT", "log1p(MAP)"],
    "coef": glm.params,
    "se": glm.bse,
    "p": glm.pvalues,
})
coeffs.to_csv(CACHE_DIR / "niche_model_coeffs.tsv", sep="\t", index=False)
print(coeffs)

# ---------------------------------------------------------------------
#  4. Project the model onto current and CMIP6 future global grids
# ---------------------------------------------------------------------
print("\nprojecting onto global 10-arcmin grids ...")


def grid_from_tif(tif: Path, band: int = 1) -> tuple[np.ndarray, rasterio.Affine]:
    with rasterio.open(tif) as src:
        arr = src.read(band).astype(np.float32)
        nodata = src.nodata
        transform = src.transform
    if nodata is not None:
        arr = np.where(np.isclose(arr, nodata), np.nan, arr)
    arr = np.where(arr < -1e30, np.nan, arr)
    return arr, transform


print("  WorldClim historical bio1, bio12")
mat_hist, _ = grid_from_tif(hist_bio1)
map_hist, _ = grid_from_tif(hist_bio12)


def predict_p(mat: np.ndarray, map_: np.ndarray) -> np.ndarray:
    """logistic(intercept + b_MAT * MAT + b_logMAP * log1p(MAP))"""
    a, b_mat, b_map = glm.params.values
    z = a + b_mat * mat + b_map * np.log1p(np.clip(map_, 0, None))
    return 1.0 / (1.0 + np.exp(-z))


def operational_score(mat: np.ndarray, map_: np.ndarray) -> np.ndarray:
    """
    P(operational keystone) = P(present | MAT, MAP)
                              × I(MAT >= THERMAL_THRESH_C)
                              × I(MAP <= ARIDITY_MAX_MM_YR)
    The thermal mask encodes the McMurdo-cold drop-out;
    the aridity mask restricts projection to the dryland envelope
    in which the cohort training data live.
    """
    p = predict_p(mat, map_)
    mask = (mat >= THERMAL_THRESH_C) & (map_ <= ARIDITY_MAX_MM_YR)
    op = np.where(mask, p, 0.0)
    op = np.where(np.isnan(mat) | np.isnan(map_), np.nan, op)
    return op


op_hist = operational_score(mat_hist, map_hist)

print("  WorldClim CMIP6 future bioclim (UKESM1-0-LL ssp3-70 2081–2100, ssp2-45 2081-2100)")
fut_files = {
    "SSP2-4.5_2100": WC_DIR / "wc2.1_10m_bioc_UKESM1-0-LL_ssp245_2081-2100.tif",
    "SSP3-7.0_2100": WC_DIR / "wc2.1_10m_bioc_UKESM1-0-LL_ssp370_2081-2100.tif",
}
op_fut: dict[str, np.ndarray] = {}
for name, tif in fut_files.items():
    with rasterio.open(tif) as src:
        # WorldClim CMIP6 multiband: bio1=band1, bio12=band12
        mat_f = src.read(1).astype(np.float32)
        map_f = src.read(12).astype(np.float32)
        nodata = src.nodata
        if nodata is not None:
            mat_f = np.where(np.isclose(mat_f, nodata), np.nan, mat_f)
            map_f = np.where(np.isclose(map_f, nodata), np.nan, map_f)
        mat_f = np.where(mat_f < -1e30, np.nan, mat_f)
        map_f = np.where(map_f < -1e30, np.nan, map_f)
    op_fut[name] = operational_score(mat_f, map_f)
    print(f"    {name}: MAT future mean = {np.nanmean(mat_f):.2f}°C, "
          f"MAP future mean = {np.nanmean(map_f):.0f} mm")

# Per-region land-area summary (10-arcmin pixel area is a function of latitude)
def pixel_area_km2(arr_shape: tuple[int, int], transform: rasterio.Affine) -> np.ndarray:
    """Approx pixel area in km² for a regular geographic grid."""
    nrows, ncols = arr_shape
    # 10-arcmin = 1/6 deg
    dlat_deg = abs(transform.e)
    dlon_deg = abs(transform.a)
    R = 6371.0  # km
    lat_top = transform.f
    lats = lat_top + dlat_deg * (np.arange(nrows) + 0.5) * np.sign(transform.e)
    # cell area = R^2 * dlon * (sin(lat+dlat/2) - sin(lat-dlat/2))
    cell_area = (R**2 * np.deg2rad(dlon_deg)
                 * (np.sin(np.deg2rad(lats + dlat_deg / 2.0))
                    - np.sin(np.deg2rad(lats - dlat_deg / 2.0))))
    cell_area = np.abs(cell_area)
    return np.broadcast_to(cell_area[:, None], (nrows, ncols))


with rasterio.open(hist_bio1) as src:
    area_grid = pixel_area_km2(src.shape, src.transform)


def total_keystone_area_km2(op: np.ndarray, threshold: float = 0.5) -> float:
    return float(np.nansum(area_grid * (op >= threshold)))


grid_summary = pd.DataFrame(
    [
        {"scenario": "Historical (1970–2000)",
         "p_oper_threshold": 0.5,
         "keystone_area_km2": total_keystone_area_km2(op_hist),
         "fraction_global_land_pct":
             100.0 * total_keystone_area_km2(op_hist)
             / float(np.nansum(area_grid * (~np.isnan(mat_hist))))},
    ]
)
for name, op in op_fut.items():
    a_fut = total_keystone_area_km2(op)
    grid_summary = pd.concat([grid_summary, pd.DataFrame([{
        "scenario": name,
        "p_oper_threshold": 0.5,
        "keystone_area_km2": a_fut,
        "fraction_global_land_pct":
            100.0 * a_fut / float(np.nansum(area_grid * (~np.isnan(mat_hist)))),
    }])], ignore_index=True)
grid_summary["delta_vs_hist_km2"] = (
    grid_summary["keystone_area_km2"] - grid_summary.loc[0, "keystone_area_km2"]
)
grid_summary["delta_vs_hist_pct"] = (
    100.0 * grid_summary["delta_vs_hist_km2"] / grid_summary.loc[0, "keystone_area_km2"]
)
grid_summary.to_csv(CACHE_DIR / "niche_grid_summary.tsv", sep="\t", index=False)
print("\n" + grid_summary.to_string(index=False))

# ---------------------------------------------------------------------
#  5. Figure
# ---------------------------------------------------------------------
print("\nbuilding figure ...")
fig = plt.figure(figsize=(13.5, 9.0))
gs = fig.add_gridspec(2, 3, height_ratios=[1, 1], width_ratios=[1, 1, 1])

cmap = plt.get_cmap("YlOrRd")
norm = mcolors.Normalize(vmin=0.0, vmax=1.0)

with rasterio.open(hist_bio1) as src:
    extent = [src.bounds.left, src.bounds.right, src.bounds.bottom, src.bounds.top]


def draw_map(ax, op_grid, title):
    im = ax.imshow(op_grid, extent=extent, origin="upper", cmap=cmap, norm=norm,
                   interpolation="nearest")
    ax.set_xlim(-180, 180)
    ax.set_ylim(-60, 75)
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_title(title, fontsize=10)
    ax.set_facecolor("#e9eef2")
    return im


ax_hist = fig.add_subplot(gs[0, 0])
im = draw_map(ax_hist, op_hist, "Current (1970–2000)")
ax_2 = fig.add_subplot(gs[0, 1])
draw_map(ax_2, op_fut["SSP2-4.5_2100"], "SSP2-4.5 · 2081–2100")
ax_3 = fig.add_subplot(gs[0, 2])
draw_map(ax_3, op_fut["SSP3-7.0_2100"], "SSP3-7.0 · 2081–2100")

# Add a colorbar in its own bottom-left axis
cbar_ax = fig.add_axes([0.07, 0.50, 0.27, 0.018])
cb = fig.colorbar(im, cax=cbar_ax, orientation="horizontal")
cb.set_label("P(operational CSP1-2 keystone)", fontsize=9)

# Δ-area panel under SSP3-7.0 minus current
ax_delta = fig.add_subplot(gs[1, 0])
delta = op_fut["SSP3-7.0_2100"] - op_hist
delta_im = ax_delta.imshow(
    delta, extent=extent, origin="upper",
    cmap="RdBu_r", vmin=-1.0, vmax=1.0, interpolation="nearest",
)
ax_delta.set_xlim(-180, 180)
ax_delta.set_ylim(-60, 75)
ax_delta.set_xticks([])
ax_delta.set_yticks([])
ax_delta.set_title("Δ P(operational keystone): SSP3-7.0 2100 − Current", fontsize=10)
ax_delta.set_facecolor("#e9eef2")
cbar_dax = fig.add_axes([0.07, 0.06, 0.27, 0.018])
cbd = fig.colorbar(delta_im, cax=cbar_dax, orientation="horizontal")
cbd.set_label("Δ P(operational)", fontsize=9)

# Niche envelope panel
ax_env = fig.add_subplot(gs[1, 1])
mat_grid = np.linspace(-15, 35, 200)
map_grid = np.linspace(0, 1500, 200)
MM, PP = np.meshgrid(mat_grid, map_grid)
op_env = operational_score(MM, PP)
ax_env.contourf(MM, PP, op_env, levels=np.linspace(0, 1, 11), cmap=cmap)
COHORT_COLORS = {
    "Empty Quarter": "#000000", "Atacama": "#984ea3", "Namib": "#1f77b4",
    "Gurbantunggut": "#2ca02c", "McMurdo": "#377eb8",
}
for desert, sub in pooled.groupby("desert"):
    ax_env.scatter(sub.MAT.mean(), sub.MAP.mean(), s=130,
                   color=COHORT_COLORS.get(desert, "k"),
                   edgecolor="white", linewidth=1.3, zorder=5,
                   label=f"{desert} ({100*sub.csp_present.mean():.0f}%)")
ax_env.axvline(THERMAL_THRESH_C, color="grey", linestyle=":", lw=1.0)
ax_env.set_xlabel("MAT (°C)")
ax_env.set_ylabel("MAP (mm/yr)")
ax_env.set_title("Climate-niche envelope · cohort observations", fontsize=10)
ax_env.legend(loc="upper right", frameon=False, fontsize=7.5)

# Area bar
ax_bar = fig.add_subplot(gs[1, 2])
xs = np.arange(len(grid_summary))
ax_bar.bar(xs, grid_summary["keystone_area_km2"] / 1e6,
           color=["#7f7f7f", "#e1a13b", "#c0392b"])
ax_bar.set_xticks(xs)
ax_bar.set_xticklabels(grid_summary["scenario"], rotation=15, fontsize=8, ha="right")
ax_bar.set_ylabel("Operational keystone area (10⁶ km²)")
ax_bar.set_title("Global suitable area", fontsize=10)
for x, v, dv in zip(xs, grid_summary["keystone_area_km2"] / 1e6,
                    grid_summary["delta_vs_hist_pct"]):
    if x == 0:
        ax_bar.text(x, v, f"{v:.1f}", ha="center", va="bottom", fontsize=8)
    else:
        ax_bar.text(x, v, f"{v:.1f}\n({dv:+.0f}%)", ha="center", va="bottom", fontsize=8)
ax_bar.spines["top"].set_visible(False)
ax_bar.spines["right"].set_visible(False)

fig.suptitle(
    "CSP1-2 climate-niche projection: current and CMIP6 SSP futures",
    fontsize=12, y=0.99,
)
plt.tight_layout(rect=(0, 0.02, 1, 0.96))
out = FIGURES_DIR / "fig_niche_projection.pdf"
plt.savefig(out, bbox_inches="tight")
print(f"wrote {out}")
print(f"\nLogistic GLM:")
print(f"  intercept = {glm.params.iloc[0]:+.3f}")
print(f"  beta_MAT  = {glm.params.iloc[1]:+.3f}  (per °C)")
print(f"  beta_logMAP = {glm.params.iloc[2]:+.3f}")
