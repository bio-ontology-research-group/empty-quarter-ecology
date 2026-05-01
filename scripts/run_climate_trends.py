"""
Per-site historical climate trend analysis (1980-2024) for the
60 Empty Quarter sites.

Computes per-site annual time series of:
  - mean annual T
  - annual P (mm)
  - count of >9.5 mm days (surface Hill threshold)
  - count of >13.3 mm days (deep Hill threshold)
  - max daily P (mm/day)
  - count of days >40 degC (extreme-heat days)

Then per-site Mann-Kendall (Kendall-tau) trend test on each metric,
with FDR control (BH) across the 60 sites per metric.

Outputs:
  cache/climate_trends_per_site.tsv
  figures/fig_climate_trends.pdf
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy import stats
from statsmodels.stats.multitest import multipletests

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "src"))
from eq import CACHE_DIR, FIGURES_DIR  # noqa: E402

GEODATA = REPO / "data" / "geodata"
HILL_SURF_MM = 9.5
HILL_DEEP_MM = 13.3
HEAT_THRESH_C = 40.0


def load_daily() -> pd.DataFrame:
    p = CACHE_DIR / "climate_historical_1995_2024.parquet"
    if not p.exists():
        raise SystemExit(
            f"missing {p} — run scripts/fetch_openmeteo_historical.py first"
        )
    df = pd.read_parquet(p)
    df["date"] = pd.to_datetime(df["date"])
    df["year"] = df["date"].dt.year
    return df


def load_site_coords() -> pd.DataFrame:
    sites = []
    for t in (1, 2, 3, 4, 5):
        x = pd.read_csv(GEODATA / f"trip{t}_geodata.tsv", sep="\t")
        x["SiteNum"] = pd.to_numeric(x["Site"], errors="coerce")
        x = x.dropna(subset=["SiteNum"]).drop_duplicates("SiteNum")
        x = x[(x.SiteNum >= 1) & (x.SiteNum <= 60)]
        sites.append(x[["SiteNum", "Latitude", "Longitude"]])
    return (
        pd.concat(sites)
        .groupby("SiteNum")[["Latitude", "Longitude"]]
        .mean()
        .reset_index()
        .rename(columns={"SiteNum": "site",
                         "Latitude": "lat", "Longitude": "lon"})
    )


def annual_metrics(daily: pd.DataFrame) -> pd.DataFrame:
    """Per-site, per-year metrics."""
    g = daily.groupby(["site", "year"])
    annual = g.agg(
        T_mean=("T", "mean"),
        P_total=("P", "sum"),
        P_max=("P", "max"),
        n_pulse_surface=("P", lambda x: int((x > HILL_SURF_MM).sum())),
        n_pulse_deep=("P", lambda x: int((x > HILL_DEEP_MM).sum())),
        n_heat_days=("T", lambda x: int((x > HEAT_THRESH_C).sum())),
    ).reset_index()
    return annual


def mk_trend(annual: pd.DataFrame, metric: str) -> pd.DataFrame:
    """Per-site Mann-Kendall (Kendall tau) trend on `metric`."""
    rows = []
    for site, sub in annual.groupby("site"):
        sub = sub.sort_values("year")
        x = sub.year.values
        y = sub[metric].values
        ok = ~np.isnan(y)
        if ok.sum() < 10:
            continue
        tau, p = stats.kendalltau(x[ok], y[ok])
        # Sen's-slope-style decade rate via OLS (more interpretable than tau)
        slope, intercept, _, _, _ = stats.linregress(x[ok], y[ok])
        rows.append(
            dict(site=int(site), metric=metric,
                 n_years=int(ok.sum()),
                 tau=tau, tau_p=p,
                 ols_slope_per_yr=slope,
                 ols_slope_per_decade=slope * 10.0,
                 baseline_mean=float(np.mean(y[ok])))
        )
    out = pd.DataFrame(rows)
    if not out.empty:
        out["q_BH"] = multipletests(out["tau_p"].values, method="fdr_bh")[1]
        out["credible_q05"] = out["q_BH"] < 0.05
    return out


def main() -> int:
    print("loading daily climate ...")
    daily = load_daily()
    print(f"  rows: {len(daily):,}; sites: {daily.site.nunique()}; "
          f"years: {daily.year.min()}–{daily.year.max()}")
    coords = load_site_coords()

    print("computing annual metrics ...")
    annual = annual_metrics(daily)
    print(f"  annual rows: {len(annual)}")

    metrics = ["T_mean", "P_total", "P_max",
               "n_pulse_surface", "n_pulse_deep", "n_heat_days"]
    trends_all = []
    for m in metrics:
        t = mk_trend(annual, m)
        if not t.empty:
            trends_all.append(t)
    trends = pd.concat(trends_all, ignore_index=True)
    trends = trends.merge(coords, on="site", how="left")
    out_tsv = CACHE_DIR / "climate_trends_per_site.tsv"
    trends.to_csv(out_tsv, sep="\t", index=False)
    print(f"\nwrote {out_tsv}  ({len(trends)} site×metric rows)")

    # ----- Summary across sites per metric -----
    print("\n=== Summary (median per-site OLS trend, per decade) ===")
    summ = (
        trends.groupby("metric")
        .agg(median_decadal=("ols_slope_per_decade", "median"),
             p25=("ols_slope_per_decade", lambda x: np.percentile(x, 25)),
             p75=("ols_slope_per_decade", lambda x: np.percentile(x, 75)),
             n_sites_credibly_increasing=("credible_q05", lambda x: int((x & (trends.loc[x.index, "tau"] > 0)).sum())),
             n_sites_credibly_decreasing=("credible_q05", lambda x: int((x & (trends.loc[x.index, "tau"] < 0)).sum())))
        .reset_index()
    )
    print(summ.to_string(index=False))

    # ----- Figure -----
    print("\nbuilding figure ...")
    fig, axes = plt.subplots(2, 3, figsize=(13.5, 7.5))
    pretty = {
        "T_mean":         ("Mean annual T",            "°C/decade"),
        "P_total":        ("Annual rainfall",          "mm/decade"),
        "P_max":          ("Max daily rainfall",       "mm/decade"),
        "n_pulse_surface":("Pulse days >9.5 mm (surface Hill)", "days/decade"),
        "n_pulse_deep":   ("Pulse days >13.3 mm (deep Hill)",   "days/decade"),
        "n_heat_days":    ("Days >40°C",               "days/decade"),
    }
    flat = axes.flatten()
    from scipy.stats import binomtest
    for ax, m in zip(flat, metrics):
        sub = trends[trends.metric == m].sort_values("lon")
        # bar coloured by raw p<0.05 (per-site significance)
        raw_sig_pos = (sub.tau_p < 0.05) & (sub.tau > 0)
        raw_sig_neg = (sub.tau_p < 0.05) & (sub.tau < 0)
        colors = np.where(raw_sig_pos, "#c0392b",
                          np.where(raw_sig_neg, "#2980b9", "#bdc3c7"))
        ax.bar(np.arange(len(sub)), sub.ols_slope_per_decade,
               color=colors, edgecolor="white", linewidth=0.4, width=0.85)
        ax.axhline(0, color="k", lw=0.4)
        med = sub.ols_slope_per_decade.median()
        ax.axhline(med, color="k", lw=0.8, ls="--",
                   label=f"median = {med:+.2f}")
        ax.set_title(pretty[m][0], fontsize=10)
        ax.set_ylabel(pretty[m][1], fontsize=9)
        ax.set_xlabel("60 sites (sorted by longitude, west $\\rightarrow$ east)",
                      fontsize=8)
        ax.set_xticks([])
        ax.spines["top"].set_visible(False); ax.spines["right"].set_visible(False)
        n_up = int(raw_sig_pos.sum())
        n_dn = int(raw_sig_neg.sum())
        # sign test across sites — how unusual is the directional consistency?
        n_pos = int((sub.tau > 0).sum())
        sign_p = binomtest(n_pos, len(sub), 0.5, alternative="two-sided").pvalue
        ax.text(0.02, 0.97,
                f"{n_up}/{len(sub)} sites $p<0.05$ $\\uparrow$\n"
                f"{n_dn}/{len(sub)} sites $p<0.05$ $\\downarrow$\n"
                f"sign test: {n_pos}/{len(sub)} pos ($p$={sign_p:.1e})",
                transform=ax.transAxes, va="top", ha="left", fontsize=7)
        ax.legend(loc="lower right", frameon=False, fontsize=8)

    fig.suptitle(
        "Per-site climate trends 1995–2024 across the 60 Empty Quarter sites",
        fontsize=12, y=1.00,
    )
    plt.tight_layout(rect=(0, 0, 1, 0.97))
    out = FIGURES_DIR / "fig_climate_trends.pdf"
    plt.savefig(out, bbox_inches="tight")
    print(f"wrote {out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
