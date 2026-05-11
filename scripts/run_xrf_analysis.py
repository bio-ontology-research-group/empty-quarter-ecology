#!/usr/bin/env python3
"""Re-run all XRF-dependent descriptive analyses with the all-trips panel.

Outputs (under ``cache/``):
    per_element_shannon.tsv      site-mean per-element XRF x Shannon
                                 (consumed by fig_main1_overview.qmd
                                 panel C and 99_audit.qmd).
    xrf_per_compartment.tsv      per-compartment XRF x Shannon
                                 (consumed by 17_xrf_compartment.qmd
                                 and Suppl S30).
    xrf_summary_all_trips.tsv    coverage summary by trip x compartment.
    xrf_lithology_pca.tsv        coordinates of the all-trips XRF PCA
                                 (replaces the Trip-5-only S4 figure).
    xrf_chemodiversity.tsv       per-sample elemental Shannon entropy.

Figures (under ``figures/``):
    figS30_xrf_per_compartment.pdf  refreshed with all-trip n.
    figS_xrf_pca_alltrips.pdf       PCA of CLR-XRF colored by
                                    compartment, all trips.

Run from the repository root:
    python scripts/run_xrf_analysis.py
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy import stats
from statsmodels.stats.multitest import multipletests
from sklearn.decomposition import PCA

REPO = Path(__file__).resolve().parents[1]
CACHE = REPO / "cache"
FIG = REPO / "figures"
DATA = REPO / "data"
CACHE.mkdir(exist_ok=True)
FIG.mkdir(exist_ok=True)

XRF_PATH = DATA / "geochemistry" / "xrf_lab_table_all_trips.tsv"
FT_PATH = CACHE / "feature_table.parquet"
META_PATH = CACHE / "metadata.parquet"

ELEMENT_COLS = ["Al", "Ba", "Br", "Ca", "Ce", "Cl", "Co", "Cr", "Cs", "Cu",
                "Eu", "Fe", "Ga", "Gd", "K", "La", "Mg", "Mn", "Mo", "Na",
                "Nd", "Ni", "P", "Pr", "S", "Sc", "Si", "Sm", "Sr", "Ti",
                "V", "Zn", "Zr"]
FOCAL = ["S", "Cl", "Na", "P", "Fe", "V", "Mn", "Ca", "K", "Si"]


def shannon_per_sample(ft: pd.DataFrame) -> pd.Series:
    """Per-column (per-sample) Shannon diversity from raw counts."""
    def _h(col: pd.Series) -> float:
        x = col[col > 0].astype(float)
        if x.empty:
            return float("nan")
        p = x / x.sum()
        return float(-(p * np.log(p)).sum())
    return ft.apply(_h, axis=0).rename("Shannon")


def clr(df: pd.DataFrame, eps: float = 1e-6) -> pd.DataFrame:
    """Centred log-ratio along columns (samples)."""
    x = df.replace(0, eps).astype(float)
    log_x = np.log(x)
    return log_x.sub(log_x.mean(axis=0), axis=1)


def load_panel() -> tuple[pd.DataFrame, pd.Series]:
    """Return (xrf, shannon_per_sample) aligned by 16S sample index."""
    xrf = pd.read_csv(XRF_PATH, sep="\t")
    ft = pd.read_parquet(FT_PATH)
    meta = pd.read_parquet(META_PATH)
    sh = shannon_per_sample(ft)
    sh.index.name = "sample"
    meta_sh = meta.join(sh).reset_index()
    return xrf, meta_sh


def site_mean_panel(xrf: pd.DataFrame, meta_sh: pd.DataFrame) -> pd.DataFrame:
    """Aggregate XRF and Shannon to per-site means (across trips and compartments).

    Used for the main-text Fig 1c forest plot.
    """
    xrf_n = (xrf.set_index("SampleID")[ELEMENT_COLS]
             .apply(pd.to_numeric, errors="coerce")
             .fillna(0))
    xrf_n["site"] = xrf.set_index("SampleID")["site"].values
    xrf_site = xrf_n.groupby("site").mean()

    sh_site = (meta_sh.dropna(subset=["site", "Shannon"])
               .groupby("site")["Shannon"].mean())

    common = xrf_site.index.intersection(sh_site.index)
    return xrf_site.loc[common].assign(Shannon=sh_site.loc[common])


def site_compartment_panel(xrf: pd.DataFrame,
                           meta_sh: pd.DataFrame) -> pd.DataFrame:
    """Aggregate XRF and Shannon to per-(site, compartment) means.

    Used for Suppl S30 per-compartment heat-map.
    """
    xrf_long = xrf.copy()
    for c in ELEMENT_COLS:
        xrf_long[c] = pd.to_numeric(xrf_long[c], errors="coerce").fillna(0)
    xrf_long["compartment"] = xrf_long["compartment"].str.lower()
    xrf_sc = (xrf_long.dropna(subset=["site", "compartment"])
              .groupby(["site", "compartment"])[ELEMENT_COLS].mean()
              .reset_index())

    sh_sc = (meta_sh.dropna(subset=["site", "compartment", "Shannon"])
             .groupby(["site", "compartment"])["Shannon"].mean()
             .reset_index())

    return xrf_sc.merge(sh_sc, on=["site", "compartment"], how="inner")


def per_element_shannon(panel: pd.DataFrame) -> pd.DataFrame:
    """Site-mean Spearman ρ for each element vs Shannon, BH-FDR."""
    xrf_clr_T = clr(panel[ELEMENT_COLS].T)  # transpose so samples=cols then take CLR over rows
    # Above call applies CLR across rows (elements) within each sample column.
    # We want CLR per sample across elements → transpose, CLR, transpose back.
    xrf_for_clr = panel[ELEMENT_COLS].T.copy()
    xrf_clr = clr(xrf_for_clr).T  # rows = sites, cols = elements
    rows = []
    for elem in ELEMENT_COLS:
        x = xrf_clr[elem]
        common = x.dropna().index.intersection(panel["Shannon"].dropna().index)
        if len(common) < 10:
            continue
        rho, p = stats.spearmanr(x.loc[common], panel.loc[common, "Shannon"])
        rows.append({"element": elem, "rho": rho, "p": p, "n": len(common)})
    df = pd.DataFrame(rows)
    if df.empty:
        return df
    _, q, _, _ = multipletests(df["p"].fillna(1).values, method="fdr_bh")
    df["q"] = q
    df["significant"] = df["q"] < 0.05
    return df.sort_values("rho", key=abs, ascending=False).reset_index(drop=True)


def per_compartment_corr(panel_sc: pd.DataFrame) -> pd.DataFrame:
    """Per-(element, compartment) Spearman ρ vs Shannon, BH-FDR."""
    rows = []
    for elem in FOCAL:
        for comp in ("surface", "deep", "rhizosphere"):
            sub = panel_sc[panel_sc["compartment"] == comp].dropna(
                subset=[elem, "Shannon"])
            if len(sub) < 8:
                rows.append({"element": elem, "compartment": comp, "rho": np.nan,
                             "p": np.nan, "n": len(sub)})
                continue
            rho, p = stats.spearmanr(sub[elem], sub["Shannon"])
            rows.append({"element": elem, "compartment": comp, "rho": rho,
                         "p": p, "n": len(sub)})
    df = pd.DataFrame(rows)
    valid = df.dropna(subset=["p"])
    _, q, _, _ = multipletests(valid["p"].values, method="fdr_bh")
    df.loc[valid.index, "q"] = q
    df["significant"] = df["q"] < 0.05
    return df


def fig_per_compartment(df: pd.DataFrame, out: Path) -> None:
    pivot_rho = (df.pivot(index="element", columns="compartment", values="rho")
                 .reindex(FOCAL)[["surface", "deep", "rhizosphere"]])
    pivot_q = (df.pivot(index="element", columns="compartment", values="q")
               .reindex(FOCAL)[["surface", "deep", "rhizosphere"]])
    pivot_n = (df.pivot(index="element", columns="compartment", values="n")
               .reindex(FOCAL)[["surface", "deep", "rhizosphere"]])

    fig, ax = plt.subplots(figsize=(4.8, 5.4))
    vmax = float(np.nanmax(np.abs(pivot_rho.values))) or 0.5
    im = ax.imshow(pivot_rho.values, cmap="RdBu_r", vmin=-vmax, vmax=vmax,
                   aspect="auto")
    for i, elem in enumerate(FOCAL):
        for j, comp in enumerate(["surface", "deep", "rhizosphere"]):
            rho = pivot_rho.iloc[i, j]
            q = pivot_q.iloc[i, j]
            n = pivot_n.iloc[i, j]
            if not np.isnan(rho):
                mark = ""
                if not np.isnan(q):
                    if q < 0.001: mark = "***"
                    elif q < 0.01: mark = "**"
                    elif q < 0.05: mark = "*"
                ax.text(j, i, f"{rho:+.2f}{mark}\nn={int(n)}",
                        ha="center", va="center",
                        fontsize=7,
                        color="white" if abs(rho) > 0.45 else "black")
    ax.set_xticks(range(3)); ax.set_xticklabels(["Surf", "Deep", "Rhiz"])
    ax.set_yticks(range(len(FOCAL))); ax.set_yticklabels(FOCAL)
    ax.set_xlabel("Compartment"); ax.set_ylabel("Element")
    ax.set_title("XRF × Shannon (per compartment, all trips)\n"
                 "Spearman ρ; BH-FDR q <0.05 marked *", fontsize=9)
    cb = plt.colorbar(im, ax=ax, shrink=0.7); cb.set_label("Spearman ρ")
    fig.tight_layout()
    fig.savefig(out, bbox_inches="tight")
    plt.close(fig)


def fig_xrf_pca(panel_sc: pd.DataFrame, out: Path) -> tuple[pd.DataFrame, dict]:
    """PCA of CLR-XRF for the all-trip site x compartment panel."""
    X = panel_sc[ELEMENT_COLS]
    Xc = clr(X.T).T  # CLR per sample
    pca = PCA(n_components=4).fit(Xc.values)
    coords = pca.transform(Xc.values)
    out_df = panel_sc[["site", "compartment"]].copy()
    out_df[["PC1", "PC2", "PC3", "PC4"]] = coords[:, :4]
    out_df["Shannon"] = panel_sc["Shannon"].values

    fig, ax = plt.subplots(figsize=(5.2, 4.4))
    palette = {"surface": "#d1410c", "deep": "#1d4ed8", "rhizosphere": "#16a34a"}
    for comp, color in palette.items():
        sub = out_df[out_df["compartment"] == comp]
        ax.scatter(sub["PC1"], sub["PC2"], s=22, c=color, edgecolor="black",
                   linewidth=0.3, alpha=0.7, label=f"{comp} (n={len(sub)})")
    ev = pca.explained_variance_ratio_
    ax.set_xlabel(f"PC1 ({100*ev[0]:.1f}%)")
    ax.set_ylabel(f"PC2 ({100*ev[1]:.1f}%)")
    ax.legend(frameon=False, fontsize=8)
    ax.set_title(f"CLR-XRF PCA, all trips (n={len(out_df)} site×compartment cells)")
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    fig.savefig(out, bbox_inches="tight")
    plt.close(fig)
    return out_df, {"explained_variance_ratio": ev.tolist()}


def chemodiversity(xrf: pd.DataFrame) -> pd.DataFrame:
    """Per-sample elemental Shannon entropy on CLR-positive proportions."""
    rows = []
    for _, r in xrf.iterrows():
        vec = pd.to_numeric(r[ELEMENT_COLS], errors="coerce").fillna(0)
        x = vec[vec > 0]
        if x.empty:
            entropy = float("nan")
        else:
            p = x / x.sum()
            entropy = float(-(p * np.log(p)).sum())
        rows.append({"SampleID": r["SampleID"],
                     "trip": r.get("trip"),
                     "site": r.get("site"),
                     "compartment": r.get("compartment"),
                     "chemodiversity_H": entropy})
    return pd.DataFrame(rows)


def main() -> None:
    print(f"loading {XRF_PATH}")
    xrf, meta_sh = load_panel()
    print(f"  XRF rows: {len(xrf)}  16S samples with Shannon: {meta_sh['Shannon'].notna().sum()}")

    # Coverage summary
    cov = (xrf.dropna(subset=["trip", "compartment"])
           .groupby(["trip", "compartment"]).size().unstack("compartment")
           .fillna(0).astype(int))
    cov.to_csv(CACHE / "xrf_summary_all_trips.tsv", sep="\t")
    print("Per-trip x compartment coverage:")
    print(cov.to_string())

    # 1) site-mean per-element forest panel
    site_panel = site_mean_panel(xrf, meta_sh)
    pe = per_element_shannon(site_panel)
    pe.to_csv(CACHE / "per_element_shannon.tsv", sep="\t", index=False)
    print(f"\nper_element_shannon: {len(pe)} elements; n sites = {site_panel.shape[0]}")
    print(pe.to_string(index=False))

    # 2) per-(site, compartment) panel + per-compartment correlations
    sc = site_compartment_panel(xrf, meta_sh)
    sc.to_csv(CACHE / "xrf_site_compartment_panel.tsv", sep="\t", index=False)
    pc = per_compartment_corr(sc)
    pc.to_csv(CACHE / "xrf_per_compartment.tsv", sep="\t", index=False)
    fig_per_compartment(pc, FIG / "figS30_xrf_per_compartment.pdf")
    print(f"\nxrf_per_compartment: {len(pc)} cells, "
          f"{int(pc['significant'].sum())} q<0.05")
    print(pc[pc["significant"]].to_string(index=False))

    # 3) lithology PCA refresh
    pca_df, pca_meta = fig_xrf_pca(sc, FIG / "figS_xrf_pca_alltrips.pdf")
    pca_df.to_csv(CACHE / "xrf_lithology_pca.tsv", sep="\t", index=False)
    print("\nXRF PCA explained variance (PC1-4):",
          [f"{v*100:.1f}%" for v in pca_meta['explained_variance_ratio']])

    # 4) chemodiversity per sample (used in supp S15)
    chem = chemodiversity(xrf)
    chem.to_csv(CACHE / "xrf_chemodiversity.tsv", sep="\t", index=False)
    by_comp = chem.groupby("compartment")["chemodiversity_H"].agg(["mean", "std", "count"])
    print("\nElemental Shannon entropy by compartment (all trips):")
    print(by_comp.to_string())


if __name__ == "__main__":
    main()
