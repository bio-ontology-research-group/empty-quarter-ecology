#!/usr/bin/env python3
"""Re-render the XRF-dependent supplementary figures with all-trip data.

Refreshes:
    figures/supp_RQ10_Lithology_fig1_xrf_pca_compartment.pdf
    figures/supp_RQ10_Lithology_fig2_si_ca_vs_shannon.pdf
    figures/supp_RQ10_Lithology_fig3_element_diversity_heatmap.pdf
    figures/supp_RQ16_Chemodiversity_fig1_chemodiv_vs_shannon.pdf
    figures/supp_RQ11_DesertVarnish_fig1_mn_vs_radiation_resistant.pdf
    figures/supp_RQ12_SulfurCycle_fig1_s_vs_srb_sob.pdf

Also writes a summary text file with the new statistics
(``cache/supplement_xrf_stats.txt``) so the supplement prose can be
updated.

Inputs:
    data/geochemistry/xrf_lab_table_all_trips.tsv
    cache/feature_table.parquet, cache/taxonomy.parquet,
    cache/metadata.parquet
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy import stats
from sklearn.decomposition import PCA
from sklearn.cluster import KMeans
from statsmodels.stats.multitest import multipletests
from skbio.stats.distance import permanova, DistanceMatrix
from scipy.spatial.distance import pdist, squareform

REPO = Path(__file__).resolve().parents[1]
CACHE = REPO / "cache"
FIG = REPO / "figures"
DATA = REPO / "data"
FIG.mkdir(exist_ok=True)

ELEMENT_COLS = ["Al", "Ba", "Br", "Ca", "Ce", "Cl", "Co", "Cr", "Cs", "Cu",
                "Eu", "Fe", "Ga", "Gd", "K", "La", "Mg", "Mn", "Mo", "Na",
                "Nd", "Ni", "P", "Pr", "S", "Sc", "Si", "Sm", "Sr", "Ti",
                "V", "Zn", "Zr"]
RAD_RESISTANT = ["Rubrobacter", "Deinococcus", "Geodermatophilus",
                 "Modestobacter", "Blastococcus"]
SRB_GENERA = ["Desulfovibrio", "Desulfotomaculum", "Desulfobacter",
              "Desulfobulbus", "Desulfomicrobium", "Desulfosporosinus",
              "Desulfomonile", "Desulfococcus"]
SOB_GENERA = ["Thiobacillus", "Acidithiobacillus", "Sulfurimonas",
              "Sulfurovum", "Beggiatoa", "Thiothrix", "Thiomicrospira"]

stats_buf: list[str] = []


def log(msg: str) -> None:
    print(msg)
    stats_buf.append(msg)


def shannon(col: pd.Series) -> float:
    x = col[col > 0].astype(float)
    if x.empty:
        return float("nan")
    p = x / x.sum()
    return float(-(p * np.log(p)).sum())


def clr_per_sample(df: pd.DataFrame, eps: float = 1e-6) -> pd.DataFrame:
    x = df.replace(0, eps).astype(float)
    log_x = np.log(x)
    return log_x.sub(log_x.mean(axis=1), axis=0)


def load_inputs():
    xrf = pd.read_csv(DATA / "geochemistry" / "xrf_lab_table_all_trips.tsv", sep="\t")
    xrf["compartment"] = xrf["compartment"].str.lower()
    ft = pd.read_parquet(CACHE / "feature_table.parquet")
    tax = pd.read_parquet(CACHE / "taxonomy.parquet")
    meta = pd.read_parquet(CACHE / "metadata.parquet")
    sh = ft.apply(shannon, axis=0).rename("Shannon")
    sh.index.name = "sample"
    meta_sh = meta.join(sh).reset_index()
    return xrf, ft, tax, meta, meta_sh


def site_compartment_panel(xrf, meta_sh):
    xrf_sc = (xrf.dropna(subset=["site", "compartment"])
              .groupby(["site", "compartment"])[ELEMENT_COLS].mean()
              .reset_index())
    sh_sc = (meta_sh.dropna(subset=["site", "compartment", "Shannon"])
             .groupby(["site", "compartment"])["Shannon"].mean()
             .reset_index())
    return xrf_sc.merge(sh_sc, on=["site", "compartment"], how="inner")


def lithology_pca(panel_sc):
    """S4 fig 1: PCA on CLR-XRF colored by compartment, with K-means clusters."""
    Xc = clr_per_sample(panel_sc[ELEMENT_COLS])
    pca = PCA(n_components=4).fit(Xc.values)
    coords = pca.transform(Xc.values)
    panel_sc = panel_sc.copy()
    panel_sc[["PC1", "PC2"]] = coords[:, :2]

    km = KMeans(n_clusters=3, n_init=10, random_state=0).fit(coords[:, :3])
    panel_sc["xrf_cluster"] = km.labels_

    # PERMANOVA: cluster + compartment vs CLR-XRF distance
    dm = DistanceMatrix(squareform(pdist(Xc.values, "euclidean")),
                         ids=[str(i) for i in range(len(panel_sc))])
    perm_cluster = permanova(dm, panel_sc["xrf_cluster"].astype(str).values,
                             permutations=999)
    perm_comp = permanova(dm, panel_sc["compartment"].values, permutations=999)

    # Si:Ca ratio
    sc = panel_sc.copy()
    sc["si_ca"] = sc["Si"] / sc["Ca"].replace(0, np.nan)
    si_ca_clean = sc["si_ca"].dropna()
    log(f"\n[S4 Lithology] n = {len(panel_sc)} site x compartment cells "
        f"(was 64 in Trip-5-only).")
    log(f"  PCA variance: PC1 = {100*pca.explained_variance_ratio_[0]:.1f}%, "
        f"PC2 = {100*pca.explained_variance_ratio_[1]:.1f}%")
    log(f"  Si:Ca: range {si_ca_clean.min():.2f} to {si_ca_clean.max():.2f}, "
        f"mean {si_ca_clean.mean():.2f}; "
        f"siliceous (Si:Ca > 1) = {int((si_ca_clean > 1).sum())}/{len(si_ca_clean)}")
    log(f"  PERMANOVA cluster: F = {perm_cluster['test statistic']:.2f}, "
        f"p = {perm_cluster['p-value']:.3f}")
    log(f"  PERMANOVA compartment: F = {perm_comp['test statistic']:.2f}, "
        f"p = {perm_comp['p-value']:.3f}")

    fig, ax = plt.subplots(figsize=(6, 4.5))
    palette = {"surface": "#d1410c", "deep": "#1d4ed8", "rhizosphere": "#16a34a"}
    for comp, color in palette.items():
        sub = panel_sc[panel_sc["compartment"] == comp]
        ax.scatter(sub["PC1"], sub["PC2"], s=22, c=color, edgecolor="black",
                   linewidth=0.3, alpha=0.7, label=f"{comp} (n={len(sub)})")
    ev = pca.explained_variance_ratio_
    ax.set_xlabel(f"PC1 ({100*ev[0]:.1f}%)")
    ax.set_ylabel(f"PC2 ({100*ev[1]:.1f}%)")
    ax.legend(frameon=False, fontsize=8)
    ax.set_title(f"CLR-XRF PCA, all trips (n={len(panel_sc)} site×compartment)")
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    out = FIG / "supp_RQ10_Lithology_fig1_xrf_pca_compartment.pdf"
    fig.savefig(out, bbox_inches="tight")
    plt.close(fig)
    return panel_sc, {"pc1": ev[0], "pc2": ev[1],
                      "perm_cluster_F": perm_cluster["test statistic"],
                      "perm_cluster_p": perm_cluster["p-value"],
                      "perm_comp_F": perm_comp["test statistic"],
                      "perm_comp_p": perm_comp["p-value"],
                      "n_pairs": len(panel_sc),
                      "siliceous_n": int((si_ca_clean > 1).sum()),
                      "si_ca_range": (float(si_ca_clean.min()), float(si_ca_clean.max())),
                      "si_ca_mean": float(si_ca_clean.mean())}


def si_ca_vs_shannon(panel_sc):
    """S4 fig 2: Si:Ca ratio vs Shannon by compartment."""
    sc = panel_sc.copy()
    sc["si_ca"] = sc["Si"] / sc["Ca"].replace(0, np.nan)
    sc = sc.dropna(subset=["si_ca", "Shannon"])

    fig, axes = plt.subplots(1, 3, figsize=(11, 3.5), sharey=True)
    palette = {"surface": "#d1410c", "deep": "#1d4ed8", "rhizosphere": "#16a34a"}
    rows = []
    for ax, comp in zip(axes, ["surface", "deep", "rhizosphere"]):
        sub = sc[sc["compartment"] == comp]
        ax.scatter(sub["si_ca"], sub["Shannon"], s=22, c=palette[comp],
                   edgecolor="black", linewidth=0.3, alpha=0.75)
        if len(sub) >= 8:
            rho, p = stats.spearmanr(sub["si_ca"], sub["Shannon"])
        else:
            rho, p = float("nan"), float("nan")
        ax.set_title(f"{comp}: ρ = {rho:+.2f}, p = {p:.3f}, n = {len(sub)}")
        ax.set_xlabel("Si : Ca")
        ax.set_xscale("log")
        ax.spines[["top", "right"]].set_visible(False)
        rows.append({"compartment": comp, "rho": rho, "p": p, "n": len(sub)})
    axes[0].set_ylabel("Shannon")
    fig.tight_layout()
    out = FIG / "supp_RQ10_Lithology_fig2_si_ca_vs_shannon.pdf"
    fig.savefig(out, bbox_inches="tight")
    plt.close(fig)
    log("\n[S4 Si:Ca vs Shannon, all trips]")
    for r in rows:
        log(f"  {r['compartment']:11s} rho = {r['rho']:+.3f}, p = {r['p']:.3g}, n = {r['n']}")
    return pd.DataFrame(rows)


def element_diversity_heatmap(panel_sc):
    """S4 fig 3: per-element x compartment Spearman heat-map."""
    rows = []
    for elem in ELEMENT_COLS:
        for comp in ("surface", "deep", "rhizosphere"):
            sub = panel_sc[panel_sc["compartment"] == comp].dropna(subset=[elem, "Shannon"])
            if len(sub) < 8 or sub[elem].nunique() < 3:
                rho, p = float("nan"), float("nan")
            else:
                rho, p = stats.spearmanr(sub[elem], sub["Shannon"])
            rows.append({"element": elem, "compartment": comp, "rho": rho, "p": p, "n": len(sub)})
    df = pd.DataFrame(rows)
    valid = df.dropna(subset=["p"])
    _, q, _, _ = multipletests(valid["p"].values, method="fdr_bh")
    df.loc[valid.index, "q"] = q

    pivot_rho = (df.pivot(index="element", columns="compartment", values="rho")
                 .reindex(ELEMENT_COLS)[["surface", "deep", "rhizosphere"]])
    pivot_q = (df.pivot(index="element", columns="compartment", values="q")
               .reindex(ELEMENT_COLS)[["surface", "deep", "rhizosphere"]])

    fig, ax = plt.subplots(figsize=(5, 8))
    vmax = float(np.nanmax(np.abs(pivot_rho.values))) or 0.5
    im = ax.imshow(pivot_rho.values, cmap="RdBu_r", vmin=-vmax, vmax=vmax,
                   aspect="auto")
    for i, elem in enumerate(ELEMENT_COLS):
        for j, comp in enumerate(["surface", "deep", "rhizosphere"]):
            rho = pivot_rho.iloc[i, j]
            q = pivot_q.iloc[i, j]
            if not np.isnan(rho):
                mark = ""
                if not np.isnan(q):
                    if q < 0.001: mark = "***"
                    elif q < 0.01: mark = "**"
                    elif q < 0.05: mark = "*"
                ax.text(j, i, f"{rho:+.2f}{mark}", ha="center", va="center",
                        fontsize=7,
                        color="white" if abs(rho) > 0.45 else "black")
    ax.set_xticks(range(3)); ax.set_xticklabels(["Surf", "Deep", "Rhiz"])
    ax.set_yticks(range(len(ELEMENT_COLS))); ax.set_yticklabels(ELEMENT_COLS)
    ax.set_xlabel("Compartment"); ax.set_ylabel("Element")
    ax.set_title("Element-diversity correlations (all trips)\n"
                 "Spearman ρ; * = BH-FDR q < 0.05",
                 fontsize=9)
    cb = plt.colorbar(im, ax=ax, shrink=0.5); cb.set_label("Spearman ρ")
    fig.tight_layout()
    out = FIG / "supp_RQ10_Lithology_fig3_element_diversity_heatmap.pdf"
    fig.savefig(out, bbox_inches="tight")
    plt.close(fig)

    n_sig = int((df["q"] < 0.05).sum())
    log(f"\n[S4 Element-diversity heatmap]")
    log(f"  {n_sig} of {len(df)} (element x compartment) cells significant at FDR q<0.05")
    return df


def chemodiversity_vs_shannon(xrf, meta_sh):
    """S15: per-sample chemodiversity vs taxonomic Shannon by compartment."""
    chem_rows = []
    for _, r in xrf.iterrows():
        vec = pd.to_numeric(r[ELEMENT_COLS], errors="coerce").fillna(0)
        x = vec[vec > 0]
        if x.empty:
            entropy = float("nan")
        else:
            p = x / x.sum()
            entropy = float(-(p * np.log(p)).sum())
        chem_rows.append({"trip": r.get("trip"), "site": r.get("site"),
                          "compartment": r.get("compartment"),
                          "chemodiversity_H": entropy})
    chem = pd.DataFrame(chem_rows)
    chem_cell = (chem.dropna(subset=["trip", "site", "compartment"])
                 .groupby(["trip", "site", "compartment"])["chemodiversity_H"]
                 .mean().reset_index())

    sh_cell = (meta_sh.dropna(subset=["trip", "site", "compartment", "Shannon"])
               .groupby(["trip", "site", "compartment"])["Shannon"].mean()
               .reset_index())
    sh_cell["compartment"] = sh_cell["compartment"].str.lower()
    merged = chem_cell.merge(sh_cell, on=["trip", "site", "compartment"], how="inner")

    log(f"\n[S15 Chemodiversity vs Shannon, all trips, n={len(merged)} cells]")
    by_comp = (merged.groupby("compartment")["chemodiversity_H"]
               .agg(["mean", "std", "count"]))
    log("  Per-compartment elemental Shannon entropy:")
    log(by_comp.to_string())

    fig, axes = plt.subplots(1, 3, figsize=(11, 3.5), sharey=True)
    palette = {"surface": "#d1410c", "deep": "#1d4ed8", "rhizosphere": "#16a34a"}
    rows = []
    for ax, comp in zip(axes, ["surface", "deep", "rhizosphere"]):
        sub = merged[merged["compartment"] == comp]
        ax.scatter(sub["chemodiversity_H"], sub["Shannon"], s=22,
                   c=palette[comp], edgecolor="black", linewidth=0.3, alpha=0.75)
        if len(sub) >= 8:
            rho, p = stats.spearmanr(sub["chemodiversity_H"], sub["Shannon"])
        else:
            rho, p = float("nan"), float("nan")
        ax.set_title(f"{comp}: ρ = {rho:+.2f}, p = {p:.3g}, n = {len(sub)}")
        ax.set_xlabel("Chemodiversity (H')")
        ax.spines[["top", "right"]].set_visible(False)
        rows.append({"compartment": comp, "rho": rho, "p": p, "n": len(sub)})
    axes[0].set_ylabel("Microbial Shannon")
    fig.tight_layout()
    out = FIG / "supp_RQ16_Chemodiversity_fig1_chemodiv_vs_shannon.pdf"
    fig.savefig(out, bbox_inches="tight")
    plt.close(fig)
    for r in rows:
        log(f"  {r['compartment']:11s} rho = {r['rho']:+.3f}, p = {r['p']:.3g}, n = {r['n']}")
    return by_comp, pd.DataFrame(rows)


def varnish_mn_radiation(xrf, ft, tax, meta_sh):
    """S16: Mn vs radiation-resistant genera abundance."""
    # genus relab per 16S sample
    if "genus" not in tax.columns and "Genus" in tax.columns:
        tax = tax.rename(columns={"Genus": "genus"})
    asv_to_genus = tax["genus"].reindex(ft.index).fillna("Unclassified")
    rel = ft.div(ft.sum(axis=0), axis=1)
    rad_rows = rel.index.intersection(asv_to_genus.index[asv_to_genus.isin(RAD_RESISTANT)])
    rad_relab = rel.loc[rad_rows].sum(axis=0).rename("rad_relab")
    rad_per_genus = (rel.assign(genus=asv_to_genus.values).groupby("genus").sum())

    # match XRF cells to 16S samples
    xrf_cell = (xrf.dropna(subset=["trip", "site", "compartment"])
                .groupby(["trip", "site", "compartment"])["Mn"].mean().reset_index())
    fe_cell = (xrf.dropna(subset=["trip", "site", "compartment"])
               .groupby(["trip", "site", "compartment"])["Fe"].mean().reset_index())
    panel = (meta_sh.assign(compartment=meta_sh["compartment"].str.lower())
             .merge(xrf_cell, on=["trip", "site", "compartment"], how="inner")
             .merge(fe_cell, on=["trip", "site", "compartment"], how="inner"))
    panel["rad_relab"] = rad_relab.reindex(panel["sample"]).values

    n_total = len(panel)
    n_mn_below = int((panel["Mn"] == 0).sum())
    log(f"\n[S16 Desert varnish, all trips]")
    log(f"  Mn = 0 (below detection): {n_mn_below}/{n_total} "
        f"({100*n_mn_below/n_total:.1f}%) (was 79.7% in Trip-5-only)")
    for g in RAD_RESISTANT:
        if g in rad_per_genus.index:
            n_pos = int((rad_per_genus.loc[g, panel["sample"]] > 0).sum())
        else:
            n_pos = 0
        log(f"  {g:18s} detected in {n_pos}/{n_total} samples")

    fig, axes = plt.subplots(1, 3, figsize=(11, 3.6), sharey=True)
    rows = []
    palette = {"surface": "#d1410c", "deep": "#1d4ed8", "rhizosphere": "#16a34a"}
    for ax, comp in zip(axes, ["surface", "deep", "rhizosphere"]):
        sub = panel[panel["compartment"] == comp].dropna(subset=["Mn", "rad_relab"])
        ax.scatter(sub["Mn"], sub["rad_relab"], s=18, c=palette[comp],
                   edgecolor="black", linewidth=0.2, alpha=0.6)
        if len(sub) >= 8 and sub["Mn"].nunique() > 1:
            rho, p = stats.spearmanr(sub["Mn"], sub["rad_relab"])
        else:
            rho, p = float("nan"), float("nan")
        ax.set_title(f"{comp}: ρ = {rho:+.2f}, p = {p:.3g}, n = {len(sub)}")
        ax.set_xlabel("Mn (% dry mass)")
        ax.spines[["top", "right"]].set_visible(False)
        rows.append({"compartment": comp, "rho": rho, "p": p, "n": len(sub)})
    axes[0].set_ylabel("Radiation-resistant genera (rel abund.)")
    fig.tight_layout()
    out = FIG / "supp_RQ11_DesertVarnish_fig1_mn_vs_radiation_resistant.pdf"
    fig.savefig(out, bbox_inches="tight")
    plt.close(fig)
    for r in rows:
        log(f"  Mn x rad-relab ({r['compartment']:11s}): rho = {r['rho']:+.3f}, "
            f"p = {r['p']:.3g}, n = {r['n']}")
    return panel, rows


def sulfur_cycle_taxa(xrf, ft, tax, meta_sh):
    """S17: XRF sulfur vs SRB/SOB genera abundance."""
    if "genus" not in tax.columns and "Genus" in tax.columns:
        tax = tax.rename(columns={"Genus": "genus"})
    asv_to_genus = tax["genus"].reindex(ft.index).fillna("Unclassified")
    rel = ft.div(ft.sum(axis=0), axis=1)
    rel_genus = rel.assign(genus=asv_to_genus.values).groupby("genus").sum()

    srb_rows = [g for g in SRB_GENERA if g in rel_genus.index]
    sob_rows = [g for g in SOB_GENERA if g in rel_genus.index]
    srb_relab = rel_genus.loc[srb_rows].sum(axis=0) if srb_rows else pd.Series(0.0, index=rel.columns)
    sob_relab = rel_genus.loc[sob_rows].sum(axis=0) if sob_rows else pd.Series(0.0, index=rel.columns)

    s_cell = (xrf.dropna(subset=["trip", "site", "compartment"])
              .groupby(["trip", "site", "compartment"])["S"].mean().reset_index())
    panel = (meta_sh.assign(compartment=meta_sh["compartment"].str.lower())
             .merge(s_cell, on=["trip", "site", "compartment"], how="inner"))
    panel["srb_relab"] = srb_relab.reindex(panel["sample"]).values
    panel["sob_relab"] = sob_relab.reindex(panel["sample"]).values

    log(f"\n[S17 Sulfur cycle, all trips, n={len(panel)} samples with XRF]")
    log(f"  SRB genera tracked: {srb_rows}")
    log(f"  SOB genera tracked: {sob_rows}")
    for g in SRB_GENERA + SOB_GENERA:
        if g in rel_genus.index:
            n_pos = int((rel_genus.loc[g, panel["sample"]] > 0).sum())
            log(f"  {g:18s} detected in {n_pos}/{len(panel)}")

    fig, axes = plt.subplots(1, 2, figsize=(9, 3.6), sharey=False)
    for ax, (col, lab) in zip(axes,
                               [("srb_relab", "Sulfate-reducing bacteria"),
                                ("sob_relab", "Sulfur-oxidising bacteria")]):
        ax.scatter(panel["S"], panel[col], s=14, c="#444444",
                   edgecolor="black", linewidth=0.15, alpha=0.5)
        sub = panel.dropna(subset=["S", col])
        if len(sub) > 8 and sub[col].nunique() > 2:
            rho, p = stats.spearmanr(sub["S"], sub[col])
        else:
            rho, p = float("nan"), float("nan")
        ax.set_title(f"{lab}: ρ = {rho:+.2f}, p = {p:.3g}, n = {len(sub)}")
        ax.set_xlabel("XRF S (% dry mass)")
        ax.set_ylabel("Genus rel. abund.")
        ax.spines[["top", "right"]].set_visible(False)
        log(f"  S x {col} rho = {rho:+.3f}, p = {p:.3g}, n = {len(sub)}")
    fig.tight_layout()
    out = FIG / "supp_RQ12_SulfurCycle_fig1_s_vs_srb_sob.pdf"
    fig.savefig(out, bbox_inches="tight")
    plt.close(fig)
    return panel


def main() -> None:
    xrf, ft, tax, meta, meta_sh = load_inputs()
    panel_sc = site_compartment_panel(xrf, meta_sh)
    panel_with_pca, _ = lithology_pca(panel_sc)
    si_ca_vs_shannon(panel_sc)
    element_diversity_heatmap(panel_sc)
    chemodiversity_vs_shannon(xrf, meta_sh)
    varnish_mn_radiation(xrf, ft, tax, meta_sh)
    sulfur_cycle_taxa(xrf, ft, tax, meta_sh)

    out = CACHE / "supplement_xrf_stats.txt"
    with open(out, "w") as fh:
        fh.write("\n".join(stats_buf))
    print(f"\nwrote {out}")


if __name__ == "__main__":
    main()
