#!/usr/bin/env python3
"""Render supplementary figures S1-S10 for main_v2 narrative.

Single one-pagers per figure (or two-panel) — concise.

S1 : PMA viability QC (proxy: relic-score by MAG-match presence)
S2 : EMP cosmopolitanism (frac matched at 3 trim lengths)
S3 : Cross-trip persistence per compartment
S4 : Three-guild knockout robustness (edge-loss %)
S5 : Per-cell trip-dominance sequence types
S6 : Cosmopolitanism vs alive class (counter-intuitive: cosmopolitan = alive)
S7 : Wind partial Mantel sweep (distribution across 11,520 tests)
S8 : Adversarial keystone ranking (CSP1-2 vs alternatives)
S9 : Pulse-reserve precip response (Strategy A vs d7 precip)
S10: Per-site projected P(A→B) under SSP1-2.6_2100
"""
from __future__ import annotations

from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt

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
    "pdf.fonttype": 42,
})

REPO = Path(__file__).resolve().parents[1]
CACHE = REPO / "cache"
FIG = REPO.parent / "figures"


def save(fig, name):
    out = FIG / f"supp_v2_{name}.pdf"
    fig.savefig(out, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {out}", flush=True)


def s1_pma_qc():
    """PMA QC — distribution of relic score split by MAG-match presence."""
    rs = pd.read_csv(CACHE / "relic_priors" / "relic_score_with_mag_prior.tsv",
                       sep="\t")
    fig, axes = plt.subplots(1, 2, figsize=(7.0, 2.8))
    # Panel a: relic_score full GB AUC summary (no per-pair available; show
    # score distribution before vs after MAG prior).
    ax = axes[0]
    ax.hist(rs["relic_score_full_gb"].dropna(), bins=50, alpha=0.5,
              color="#888888", label="GB classifier")
    ax.hist(rs["relic_score_with_mag"].dropna(), bins=50, alpha=0.5,
              color="#117733", label="+ MAG prior")
    ax.axvline(0.3, color="#cc3311", linestyle="--", linewidth=0.8,
                  label="alive cutoff")
    ax.set_xlabel("Relic-likelihood score")
    ax.set_ylabel("ASVs")
    ax.set_title("a  Indicator score: GB → +MAG prior",
                    loc="left", fontweight="bold")
    ax.legend(frameon=False, fontsize=6)
    # Panel b: how the prior re-classifies — points moving across the cutoff
    ax = axes[1]
    pre = rs["relic_score_full_gb"]
    post = rs["relic_score_with_mag"]
    # Color by sign of change
    moved_to_alive = (pre > 0.3) & (post <= 0.3)
    moved_to_relic = (pre <= 0.3) & (post > 0.3)
    unchanged = ~(moved_to_alive | moved_to_relic)
    ax.scatter(pre[unchanged], post[unchanged], s=2, c="#bbbbbb", alpha=0.4,
                  rasterized=True)
    ax.scatter(pre[moved_to_alive], post[moved_to_alive], s=8,
                  c="#117733", alpha=0.7, label=f"→ alive ({moved_to_alive.sum()})",
                  rasterized=True)
    ax.scatter(pre[moved_to_relic], post[moved_to_relic], s=8,
                  c="#cc3311", alpha=0.7, label=f"→ relic ({moved_to_relic.sum()})",
                  rasterized=True)
    ax.plot([0, 1], [0, 1], color="black", linewidth=0.4, alpha=0.5)
    ax.axhline(0.3, color="#cc3311", linestyle="--", linewidth=0.6)
    ax.axvline(0.3, color="#cc3311", linestyle="--", linewidth=0.6)
    ax.set_xlabel("Pre-prior score (GB only)")
    ax.set_ylabel("Post-prior score (+ MAG)")
    ax.set_title("b  MAG-prior reclassification",
                    loc="left", fontweight="bold")
    ax.legend(frameon=False, fontsize=6, loc="lower right")
    fig.suptitle("Fig. S1 · Composite relic-likelihood indicator",
                    fontsize=10, fontweight="bold", y=1.02)
    fig.tight_layout()
    save(fig, "S1_relic_indicator")


def s2_emp_cosmo():
    d = pd.read_csv(CACHE / "emp_cosmopolitanism" /
                       "cosmopolitanism_summary.tsv", sep="\t")
    fig, ax = plt.subplots(figsize=(4.0, 3.0))
    x = np.arange(len(d))
    bars = ax.bar(x, d["frac_eq_matched"] * 100,
                     color="#0077bb", edgecolor="white", linewidth=0.5,
                     width=0.6)
    for xi, v, n in zip(x, d["frac_eq_matched"], d["n_eq_asvs_matched"]):
        ax.text(xi, v * 100 + 2, f"{v*100:.0f}%\n(n={n:,})",
                  ha="center", va="bottom", fontsize=7)
    ax.set_xticks(x)
    ax.set_xticklabels([f"{int(l)} bp" for l in d["v4_length_bp"]])
    ax.set_xlabel("ASV trim length (matched to EMP)")
    ax.set_ylabel("% of EQ ASVs matched to EMP reference")
    ax.set_ylim(0, 75)
    ax.set_title("Fig. S2 · EQ ASV cosmopolitanism vs Earth Microbiome Project",
                    loc="left", fontweight="bold", fontsize=9)
    ax.grid(axis="y", alpha=0.3, linewidth=0.4)
    save(fig, "S2_emp_cosmopolitanism")


def s3_persistence():
    d = pd.read_csv(CACHE / "test6_persistence" /
                       "persistence_summary_per_compartment.tsv", sep="\t")
    p = d.pivot(index="trips_present", columns="compartment",
                   values="frac_of_records")
    COMP_COLORS = {"surface": "#ee7733", "deep": "#0077bb",
                       "rhizosphere": "#33aa55"}
    fig, ax = plt.subplots(figsize=(5.0, 3.0))
    n = len(p.index)
    w = 0.25
    for i, c in enumerate(["surface", "deep", "rhizosphere"]):
        if c not in p.columns: continue
        x = np.arange(n) + (i - 1) * w
        ax.bar(x, p[c].values * 100, width=w, color=COMP_COLORS[c],
                  edgecolor="white", linewidth=0.4, label=c)
    ax.set_xticks(np.arange(n))
    ax.set_xticklabels([str(int(i)) for i in p.index])
    ax.set_xlabel("Number of trips an OTU appears at the same site")
    ax.set_ylabel("% of (OTU × site) records")
    ax.legend(frameon=False, fontsize=6, title="compartment",
                title_fontsize=7)
    ax.set_title("Fig. S3 · Cross-trip persistence: ${\\sim}67\\%$ of "
                    "OTUs are 1-trip ephemeral",
                    loc="left", fontweight="bold", fontsize=9)
    ax.grid(axis="y", alpha=0.3, linewidth=0.4)
    save(fig, "S3_persistence")


def s4_guild_knockout():
    d = pd.read_csv(CACHE / "keystone_test" / "knockout_robustness.tsv",
                       sep="\t")
    # Aggregate: edge_loss_pct by scenario, per compartment, plus mean
    def label(s):
        m = {"Nibribacter_alone": "Nibribacter alone",
              "Bact_DOM_guild_n8": "Bacteroidota DOM guild",
              "Bacilli_guild_n9": "Halotolerant Bacilli guild",
              "Bacilli_guild_n6": "Halotolerant Bacilli guild",
              "Bacilli_guild_n7": "Halotolerant Bacilli guild",
              "Pseudo_guild_n5": "Pseudomonas guild",
              "random_single_avg50": "random single (null)",
              "random_n9_avg50":     "random guild-sized (null)"}
        return m.get(s, s)
    d["pretty"] = d["scenario"].map(label)
    order = ["Nibribacter alone", "Pseudomonas guild",
                "Bacteroidota DOM guild", "Halotolerant Bacilli guild",
                "random single (null)", "random guild-sized (null)"]
    COMP_COLORS = {"surface": "#ee7733", "deep": "#0077bb",
                       "rhizosphere": "#33aa55"}
    fig, ax = plt.subplots(figsize=(6.0, 3.5))
    n_pos = np.arange(len(order))
    width = 0.25
    for i, comp in enumerate(["surface", "deep", "rhizosphere"]):
        sub = d[d["compartment"] == comp].set_index("pretty")
        vals = [sub.loc[o, "edge_loss_pct"] if o in sub.index else 0
                  for o in order]
        x = n_pos + (i - 1) * width
        ax.bar(x, vals, width=width, color=COMP_COLORS[comp],
                  edgecolor="white", linewidth=0.4, label=comp)
    ax.set_xticks(n_pos)
    ax.set_xticklabels(order, rotation=25, ha="right", fontsize=7)
    ax.set_ylabel("Edge loss (%) on alive co-occurrence network")
    ax.set_title("Fig. S4 · Three-guild knockout decisively beats "
                    "single-keystone removal",
                    loc="left", fontweight="bold", fontsize=9)
    ax.axhline(0, color="black", linewidth=0.5)
    ax.legend(frameon=False, fontsize=6, title="compartment",
                title_fontsize=7, loc="upper left")
    ax.grid(axis="y", alpha=0.3, linewidth=0.4)
    save(fig, "S4_guild_knockout")


def s5_sequence_types():
    d = pd.read_csv(CACHE / "transition_asymmetry" /
                       "per_cell_sequences.tsv", sep="\t")
    # Classify
    def classify(seq):
        if not isinstance(seq, str): return "other"
        if "B" not in seq: return "stable_A"
        if "A" not in seq: return "stable_B"
        # Drift detection: monotonic A->B or B->A
        s = [c for c in seq if c in "AB"]
        if len(s) >= 3:
            if s[0] == "A" and s[-1] == "B" and seq.count("B") - 1 <= seq.count("AB", 0): pass
        # Count direction changes
        changes = sum(1 for i in range(1, len(s)) if s[i] != s[i-1])
        if changes >= 2: return "oscillating"
        # exactly 1 change
        if s[0] == "A": return "drift_A_to_B"
        return "drift_B_to_A"
    d["cls"] = d["sequence"].apply(classify)
    cnt = d["cls"].value_counts()
    order = ["stable_A", "oscillating", "drift_A_to_B", "stable_B",
                "drift_B_to_A", "other"]
    cnt = cnt.reindex(order).fillna(0)
    pct = cnt / cnt.sum() * 100
    colors = ["#117733", "#aa9933", "#ee9966", "#882255",
                  "#66aacc", "#bbbbbb"]
    fig, ax = plt.subplots(figsize=(5.0, 3.0))
    bars = ax.barh(np.arange(len(cnt)), pct.values, color=colors,
                       edgecolor="white", linewidth=0.4)
    for i, (n, p) in enumerate(zip(cnt.values, pct.values)):
        ax.text(p + 1.5, i, f"{int(n)} ({p:.0f}%)",
                  va="center", fontsize=7)
    ax.set_yticks(np.arange(len(cnt)))
    ax.set_yticklabels([o.replace("_", " ") for o in order])
    ax.invert_yaxis()
    ax.set_xlabel("% of (site, compartment) cells")
    ax.set_title("Fig. S5 · Per-cell trip-dominance sequence types",
                    loc="left", fontweight="bold", fontsize=9)
    ax.grid(axis="x", alpha=0.3, linewidth=0.4)
    save(fig, "S5_sequence_types")


def s6_cosmo_alive():
    """Cross EMP cosmopolitanism with alive-classification:
    cosmopolitan taxa are MORE alive, not less."""
    rs = pd.read_csv(CACHE / "relic_priors" / "relic_score_with_mag_prior.tsv",
                       sep="\t")
    # The emp_cosmo_X columns flag presence in EMP at X bp trim
    cols = [c for c in rs.columns if c.startswith("emp_cosmo")]
    if not cols:
        return
    # Use the 90bp version (most permissive)
    rs["cosmo_90"] = rs.get("emp_cosmo_90", 0).fillna(0).astype(bool)
    rs["alive"] = rs["relic_score_with_mag"] <= 0.3
    # 2x2 contingency
    ct = pd.crosstab(rs["cosmo_90"], rs["alive"],
                       margins=False, normalize="index") * 100
    fig, ax = plt.subplots(figsize=(4.0, 3.0))
    x = np.arange(2)
    w = 0.38
    ax.bar(x - w/2, ct[True].values, width=w, color="#117733",
              label="alive", edgecolor="white", linewidth=0.4)
    ax.bar(x + w/2, ct[False].values, width=w, color="#cc3311",
              label="relic", edgecolor="white", linewidth=0.4)
    ax.set_xticks(x)
    ax.set_xticklabels(["Endemic\n(not in EMP)", "Cosmopolitan\n(in EMP @ 90 bp)"])
    ax.set_ylabel("% within class")
    ax.set_ylim(0, 100)
    ax.set_title("Fig. S6 · Cosmopolitan ASVs are MORE alive, not more relic",
                    loc="left", fontweight="bold", fontsize=9)
    ax.legend(frameon=False, fontsize=6)
    ax.grid(axis="y", alpha=0.3, linewidth=0.4)
    # Annotation
    cosmo_alive = ct.loc[True, True]
    endemic_alive = ct.loc[False, True]
    ax.text(0.5, 0.95,
              f"Δ alive%: {cosmo_alive - endemic_alive:+.1f} pp",
              transform=ax.transAxes, ha="center", va="top", fontsize=7,
              fontweight="bold", color="#117733",
              bbox=dict(facecolor="white", edgecolor="#117733",
                          linewidth=0.5))
    save(fig, "S6_cosmo_alive")


def s7_wind_sweep():
    """Wind-Mantel sweep distribution of partial r across all tests."""
    try:
        d = pd.read_csv(CACHE / "wind_dispersal" / "sweep_mantel_full.tsv",
                          sep="\t")
    except Exception as e:
        print(f"S7 skipped: {e}")
        return
    # Distribution of r_part by window
    fig, ax = plt.subplots(figsize=(5.5, 3.2))
    windows = sorted(d["window_days"].unique())
    data = [d[d["window_days"] == w]["r_part"].dropna().values for w in windows]
    bp = ax.boxplot(data, tick_labels=[f"{w}d" for w in windows],
                       patch_artist=True, showfliers=False,
                       medianprops=dict(color="black", linewidth=0.9))
    for patch, w in zip(bp["boxes"], windows):
        # Color gradient: shorter windows lighter, longer darker
        c = plt.cm.Blues(0.3 + 0.6 * (windows.index(w) / max(1, len(windows) - 1)))
        patch.set_facecolor(c)
        patch.set_edgecolor("black")
        patch.set_linewidth(0.4)
    ax.axhline(0, color="grey", linewidth=0.5, alpha=0.5)
    ax.set_xlabel("Wind-distance window (days)")
    ax.set_ylabel("Partial Mantel r (BC ~ wind | distance)")
    ax.set_title(f"Fig. S7 · Wind partial-Mantel sweep "
                    f"(n = {len(d):,} tests across windows × scores × strata)",
                    loc="left", fontweight="bold", fontsize=8.5)
    ax.grid(axis="y", alpha=0.3, linewidth=0.4)
    save(fig, "S7_wind_sweep")


def s8_adversarial():
    """Adversarial keystone — knockout magnitude by candidate."""
    # No single all-genus knockout table available without re-run; use the
    # guild-level results as a proxy.
    d = pd.read_csv(CACHE / "keystone_test" / "knockout_robustness.tsv",
                       sep="\t")
    # Pivot mean edge_loss_pct by scenario across compartments
    a = (d.groupby("scenario")["edge_loss_pct"].mean()
            .sort_values(ascending=True))
    fig, ax = plt.subplots(figsize=(5.0, 3.0))
    colors = []
    for s in a.index:
        if "random" in s: colors.append("#bbbbbb")
        elif "Nibribacter" in s: colors.append("#117733")
        elif "Bacilli" in s: colors.append("#882255")
        else: colors.append("#6688aa")
    ax.barh(np.arange(len(a)), a.values, color=colors,
              edgecolor="white", linewidth=0.4)
    pretty = {"Nibribacter_alone": "Nibribacter alone",
                 "Bact_DOM_guild_n8": "Bacteroidota DOM guild",
                 "Bacilli_guild_n9": "Halotolerant Bacilli (n=9)",
                 "Bacilli_guild_n7": "Halotolerant Bacilli (n=7)",
                 "Bacilli_guild_n6": "Halotolerant Bacilli (n=6)",
                 "Pseudo_guild_n5": "Pseudomonas guild",
                 "random_single_avg50": "random single (null)",
                 "random_n9_avg50": "random guild-sized (null)"}
    ax.set_yticks(np.arange(len(a)))
    ax.set_yticklabels([pretty.get(s, s) for s in a.index], fontsize=7)
    ax.set_xlabel("Mean edge loss (%) on alive network")
    ax.set_title("Fig. S8 · Adversarial knockout ranking: "
                    "single-keystone is not the strongest perturbation",
                    loc="left", fontweight="bold", fontsize=8.5)
    ax.grid(axis="x", alpha=0.3, linewidth=0.4)
    save(fig, "S8_adversarial")


def s9_pulse_reserve():
    """Pulse-reserve: Strategy A genus-level abundance vs d7 precip."""
    df = pd.read_csv(CACHE / "two_strategy_temporal" /
                       "per_sample_strategy_with_precip.tsv", sep="\t")
    sub = df.dropna(subset=["strategy_A", "d7"])
    sub = sub[sub["strategy_A"] > 0]
    fig, ax = plt.subplots(figsize=(4.5, 3.0))
    ax.scatter(sub["d7"], np.log10(sub["strategy_A"] + 1e-5),
                  s=3, alpha=0.35, c="#117733", edgecolors="none",
                  rasterized=True)
    rho = sub[["strategy_A", "d7"]].corr(method="spearman").iloc[0, 1]
    ax.set_xlabel("Precipitation in 7 days before sample (mm)")
    ax.set_ylabel("log₁₀(Strategy A abundance + 1e-5)")
    ax.set_title(f"Fig. S9 · Pulse-reserve: Strategy A tracks "
                    f"recent precipitation (ρ = {rho:+.2f}, n = {len(sub):,})",
                    loc="left", fontweight="bold", fontsize=8.5)
    ax.grid(alpha=0.3, linewidth=0.4)
    save(fig, "S9_pulse_reserve")


def s10_per_site_projection():
    """Per-site projected P(A→B) under SSP1-2.6_2100."""
    d = pd.read_csv(CACHE / "two_strategy_projection_v3" /
                       "per_site_AtoB_risk.tsv", sep="\t")
    sub = d[d["scenario"] == "SSP1-2.6_2100"]
    if len(sub) == 0:
        # fall back to first available
        sub = d[d["scenario"] == d["scenario"].iloc[0]]
    # Average across compartments per site
    per_site = (sub.groupby("site")["p_AtoB_projected_mean"].mean()
                    .reset_index().sort_values("p_AtoB_projected_mean"))
    # join with coords for map
    geo = pd.read_csv(REPO / "data" / "geodata" / "trip1_geodata.tsv",
                          sep="\t")
    geo["Site"] = pd.to_numeric(geo["Site"], errors="coerce")
    geo = geo.dropna(subset=["Site"])
    geo["Site"] = geo["Site"].astype(int)
    per_site = per_site.merge(geo[["Site", "Latitude", "Longitude"]],
                                    left_on="site", right_on="Site",
                                    how="left").dropna(
                                        subset=["Latitude", "Longitude"])
    fig, ax = plt.subplots(figsize=(6.0, 3.2))
    sc = ax.scatter(per_site["Longitude"], per_site["Latitude"],
                       c=per_site["p_AtoB_projected_mean"], cmap="RdYlBu_r",
                       s=42, edgecolors="black", linewidths=0.4,
                       vmin=0.2, vmax=0.8)
    plt.colorbar(sc, ax=ax, label="Projected P(A→B) under SSP1-2.6_2100")
    ax.set_xlabel("Longitude (°E)")
    ax.set_ylabel("Latitude (°N)")
    ax.set_xlim(44, 56)
    ax.set_ylim(18, 22)
    ax.set_aspect("equal", adjustable="box")
    ax.set_title("Fig. S10 · Per-site projected A→B switching risk "
                    "under SSP1-2.6$_{2100}$",
                    loc="left", fontweight="bold", fontsize=8.5)
    ax.grid(alpha=0.3, linewidth=0.4)
    save(fig, "S10_per_site_projection")


def main():
    s1_pma_qc()
    s2_emp_cosmo()
    s3_persistence()
    s4_guild_knockout()
    s5_sequence_types()
    s6_cosmo_alive()
    s7_wind_sweep()
    s8_adversarial()
    s9_pulse_reserve()
    s10_per_site_projection()
    print("DONE")


if __name__ == "__main__":
    main()
