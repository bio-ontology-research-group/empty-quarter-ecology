#!/usr/bin/env python3
"""Main Fig 4 (Transitions + climate forecast).

Panels:
  (a) Transition matrix (A→A 88%, A→B 12%, B→A 48%, B→B 52%).
  (b) Logit coefficients on P(A→B) and P(B→A).
  (c) CMIP6 π_B trajectory with bootstrap 95% CI bars and LTO range.
  (d) T-only vs P-only vs combined decomposition (SSP3-7.0_2100).

Reads:
  cache/transition_asymmetry/all_transitions.tsv
  cache/two_strategy_projection_v3/{scenario_summary_v3.tsv,
        decomposition_ssp370_2100.tsv}
  cache/two_strategy_projection_v3/uncertainty/{bootstrap_ci.tsv,
        lto_trip_pair.tsv}
"""
from __future__ import annotations

from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import StandardScaler

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
DATA = REPO / "data"
FIG = REPO.parent / "figures"

C_A = "#117733"
C_B = "#882255"


def panel_a_transition_matrix(ax):
    tr = pd.read_csv(CACHE / "transition_asymmetry" / "all_transitions.tsv",
                       sep="\t")
    counts = tr.groupby(["from_dom", "to_dom"]).size().unstack(fill_value=0)
    # Order
    counts = counts.reindex(index=["A", "B"], columns=["A", "B"], fill_value=0)
    # Row-normalise to get transition probabilities
    probs = counts.div(counts.sum(axis=1), axis=0)
    im = ax.imshow(probs.values, cmap="Blues", vmin=0, vmax=1,
                       aspect="equal")
    for i in range(2):
        for j in range(2):
            v = probs.iloc[i, j]
            n = counts.iloc[i, j]
            ax.text(j, i, f"{v*100:.0f}%\n(n={n})",
                       ha="center", va="center", fontsize=9,
                       fontweight="bold",
                       color="white" if v > 0.5 else "black")
    ax.set_xticks([0, 1]); ax.set_xticklabels(["A", "B"])
    ax.set_yticks([0, 1]); ax.set_yticklabels(["A", "B"])
    ax.set_xlabel("To trip")
    ax.set_ylabel("From trip")
    ax.set_title("a  Transition asymmetry\n(n = {:,} pairs)".format(
        len(tr)), loc="left", fontweight="bold", fontsize=8.5)
    # Annotation
    n_AA = counts.loc["A", "A"]; n_A_total = counts.loc["A"].sum()
    n_BB = counts.loc["B", "B"]; n_B_total = counts.loc["B"].sum()
    ax.text(1.65, 0.0, f"A persistence: {n_AA/n_A_total*100:.0f}%\n"
                            f"(resilient default)",
              ha="left", va="center", fontsize=7, color=C_A,
              transform=ax.transData)
    ax.text(1.65, 1.0, f"B persistence: {n_BB/n_B_total*100:.0f}%\n"
                            f"(transient)",
              ha="left", va="center", fontsize=7, color=C_B,
              transform=ax.transData)
    ax.set_xlim(-0.5, 4.5)


def _build_climate_features():
    tr = pd.read_csv(CACHE / "transition_asymmetry" / "all_transitions.tsv",
                       sep="\t")
    np_daily = pd.read_parquet(CACHE / "nasa_power_daily.parquet")
    np_daily["Date"] = pd.to_datetime(np_daily["Date"])
    geo = []
    for t in (1, 2, 3, 4, 5):
        g = pd.read_csv(DATA / "geodata" / f"trip{t}_geodata.tsv", sep="\t")
        g["trip"] = t
        g["Site"] = pd.to_numeric(g["Site"], errors="coerce")
        g = g.dropna(subset=["Site"])
        g["Site"] = g["Site"].astype(int)
        g["CenterDate"] = pd.to_datetime(g["CenterDate"])
        geo.append(g[["Site", "trip", "CenterDate"]])
    td = pd.concat(geo, ignore_index=True).drop_duplicates(["Site", "trip"])
    rows = []
    for _, r in td.iterrows():
        s = r["Site"]; tt = r["trip"]; cd = r["CenterDate"]
        sub = np_daily[np_daily["Site"] == s]
        rec = {"site": s, "trip": tt}
        for w_d, w_lab in [(30, "T_d30"), (90, "T_d90"), (365, "T_d365")]:
            w = sub[(sub["Date"] >= cd - pd.Timedelta(days=w_d)) &
                       (sub["Date"] < cd)]
            rec[w_lab] = w["TS"].mean() if len(w) else np.nan
        rows.append(rec)
    site_T = pd.DataFrame(rows)
    site_T_to = site_T.rename(columns={"trip": "to_trip",
                                          "T_d30": "to_T_d30",
                                          "T_d90": "to_T_d90",
                                          "T_d365": "to_T_d365"})
    site_T_from = site_T.rename(columns={"trip": "from_trip",
                                            "T_d30": "from_T_d30",
                                            "T_d90": "from_T_d90",
                                            "T_d365": "from_T_d365"})
    tr = tr.merge(site_T_to, on=["site", "to_trip"], how="left")
    tr = tr.merge(site_T_from, on=["site", "from_trip"], how="left")
    return tr


def panel_b_coefficients(ax):
    """Fit the two logits and plot coefficients."""
    tr = _build_climate_features()
    FEAT = ["to_d7", "to_d365", "delta_d30", "delta_d90", "delta_d180",
              "to_T_d30", "to_T_d90", "to_T_d365"]
    A = tr[tr["from_dom"] == "A"].dropna(subset=FEAT).copy()
    A["y"] = (A["transition"] == "A->B").astype(int)
    sA = StandardScaler().fit(A[FEAT].values)
    lrA = LogisticRegression(C=1.0, max_iter=2000, class_weight="balanced",
                                  random_state=42).fit(sA.transform(A[FEAT].values),
                                                         A["y"].values)
    B = tr[tr["from_dom"] == "B"].dropna(subset=FEAT).copy()
    B["y"] = (B["transition"] == "B->A").astype(int)
    sB = StandardScaler().fit(B[FEAT].values)
    lrB = LogisticRegression(C=1.0, max_iter=2000, class_weight="balanced",
                                  random_state=42).fit(sB.transform(B[FEAT].values),
                                                         B["y"].values)

    feat_pretty = {"to_d7": "precip d7",
                       "to_d365": "precip d365",
                       "delta_d30": "Δ precip d30",
                       "delta_d90": "Δ precip d90",
                       "delta_d180": "Δ precip d180",
                       "to_T_d30": "T d30",
                       "to_T_d90": "T d90",
                       "to_T_d365": "T d365"}
    y_pos = np.arange(len(FEAT))
    ax.barh(y_pos - 0.18, lrA.coef_[0], height=0.36, color=C_A,
              label="P(A→B)", edgecolor="white", linewidth=0.4)
    ax.barh(y_pos + 0.18, lrB.coef_[0], height=0.36, color=C_B,
              label="P(B→A)", edgecolor="white", linewidth=0.4)
    ax.axvline(0, color="grey", linewidth=0.6, alpha=0.6)
    ax.set_yticks(y_pos)
    ax.set_yticklabels([feat_pretty[f] for f in FEAT], fontsize=7)
    ax.invert_yaxis()
    ax.set_xlabel("Standardised logit coefficient")
    ax.set_title("b  Climate effects on transitions",
                    loc="left", fontweight="bold")
    ax.legend(frameon=False, fontsize=6, loc="lower right")
    ax.grid(axis="x", alpha=0.3, linewidth=0.4)
    # Annotate the T_d365 row
    i = FEAT.index("to_T_d365")
    ax.text(lrA.coef_[0][i] + 0.05, i - 0.18,
              f"{lrA.coef_[0][i]:+.2f}", fontsize=6, va="center",
              fontweight="bold", color=C_A)
    ax.text(lrB.coef_[0][i] - 0.05, i + 0.18,
              f"{lrB.coef_[0][i]:+.2f}", fontsize=6, va="center",
              ha="right", fontweight="bold", color=C_B)


def panel_c_projection(ax):
    """CMIP6 pi_B with bootstrap 95% CI bars + LTO range markers."""
    boot = pd.read_csv(CACHE / "two_strategy_projection_v3" / "uncertainty" /
                          "bootstrap_ci.tsv", sep="\t")
    lto = pd.read_csv(CACHE / "two_strategy_projection_v3" / "uncertainty" /
                         "lto_trip_pair.tsv", sep="\t")
    order = ["Historical", "SSP1-2.6_2100", "SSP2-4.5_2100", "SSP3-7.0_2100"]
    boot["o"] = boot["scenario"].map(lambda s: order.index(s))
    boot = boot.sort_values("o")
    x = np.arange(len(order))
    # Bootstrap CI as boxplot-like vertical
    point = boot["pi_B_point"].values
    lo = boot["pi_B_lo"].values
    hi = boot["pi_B_hi"].values
    yerr = np.vstack([point - lo, hi - point])
    colours = ["#444444", "#5588cc", "#cc6633", "#cc2222"]
    # Plot per-point with its own colour
    for xi, p, lo_i, hi_i, c in zip(x, point, lo, hi, colours):
        ax.errorbar(xi, p, yerr=[[p - lo_i], [hi_i - p]], fmt="o",
                       markersize=8, markerfacecolor=c, markeredgecolor="black",
                       markeredgewidth=0.5, ecolor=c, capsize=4,
                       linewidth=1.5, zorder=3)
    # LTO range as a small jittered scatter at SSP1-2.6_2100 position
    lto_ssp126 = lto["pi_B_ssp126_2100"].dropna().values
    j = np.random.RandomState(1).uniform(-0.18, 0.18, size=len(lto_ssp126))
    ax.scatter(np.full_like(lto_ssp126, 1) + j, lto_ssp126,
                  s=20, c="#5588cc", alpha=0.4, edgecolor="black",
                  linewidth=0.3, zorder=2,
                  label=f"LTO splits (n={len(lto_ssp126)})")
    # Annotate non-overlapping CIs
    ax.axhline(0.486, color="grey", linewidth=0.4, alpha=0.4)
    ax.axhline(0.660, color="grey", linewidth=0.4, alpha=0.4)
    ax.fill_between([-0.4, 0.4], 0.418, 0.486, color="#cccccc", alpha=0.25,
                       zorder=1)
    ax.set_xticks(x)
    ax.set_xticklabels([s.replace("_", "\n") for s in order], fontsize=7)
    ax.set_ylim(0.30, 1.05)
    ax.set_ylabel("Equilibrium B-fraction π_B")
    ax.set_title("c  CMIP6 projection · bootstrap 95% CI",
                    loc="left", fontweight="bold")
    ax.grid(axis="y", alpha=0.3, linewidth=0.4)
    # Headline annotation
    ax.text(1, 0.85,
              "SSP1-2.6\n0.46 → 0.76\n[0.66, 0.83]",
              fontsize=7, ha="center", va="bottom",
              color="#5588cc",
              bbox=dict(facecolor="white", edgecolor="#5588cc",
                          linewidth=0.6))
    ax.legend(frameon=False, fontsize=6, loc="lower right")


def panel_d_decomposition(ax):
    d = pd.read_csv(CACHE / "two_strategy_projection_v3" /
                       "decomposition_ssp370_2100.tsv", sep="\t")
    # Scenario column has Historical / T-only / P-only / Combined
    order = ["Historical", "T-only", "P-only", "Combined"]
    d["o"] = d["scenario"].map(lambda s: order.index(s) if s in order else 99)
    d = d.sort_values("o")
    x = np.arange(len(d))
    colors = ["#444444", "#cc2222", "#5588cc", "#882255"]
    bars = ax.bar(x, d["pi_B"].values, color=colors, edgecolor="white",
                      linewidth=0.5, width=0.65)
    for xi, v in zip(x, d["pi_B"].values):
        ax.text(xi, v + 0.02, f"{v:.2f}", ha="center", va="bottom",
                  fontsize=7, fontweight="bold")
    ax.set_xticks(x)
    ax.set_xticklabels(d["scenario"].values, fontsize=7)
    ax.set_ylim(0, 1.1)
    ax.set_ylabel("π_B")
    ax.set_title("d  Decomposition (SSP3-7.0_2100, +4 °C, −10% P)",
                    loc="left", fontweight="bold")
    ax.grid(axis="y", alpha=0.3, linewidth=0.4)
    ax.text(0.98, 0.96,
              "Temperature dominates ~10× over precipitation",
              transform=ax.transAxes, ha="right", va="top", fontsize=7,
              fontstyle="italic", color="#555555",
              bbox=dict(facecolor="white", edgecolor="none", alpha=0.85))


def main():
    fig = plt.figure(figsize=(7.2, 6.0))
    gs = GridSpec(2, 2, figure=fig, wspace=0.36, hspace=0.55,
                     left=0.10, right=0.96, top=0.94, bottom=0.10)
    ax_a = fig.add_subplot(gs[0, 0])
    ax_b = fig.add_subplot(gs[0, 1])
    ax_c = fig.add_subplot(gs[1, 0])
    ax_d = fig.add_subplot(gs[1, 1])
    panel_a_transition_matrix(ax_a)
    panel_b_coefficients(ax_b)
    panel_c_projection(ax_c)
    panel_d_decomposition(ax_d)
    out = FIG / "fig4_transitions_climate_v2.pdf"
    fig.savefig(out, dpi=300, bbox_inches="tight")
    print(f"wrote {out}", flush=True)


if __name__ == "__main__":
    main()
