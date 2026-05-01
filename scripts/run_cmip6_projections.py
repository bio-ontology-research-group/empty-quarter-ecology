"""
CMIP6-driven climate projections for the Empty Quarter digital twin.

Re-fits the Tier-3 NumPyro hierarchical state-space model
(notebook 11) and applies do() interventions for SSP1-2.6, SSP2-4.5,
and SSP3-7.0 at 2050 and 2100, using CMIP6 ensemble Δ(T, P) values
for the Arabian Peninsula from Almazroui et al. 2020
(Earth Syst Environ 4:611–630, doi:10.1007/s41748-020-00183-5;
31 CMIP6 GCMs, regional ensemble means relative to 1995–2014).

Two pathway formulations:
  (a) "direct" — applies CMIP6 ΔT and ΔP to the model's temp_mean_W30d
      and rain_W30d forcings only. Salinity is held fixed.
  (b) "with_S_amp" — additionally raises the S z-score by
      S_AMP_PER_DEGC × ΔT, encoding the warming-driven evaporative
      concentration of solutes in shallow groundwater / sabkha pans.

Reclamation reference: do(S −2 SD) from the same fit, for direct
comparison of climate forcing vs management intervention magnitudes.

Outputs:
  cache/cmip6_interventions.tsv
  figures/fig_cmip6_projections.pdf
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import jax
import jax.numpy as jnp
import numpyro
import numpyro.distributions as dist
from numpyro.infer import MCMC, NUTS

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO / "src"))
from eq import CACHE_DIR, FIGURES_DIR  # noqa: E402

# ---------------------------------------------------------------------
#  CMIP6 deltas for the Arabian Peninsula (Almazroui et al. 2020)
# ---------------------------------------------------------------------
#  (scenario, horizon): (ΔT °C, ΔP % annual)
#  Conservative ensemble means, mid-21st (~2040–2059) and end-21st
#  (~2080–2099) horizons relative to the 1995–2014 baseline.
CMIP6_DELTAS: dict[tuple[str, str], tuple[float, float]] = {
    ("SSP1-2.6", "2050"): (+1.2, -2.0),
    ("SSP1-2.6", "2100"): (+1.6, -3.0),
    ("SSP2-4.5", "2050"): (+1.7, -3.0),
    ("SSP2-4.5", "2100"): (+2.7, -5.0),
    ("SSP3-7.0", "2050"): (+2.1, -5.0),
    ("SSP3-7.0", "2100"): (+4.0, -10.0),
}

#  Salinity amplification (SD per °C warming).  Bounding estimate from
#  evaporative concentration of solutes in arid evaporite pans;
#  conservative — the central estimate, not a worst case.
S_AMP_PER_DEGC = 0.15

COMP_NAMES = {0: "surface", 1: "deep", 2: "rhizosphere"}

# ---------------------------------------------------------------------
#  1. Aggregate the causal panel (mirrors notebook 11)
# ---------------------------------------------------------------------
frame = pd.read_parquet(CACHE_DIR / "causal_frame_tier1.parquet")
group_cols = ["site", "trip", "compartment_num"]
agg_num = frame.groupby(group_cols).agg(
    shannon=("shannon", "mean"),
    rain_W30d=("rain_W30d", "mean"),
    temp_mean_W30d=("temp_mean_W30d", "mean"),
).reset_index()

site_geo = (
    frame.dropna(subset=["S", "P"])
    .groupby("site")[["S", "P", "Cl", "Na"]]
    .mean()
    .reset_index()
)
agg = agg_num.merge(site_geo, on="site", how="left")
agg = agg.dropna(
    subset=["shannon", "rain_W30d", "temp_mean_W30d", "S"]
).reset_index(drop=True)
print(f"panel: {len(agg)} rows; {agg.site.nunique()} sites; "
      f"{agg.compartment_num.nunique()} compartments")

# ---------------------------------------------------------------------
#  2. Fit Tier-3 hierarchical model
# ---------------------------------------------------------------------
def model(site_idx, comp_idx, forcing, y, n_sites, n_comp, n_force):
    mu_alpha = numpyro.sample("mu_alpha", dist.Normal(0, 2))
    tau_alpha = numpyro.sample("tau_alpha", dist.HalfNormal(1.0))
    with numpyro.plate("sites", n_sites):
        alpha = numpyro.sample("alpha", dist.Normal(mu_alpha, tau_alpha))
    with numpyro.plate("comps", n_comp):
        beta_comp = numpyro.sample("beta_comp", dist.Normal(0, 1))
    mu_gamma = numpyro.sample("mu_gamma", dist.Normal(jnp.zeros(n_force), 1))
    tau_gamma = numpyro.sample("tau_gamma", dist.HalfNormal(jnp.ones(n_force)))
    with numpyro.plate("g_comps", n_comp):
        with numpyro.plate("g_force", n_force):
            gamma = numpyro.sample(
                "gamma",
                dist.Normal(mu_gamma[:, None], tau_gamma[:, None]),
            )
    sigma = numpyro.sample("sigma", dist.HalfNormal(1.0))
    lin = (alpha[site_idx]
           + beta_comp[comp_idx]
           + jnp.sum(gamma[:, comp_idx].T * forcing, axis=1))
    numpyro.sample("obs", dist.Normal(lin, sigma), obs=y)


site_ids, site_inv = np.unique(agg.site, return_inverse=True)
comp_ids, comp_inv = np.unique(agg.compartment_num, return_inverse=True)
forcing_cols = ["rain_W30d", "temp_mean_W30d", "S"]
F = agg[forcing_cols].values
F_mean, F_sd = F.mean(axis=0), F.std(axis=0)
Fz = (F - F_mean) / F_sd
y_mean, y_sd = agg["shannon"].mean(), agg["shannon"].std()
yz = (agg["shannon"].values - y_mean) / y_sd

print("baseline forcings:")
for col, m, s in zip(forcing_cols, F_mean, F_sd):
    print(f"  {col:18s} mean={m:7.3f}  sd={s:7.3f}")

rng_key = jax.random.PRNGKey(0)
kernel = NUTS(model, target_accept_prob=0.9)
mcmc = MCMC(
    kernel, num_warmup=500, num_samples=1000, num_chains=2, progress_bar=False
)
mcmc.run(
    rng_key,
    site_idx=jnp.asarray(site_inv),
    comp_idx=jnp.asarray(comp_inv),
    forcing=jnp.asarray(Fz),
    y=jnp.asarray(yz),
    n_sites=len(site_ids),
    n_comp=len(comp_ids),
    n_force=len(forcing_cols),
)
samples = {k: np.asarray(v) for k, v in mcmc.get_samples().items()}

print("posterior mu_gamma (forcing coefficients on standardised Shannon):")
for k, name in enumerate(forcing_cols):
    q = np.percentile(samples["mu_gamma"][:, k], [2.5, 50, 97.5])
    print(f"  {name:18s} median={q[1]:+.3f}  95% CI=[{q[0]:+.3f}, {q[2]:+.3f}]")

# ---------------------------------------------------------------------
#  3. Apply do(climate scenario) interventions per compartment
# ---------------------------------------------------------------------
def delta_per_compartment(samples: dict, comp_id: int, delta_z: np.ndarray) -> np.ndarray:
    """
    Posterior draws of Δ(z-Shannon) for a uniform z-score shift
    applied to all forcings within one compartment.
    Because the model is linear in forcing, Δy = γ_c · Δz_F.
    Returns shape (n_draws,).
    """
    gamma = samples["gamma"]  # (n_draws, n_force, n_comp)
    return np.einsum("df,f->d", gamma[:, :, comp_id], delta_z)


records = []
for (sce, hor), (dT_C, dP_pct) in CMIP6_DELTAS.items():
    # Physical Δ → z-score Δ
    delta_phys = np.array(
        [
            (dP_pct / 100.0) * F_mean[0],  # ΔP applied to mean rainfall (mm/30d)
            dT_C,                          # ΔT additive (°C)
            0.0,                           # S unchanged in direct pathway
        ]
    )
    delta_z_direct = delta_phys / F_sd

    delta_z_amp = delta_z_direct.copy()
    delta_z_amp[2] = S_AMP_PER_DEGC * dT_C  # SD-units shift in S

    for comp_id, comp_name in COMP_NAMES.items():
        for label, dz in [("direct", delta_z_direct), ("with_S_amp", delta_z_amp)]:
            d = delta_per_compartment(samples, comp_id, dz)
            records.append(
                dict(
                    scenario=sce,
                    horizon=hor,
                    compartment=comp_name,
                    pathway=label,
                    delta_T_C=dT_C,
                    delta_P_pct=dP_pct,
                    delta_shannon_z_median=float(np.median(d)),
                    ci_lo=float(np.quantile(d, 0.025)),
                    ci_hi=float(np.quantile(d, 0.975)),
                    credible_95=bool(
                        np.quantile(d, 0.025) > 0 or np.quantile(d, 0.975) < 0
                    ),
                )
            )

# Reclamation reference: do(S −2 SD)
delta_z_recl = np.array([0.0, 0.0, -2.0])
for comp_id, comp_name in COMP_NAMES.items():
    d = delta_per_compartment(samples, comp_id, delta_z_recl)
    records.append(
        dict(
            scenario="Reclamation",
            horizon="do(S-2SD)",
            compartment=comp_name,
            pathway="intervention",
            delta_T_C=0.0,
            delta_P_pct=0.0,
            delta_shannon_z_median=float(np.median(d)),
            ci_lo=float(np.quantile(d, 0.025)),
            ci_hi=float(np.quantile(d, 0.975)),
            credible_95=bool(
                np.quantile(d, 0.025) > 0 or np.quantile(d, 0.975) < 0
            ),
        )
    )

out = pd.DataFrame(records)
out_path = CACHE_DIR / "cmip6_interventions.tsv"
out.to_csv(out_path, sep="\t", index=False)
print(f"\nwrote {out_path}")
print(out.to_string(index=False))

# ---------------------------------------------------------------------
#  4. Figure
# ---------------------------------------------------------------------
COMP_ORDER = ["surface", "deep", "rhizosphere"]
SCEN_ORDER = ["SSP1-2.6", "SSP2-4.5", "SSP3-7.0"]
HORIZONS = ["2050", "2100"]
COLORS = {"SSP1-2.6": "#3b8bba", "SSP2-4.5": "#e1a13b", "SSP3-7.0": "#c0392b"}

fig, axes = plt.subplots(1, 3, figsize=(11.5, 4.0), sharey=True)
plt.rcParams.update({"font.size": 9})

x_positions = {}
xtick_labels = []
xtick_pos = []
i = 0
for sce in SCEN_ORDER:
    for hor in HORIZONS:
        x_positions[(sce, hor)] = i
        xtick_labels.append(f"{sce}\n{hor}")
        xtick_pos.append(i)
        i += 1
i_recl = i + 0.5
x_positions[("Reclamation", "do(S-2SD)")] = i_recl
xtick_pos.append(i_recl)
xtick_labels.append("Reclam.\n(S−2 SD)")

for ax, comp in zip(axes, COMP_ORDER):
    sub = out[out.compartment == comp]
    for _, r in sub.iterrows():
        x = x_positions[(r.scenario, r.horizon)]
        if r.scenario == "Reclamation":
            ax.errorbar(
                x, r.delta_shannon_z_median,
                yerr=[[r.delta_shannon_z_median - r.ci_lo],
                      [r.ci_hi - r.delta_shannon_z_median]],
                fmt="D", color="#27ae60", markersize=8, capsize=4, lw=1.6,
                label="Reclamation",
            )
        else:
            offset = -0.18 if r.pathway == "direct" else +0.18
            marker = "o" if r.pathway == "direct" else "s"
            face = COLORS[r.scenario] if r.pathway == "direct" else "white"
            ax.errorbar(
                x + offset, r.delta_shannon_z_median,
                yerr=[[r.delta_shannon_z_median - r.ci_lo],
                      [r.ci_hi - r.delta_shannon_z_median]],
                fmt=marker, color=COLORS[r.scenario],
                markerfacecolor=face, markersize=6, capsize=3, lw=1.2,
            )
    ax.axhline(0, color="k", lw=0.5)
    ax.set_xticks(xtick_pos)
    ax.set_xticklabels(xtick_labels, fontsize=7.5)
    ax.set_title(comp.capitalize(), fontsize=10)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    if comp == "surface":
        ax.set_ylabel(r"$\Delta$ Shannon ($\sigma$ units, posterior median $\pm$ 95% CI)")

# legend on the rightmost axis
import matplotlib.lines as mlines
legend_handles = [
    mlines.Line2D([], [], marker="o", color=COLORS["SSP1-2.6"],
                  linestyle="None", markersize=6, label="SSP1-2.6 (direct)"),
    mlines.Line2D([], [], marker="o", color=COLORS["SSP2-4.5"],
                  linestyle="None", markersize=6, label="SSP2-4.5 (direct)"),
    mlines.Line2D([], [], marker="o", color=COLORS["SSP3-7.0"],
                  linestyle="None", markersize=6, label="SSP3-7.0 (direct)"),
    mlines.Line2D([], [], marker="s", color="grey", markerfacecolor="white",
                  linestyle="None", markersize=6,
                  label="+ salinity amplification"),
    mlines.Line2D([], [], marker="D", color="#27ae60",
                  linestyle="None", markersize=8, label="do(S−2 SD) reclamation"),
]
axes[-1].legend(handles=legend_handles, loc="upper right", frameon=False, fontsize=7.5)

fig.suptitle(
    r"CMIP6-driven projections of community Shannon: "
    r"climate forcing vs reclamation intervention",
    fontsize=11,
)
plt.tight_layout(rect=(0, 0, 1, 0.95))
fig_path = FIGURES_DIR / "fig_cmip6_projections.pdf"
plt.savefig(fig_path, bbox_inches="tight")
print(f"wrote {fig_path}")
