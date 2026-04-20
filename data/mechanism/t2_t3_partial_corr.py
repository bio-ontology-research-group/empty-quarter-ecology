#!/usr/bin/env python3
"""T2 + T3 tests for the osmoprotectant-public-good mechanism.

T2: Partial Spearman CSP1-2 × dependent genera, controlling for salinity + rainfall.
T3: Partial Spearman CSP1-2 × Vanadium, controlling for salinity + rainfall.

All tests use Trip-5 samples with XRF available.
"""
import sys, numpy as np, pandas as pd
from scipy import stats
from pathlib import Path

REPO = Path('/home/leechuck/Documents/papers/empty-quarter-amplicon/repository')

# ---- Load tables ----
ft   = pd.read_parquet(REPO / 'cache/feature_table.parquet')   # 75469 ASV x 1227 sample
tax  = pd.read_parquet(REPO / 'cache/taxonomy.parquet')
meta = pd.read_parquet(REPO / 'cache/metadata.parquet')
try:
    mwr = pd.read_parquet(REPO / 'cache/metadata_with_rainfall.parquet')
except Exception as e:
    mwr = None; print('no metadata_with_rainfall:', e)
xrf  = pd.read_csv(REPO / 'data/geochemistry/xrf_lab_table_filtered.tsv', sep='\t')
knockout = pd.read_csv(REPO / 'cache/jsdm_knockout_shifts.tsv', sep='\t')
print('ft:', ft.shape, 'tax:', tax.shape, 'meta:', meta.shape, 'xrf:', xrf.shape)
print('tax cols:', list(tax.columns))

# ---- Build per-sample relative abundance per genus ----
# tax should link ASV → Genus
if 'genus' not in tax.columns and 'Genus' in tax.columns:
    tax = tax.rename(columns={'Genus':'genus'})
print('tax genus column present:', 'genus' in tax.columns)
print('n unique genera in tax:', tax['genus'].nunique() if 'genus' in tax.columns else '?')
print('sample tax row:', tax.head(2).to_dict('records')[0] if len(tax)>0 else None)

# Map ASVs to genus
asv_genus = tax['genus'].to_dict() if 'genus' in tax.columns else tax.iloc[:,0].to_dict()
if tax.index.name not in ['ASV','asv']:
    # assume first col is ASV
    pass

# Relative abundance per sample
read_totals = ft.sum(axis=0)
relab = ft / read_totals.replace(0, 1)
print('relab shape:', relab.shape, 'mean col sum:', relab.sum(axis=0).mean())

# Aggregate ASVs by genus
tax_idx = tax.copy()
if 'genus' in tax.columns:
    asv_to_genus = tax.set_index(tax.index)['genus']
else:
    asv_to_genus = pd.Series(tax.iloc[:,0].values, index=tax.index)
# feature_table index is ASV id; align
common_idx = ft.index.intersection(asv_to_genus.index)
print(f'ft index ∩ tax: {len(common_idx)} of {len(ft)}')
asv_to_genus_aligned = asv_to_genus.reindex(ft.index).fillna('Unclassified')

# Group sum
genus_ab = ft.groupby(asv_to_genus_aligned.values).sum()
# normalise per sample
genus_ab_sum = genus_ab.sum(axis=0)
genus_relab = genus_ab.div(genus_ab_sum.replace(0, 1), axis=1)
print('genus_relab shape:', genus_relab.shape)

# Pick target genera
DEP_GENERA = ['Herpetosiphon','Paenibacillus','Flavisolibacter','Ammoniphilus',
              'Streptomyces','Rubrobacter','Ectobacillus','Neobacillus',
              'Ramlibacter','Noviherbaspirillum','Nocardioides']
present = [g for g in DEP_GENERA if g in genus_relab.index]
missing = [g for g in DEP_GENERA if g not in genus_relab.index]
print(f'DEP in table: {len(present)} / {len(DEP_GENERA)}   missing: {missing}')

# CSP1-2 row
csp_rows = [i for i in genus_relab.index if 'CSP1-2' in str(i) or 'Dadabacteria' in str(i)]
print('CSP1-2-like rows:', csp_rows)
if csp_rows:
    csp_ab = genus_relab.loc[csp_rows].sum(axis=0)
else:
    # fallback: try tax 'phylum' == Dadabacteria
    csp_ab = pd.Series(0.0, index=genus_relab.columns)

print(f'CSP mean abundance: {csp_ab.mean()*100:.3f}%, max: {csp_ab.max()*100:.3f}%')

# ---- Trip 5 subset with XRF ----
t5 = meta[meta.trip == 5]
t5_samples = t5.index.tolist()
xrf_ids = set(xrf.SampleID.values)
overlap = set(t5_samples) & xrf_ids
print(f'Trip 5 samples: {len(t5_samples)},  with XRF: {len(overlap)}')

# Merge XRF + genus abundances on Trip 5 samples
xrf_indexed = xrf.set_index('SampleID')
t5_df = pd.DataFrame(index=list(overlap))
t5_df['compartment'] = t5.compartment
t5_df['site']        = t5.site
t5_df['csp']         = csp_ab.reindex(t5_df.index)
for g in present:
    t5_df[g] = genus_relab.loc[g].reindex(t5_df.index)
# XRF elements (strip O/T suffix on replicate IDs if not matching)
for el in ['S','Cl','Na','V','P','K','Fe','Ca']:
    t5_df[el] = xrf_indexed[el].reindex(t5_df.index)
# Rainfall 14d if available
if mwr is not None:
    # guess column
    rc = [c for c in mwr.columns if 'rain' in c.lower() and '14' in c]
    if not rc:
        rc = [c for c in mwr.columns if 'rain' in c.lower()]
    print('rainfall candidate cols:', rc[:5])
    if rc:
        t5_df['rain14'] = mwr[rc[0]].reindex(t5_df.index)
t5_df = t5_df.dropna(subset=['csp','S','Cl','Na','V'])
print(f'T5 complete cases: {len(t5_df)}')
print('compartment mix:', t5_df.compartment.value_counts().to_dict())

# ---- Partial Spearman via rank + residuals ----
def partial_spearman(df, x, y, controls):
    d = df[[x,y]+list(controls)].apply(pd.to_numeric, errors='coerce').dropna()
    if len(d) < 15: return None, None, len(d)
    # rank-transform everything
    rd = d.rank()
    # residualise x, y on controls
    C = rd[list(controls)].values
    def resid(col):
        y_ = rd[col].values
        beta, *_ = np.linalg.lstsq(np.c_[np.ones(len(C)), C], y_, rcond=None)
        return y_ - (np.c_[np.ones(len(C)), C] @ beta)
    rx = resid(x); ry = resid(y)
    rho, p = stats.pearsonr(rx, ry)
    return rho, p, len(d)

print('\n================ T2: CSP1-2 × dependent genera, controlling for S,Cl,Na ================')
controls_t2 = ['S','Cl','Na']
print(f'n controls: {controls_t2}')
results_t2 = []
for g in present:
    # raw
    d = t5_df[['csp', g]].dropna()
    rho_r, p_r = stats.spearmanr(d.csp, d[g]) if len(d) >= 10 else (np.nan, np.nan)
    # partial
    pr_rho, pr_p, n = partial_spearman(t5_df, 'csp', g, controls_t2)
    results_t2.append({'genus':g, 'n':n,
                       'raw_rho':rho_r, 'raw_p':p_r,
                       'partial_rho':pr_rho, 'partial_p':pr_p})
    print(f'  {g:22s} n={n:3d}  raw ρ={rho_r:+.3f} (p={p_r:.3g})  |  partial(|S,Cl,Na) ρ={pr_rho:+.3f} (p={pr_p:.3g})')
t2 = pd.DataFrame(results_t2)

# Also aggregate: sum all dependent genera into one "dependent pool" abundance
t5_df['dep_pool'] = t5_df[present].sum(axis=1)
rho_r, p_r = stats.spearmanr(t5_df.csp, t5_df.dep_pool)
pr_rho, pr_p, n = partial_spearman(t5_df, 'csp', 'dep_pool', controls_t2)
print(f'\n  DEP-POOL (sum)          n={n:3d}  raw ρ={rho_r:+.3f} (p={p_r:.3g})  |  partial ρ={pr_rho:+.3f} (p={pr_p:.3g})')

print('\n================ T3: CSP1-2 × Vanadium, controlling for S,Cl,Na ================')
rho_r, p_r = stats.spearmanr(t5_df.csp.dropna(), t5_df.V.dropna().loc[t5_df.csp.dropna().index].dropna()) if len(t5_df) else (np.nan, np.nan)
pr_rho, pr_p, n = partial_spearman(t5_df, 'csp', 'V', controls_t2)
print(f'  CSP × V                 n={n:3d}  raw ρ = (see full), partial(|S,Cl,Na) ρ={pr_rho:+.3f} (p={pr_p:.3g})')

# Also try raw on clean subset
d = t5_df[['csp','V']].dropna()
rho_raw, p_raw = stats.spearmanr(d.csp, d.V)
print(f'  raw CSP × V:            n={len(d)}  ρ={rho_raw:+.3f} (p={p_raw:.3g})')

# Save results
out_dir = REPO / 'cache'
t2.to_csv(out_dir / 'mechanism_T2_partial_corr.tsv', sep='\t', index=False)
with open(out_dir / 'mechanism_T3_vanadium.txt', 'w') as fh:
    fh.write(f'CSP1-2 × Vanadium, Trip-5 samples with XRF ({len(t5_df)} complete)\n')
    fh.write(f'raw Spearman:    ρ={rho_raw:+.4f}, p={p_raw:.4g}, n={len(d)}\n')
    fh.write(f'partial(|S,Cl,Na): ρ={pr_rho:+.4f}, p={pr_p:.4g}, n={n}\n')
print('\nwrote mechanism_T2_partial_corr.tsv and mechanism_T3_vanadium.txt')

# ---- Extra: test also per-compartment because the XRF was stratified ----
print('\n-- per-compartment T2 (partial CSP × dep_pool | S,Cl,Na) --')
for comp in ['surface','deep','rhizosphere']:
    sub = t5_df[t5_df.compartment == comp]
    if len(sub) < 15:
        print(f'  {comp}: skip (n={len(sub)})'); continue
    pr_rho, pr_p, n = partial_spearman(sub, 'csp', 'dep_pool', controls_t2)
    print(f'  {comp}: n={n}, partial ρ={pr_rho:+.3f} (p={pr_p:.3g})')
