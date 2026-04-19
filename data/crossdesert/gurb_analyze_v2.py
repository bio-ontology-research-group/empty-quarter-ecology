#!/usr/bin/env python3
"""Gurbantunggut analysis v2 — matches stage3 sample-aware OTU table approach.
Computes per-sample Shannon & CSP1-2 abundance, no pandas load of 6GB table.
"""
import sys, re, subprocess, os
from pathlib import Path
from collections import defaultdict
import numpy as np
from scipy import stats

WD = Path('/data/emptyquarter/ecology-paper-runs/crossdesert')
VSEARCH = '/storage/miniforge3/envs/metagenomics/bin/vsearch'

otus = WD / 'processed/otus/Gurbantunggut_otus.fa'
csp_ref = WD / 'csp1-2_asvs.fasta'
derep_dir = WD / 'processed/derep/Gurbantunggut'

assert otus.exists()

# Pool per-sample derep files into one file (with labels preserved)
pooled = WD / 'processed/tables/Gurbantunggut_pooled.fa'
hits_b6 = WD / 'processed/tables/Gurbantunggut_hits.b6'

if not pooled.exists() or pooled.stat().st_size == 0:
    print("Pooling derep files...")
    with open(pooled, 'w') as out:
        for fa in sorted(derep_dir.glob('*.derep.fa')):
            with open(fa) as fh:
                out.write(fh.read())

# Map pooled reads back to OTUs via usearch_global
if not hits_b6.exists() or hits_b6.stat().st_size == 0:
    print("Running vsearch usearch_global pooled reads -> OTUs...")
    subprocess.run([
        VSEARCH, '--usearch_global', str(pooled),
        '--db', str(otus), '--id', '0.97', '--strand', 'both', '--threads', '8',
        '--blast6out', str(hits_b6), '--top_hits_only'
    ], check=True)
print("Mapping complete")

# Build counts[otu][sample] = size
counts = defaultdict(lambda: defaultdict(int))
samples_seen = set()
with open(hits_b6) as fh:
    for line in fh:
        f = line.rstrip('\n').split('\t')
        if len(f) < 2:
            continue
        q, t = f[0], f[1]
        # query label format: "SRR5838954_4 size=N"
        # sample ID = first underscore-separated token
        run = q.split('_')[0]
        samples_seen.add(run)
        m = re.search(r'size=(\d+)', q)
        size = int(m.group(1)) if m else 1
        counts[t][run] += size
print(f"Samples seen in hits: {len(samples_seen)}")

def shannon(counts_dict):
    vals = np.array([v for v in counts_dict.values() if v > 0], dtype=float)
    if vals.sum() == 0:
        return float('nan')
    p = vals / vals.sum()
    return float(-(p * np.log(p)).sum())

# Invert: per-sample OTU counts
per_sample_counts = defaultdict(dict)
for otu, sample_counts in counts.items():
    for run, c in sample_counts.items():
        per_sample_counts[run][otu] = c

# Per-sample Shannon
sh_per = {run: shannon(d) for run, d in per_sample_counts.items()}

# CSP1-2 search per sample against derep files
print("CSP1-2 search per sample...")
def csp_search(derep_fa, pid):
    r = subprocess.run([
        VSEARCH, '--usearch_global', str(derep_fa), '--db', str(csp_ref),
        '--id', str(pid), '--strand', 'both', '--threads', '2',
        '--blast6out', '/dev/stdout'
    ], capture_output=True, text=True)
    abund = 0
    for line in r.stdout.splitlines():
        f = line.split('\t')
        if len(f) < 2:
            continue
        m = re.search(r'size=(\d+)', f[0])
        abund += int(m.group(1)) if m else 1
    return abund

rows = []
derep_fas = sorted(derep_dir.glob('*.derep.fa'))
for fa in derep_fas:
    run = fa.stem.replace('.derep','')
    csp97 = csp_search(fa, 0.97)
    csp85 = csp_search(fa, 0.85)
    total = 0
    with open(fa) as fh:
        for line in fh:
            if line.startswith('>'):
                m = re.search(r'size=(\d+)', line)
                total += int(m.group(1)) if m else 1
    rows.append({
        'run_accession': run,
        'desert': 'Gurbantunggut',
        'shannon_otu97': sh_per.get(run, float('nan')),
        'total_reads_postqc': total,
        'n_otus_observed': sum(1 for v in per_sample_counts.get(run,{}).values() if v>0),
        'csp_abundance_97id': csp97,
        'csp_abundance_85id': csp85,
        'csp_rel_97': csp97/total if total else float('nan'),
        'csp_rel_85': csp85/total if total else float('nan'),
    })

# Write
import csv
cols = list(rows[0].keys())
with open(WD / 'gurbantunggut_per_sample.tsv', 'w') as fh:
    w = csv.DictWriter(fh, fieldnames=cols, delimiter='\t')
    w.writeheader()
    w.writerows(rows)

# Stats
n = len(rows)
sh_vals = np.array([r['shannon_otu97'] for r in rows if not np.isnan(r['shannon_otu97'])])
csp85_rel = np.array([r['csp_rel_85'] for r in rows])
csp97_rel = np.array([r['csp_rel_97'] for r in rows])
csp85_abs = np.array([r['csp_abundance_85id'] for r in rows])
csp97_abs = np.array([r['csp_abundance_97id'] for r in rows])

n_csp85 = (csp85_abs > 0).sum()
n_csp97 = (csp97_abs > 0).sum()
print(f"\nGurbantunggut summary:")
print(f"  Total samples: {n}")
print(f"  Samples with Shannon computed: {len(sh_vals)} (median {np.median(sh_vals):.2f})")
print(f"  Samples with CSP1-2 at 97% V4: {n_csp97} ({100*n_csp97/n:.0f}%)")
print(f"  Samples with CSP1-2 at 85% V4: {n_csp85} ({100*n_csp85/n:.0f}%)")
print(f"  Mean CSP rel abund 85%: {csp85_rel.mean()*100:.4f}%")
print(f"  Max  CSP rel abund 85%: {csp85_rel.max()*100:.4f}%")

# Spearman CSP×Shannon on CSP-positive subset
sh_col = np.array([r['shannon_otu97'] for r in rows])
valid = ~np.isnan(sh_col) & (csp85_abs > 0)
if valid.sum() >= 10:
    rho_pos, p_pos = stats.spearmanr(csp85_rel[valid], sh_col[valid])
    print(f"  Spearman (CSP>0 only): n={valid.sum()}, rho={rho_pos:.3f}, p={p_pos:.4g}")
else:
    rho_pos = p_pos = float('nan')
    print(f"  Insufficient CSP-positive samples for correlation")

# All-samples version
valid_all = ~np.isnan(sh_col)
rho_all, p_all = stats.spearmanr(csp85_rel[valid_all], sh_col[valid_all])
print(f"  Spearman (all samples):  n={valid_all.sum()}, rho={rho_all:.3f}, p={p_all:.4g}")

with open(WD / 'gurbantunggut_final_stats.txt', 'w') as fh:
    fh.write(f"Gurbantunggut Desert (SRP112798, 120 samples, V4 515F/806R)\n")
    fh.write(f"Total samples: {n}\n")
    fh.write(f"Samples with computed Shannon (OTU-97): {len(sh_vals)}\n")
    fh.write(f"Median Shannon: {np.median(sh_vals):.3f}\n")
    fh.write(f"CSP1-2 prevalence at 85% V4: {n_csp85}/{n} ({100*n_csp85/n:.1f}%)\n")
    fh.write(f"CSP1-2 prevalence at 97% V4: {n_csp97}/{n} ({100*n_csp97/n:.1f}%)\n")
    fh.write(f"Mean CSP rel abundance at 85%: {csp85_rel.mean()*100:.4f}%\n")
    fh.write(f"Max CSP rel abundance at 85%: {csp85_rel.max()*100:.4f}%\n")
    if not np.isnan(rho_pos):
        fh.write(f"Spearman(CSP85, Shannon)  CSP>0 subset: n={valid.sum()}, rho={rho_pos:.4f}, p={p_pos:.4g}\n")
    fh.write(f"Spearman(CSP85, Shannon)  all samples: n={valid_all.sum()}, rho={rho_all:.4f}, p={p_all:.4g}\n")
print(f"\nwrote {WD/'gurbantunggut_final_stats.txt'} and gurbantunggut_per_sample.tsv")
print("[gurb-analyze-v2] done")
