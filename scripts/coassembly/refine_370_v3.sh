#!/usr/bin/env bash
#SBATCH --job-name=csp370-refine3
#SBATCH --output=/data/emptyquarter/ecology-paper-runs/logs/csp370_refine3.%j.out
#SBATCH --error=/data/emptyquarter/ecology-paper-runs/logs/csp370_refine3.%j.err
#SBATCH --partition=debug
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=2:00:00
#SBATCH --exclude=node001,node004,node007

set -u -o pipefail
source /storage/miniforge3/etc/profile.d/conda.sh

WD=/data/emptyquarter/ecology-paper-runs/csp_mag/refine_370_v3
mkdir -p $WD && cd $WD

# Reuse computed coverage from v2
cp /data/emptyquarter/ecology-paper-runs/csp_mag/refine_370_v2/contig_stats.tsv .
cp /data/emptyquarter/ecology-paper-runs/csp_mag/refine_370_v2/merged_stats.tsv .

echo "=== Try multiple thresholds and also TNF-based outlier detection ==="
python3 << 'PYEOF'
import pandas as pd, numpy as np
from pathlib import Path

BIN = "/data/emptyquarter/ecology-paper-runs/csp_mag/bins/coassembly_clean/SemiBin_370.fa"

# Compute tetranucleotide frequency per contig
def tnf(seq, k=4):
    counts = {}
    seq = seq.upper()
    for i in range(len(seq)-k+1):
        kmer = seq[i:i+k]
        if 'N' in kmer: continue
        counts[kmer] = counts.get(kmer, 0) + 1
    total = sum(counts.values())
    return {k: v/total for k, v in counts.items()} if total else {}

contigs = {}
name, seq = None, []
with open(BIN) as fh:
    for line in fh:
        if line.startswith('>'):
            if name: contigs[name] = ''.join(seq)
            name = line[1:].split()[0].strip(); seq = []
        else:
            seq.append(line.strip())
if name: contigs[name] = ''.join(seq)

print(f"contigs: {len(contigs)}")

tnf_vectors = {}
all_kmers = set()
for name, s in contigs.items():
    t = tnf(s)
    tnf_vectors[name] = t
    all_kmers.update(t.keys())
kmers = sorted(all_kmers)
# Build matrix
mat = np.array([[tnf_vectors[n].get(k, 0) for k in kmers] for n in contigs])
# Normalise rows
mat = mat / (mat.sum(axis=1, keepdims=True) + 1e-9)

# Compute centroid
centroid = mat.mean(axis=0)
# Cosine distance to centroid
def cos_dist(row, c):
    num = np.dot(row, c)
    den = (np.linalg.norm(row) * np.linalg.norm(c)) + 1e-9
    return 1 - num/den

cdist = np.array([cos_dist(r, centroid) for r in mat])
print(f"TNF cosine-dist to centroid: median={np.median(cdist):.4f}, p95={np.percentile(cdist, 95):.4f}")

stats = pd.read_csv("merged_stats.tsv", sep="\t")
stats = stats.set_index("contig").reindex(list(contigs.keys()))
stats["tnf_cdist"] = cdist
stats.reset_index().to_csv("stats_with_tnf.tsv", sep="\t", index=False)

# Try different refinement strategies
strategies = [
    # (name, gc_z_max, cov_z_max, min_len, tnf_p_cut)
    ("v3a_gc25_len2k",       2.5, 3.0, 2000, None),
    ("v3b_gc20_len2k",       2.0, 3.0, 2000, None),
    ("v3c_gc25_cov25_len2k", 2.5, 2.5, 2000, None),
    ("v3d_gc25_tnf95",       2.5, 3.0, 1500, 0.95),
    ("v3e_gc25_tnf90",       2.5, 3.0, 1500, 0.90),
    ("v3f_gc20_cov20_tnf90_len3k", 2.0, 2.0, 3000, 0.90),
]

# Use MAD
gc_med = stats.gc.median()
gc_mad = 1.4826 * np.median(np.abs(stats.gc - gc_med))
cov_log_med = stats.cov_log.median()
cov_log_mad = 1.4826 * np.median(np.abs(stats.cov_log - cov_log_med))

results = []
for name, gc_z, cov_z, min_len, tnf_p in strategies:
    gc_outlier = abs((stats.gc - gc_med) / max(gc_mad, 0.005)) > gc_z
    cov_outlier = abs((stats.cov_log - cov_log_med) / max(cov_log_mad, 0.1)) > cov_z
    short = stats.length < min_len
    if tnf_p is not None:
        tnf_cut = stats.tnf_cdist.quantile(tnf_p)
        tnf_outlier = stats.tnf_cdist > tnf_cut
    else:
        tnf_outlier = pd.Series(False, index=stats.index)
    drop = gc_outlier | cov_outlier | short | tnf_outlier
    keep = stats[~drop].index.tolist()
    kept_size = stats[~drop].length.sum()
    with open(f"keep_{name}.txt","w") as fh:
        for k in keep: fh.write(k+"\n")
    results.append((name, len(keep), kept_size, drop.sum(), stats[drop].length.sum()))
    print(f"  {name}: kept {len(keep)} contigs ({kept_size/1e6:.2f} Mb), dropped {drop.sum()} ({stats[drop].length.sum()/1e3:.0f} kb)")

# Write refined FASTAs for each strategy
for name, _, _, _, _ in results:
    keep = set(open(f"keep_{name}.txt").read().splitlines())
    with open(BIN) as fh, open(f"refined_{name}.fa","w") as out:
        w = False
        for line in fh:
            if line.startswith(">"):
                nm = line[1:].split()[0].strip()
                w = nm in keep
            if w: out.write(line)
    print(f"  wrote refined_{name}.fa")
PYEOF

# CheckM2 on all refined bins at once
echo
echo "=== CheckM2 on all refined variants ==="
mkdir -p checkm2_in
for f in refined_*.fa; do
  cp $f checkm2_in/$f
done
conda activate checkm2-v3
CHECKM2=/storage/miniforge3/envs/checkm2-v3/bin/checkm2
CHECKM_DB=/data/emptyquarter/ecology-paper-runs/checkm2_db/uniref100.KO.1.dmnd
mkdir -p checkm2_refined
$CHECKM2 predict --input checkm2_in -x fa --output-directory checkm2_refined \
  --threads 8 --database_path $CHECKM_DB --force 2>&1 | tail -5

echo
echo "=== Comparison ==="
printf "%-40s %8s %8s %10s\n" "variant" "C%" "X%" "size"
awk 'NR>1 {printf "%-40s %8.2f %8.2f %10d\n", $1, $2, $3, $9}' checkm2_refined/quality_report.tsv | sort -k2,2 -rn

echo
echo "=== Before (original SemiBin_370): C=93.37% X=6.53% size=4,082,165 ==="
echo "=== v2 refinement:                  C=93.14% X=5.15% size=3,991,002 ==="
echo "=== DONE ==="
