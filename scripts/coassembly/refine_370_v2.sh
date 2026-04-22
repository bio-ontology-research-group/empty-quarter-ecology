#!/usr/bin/env bash
#SBATCH --job-name=csp370-refine2
#SBATCH --output=/data/emptyquarter/ecology-paper-runs/logs/csp370_refine2.%j.out
#SBATCH --error=/data/emptyquarter/ecology-paper-runs/logs/csp370_refine2.%j.err
#SBATCH --partition=debug
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=4:00:00
#SBATCH --exclude=node001,node004,node007

set -u -o pipefail
source /storage/miniforge3/etc/profile.d/conda.sh
conda activate metagenomics

WD=/data/emptyquarter/ecology-paper-runs/csp_mag/refine_370_v2
mkdir -p $WD && cd $WD
BIN=/data/emptyquarter/ecology-paper-runs/csp_mag/bins/coassembly_clean/SemiBin_370.fa

echo "=== contig stats (GC + length) ==="
python3 << 'PYEOF'
from pathlib import Path
bin_fa = Path("/data/emptyquarter/ecology-paper-runs/csp_mag/bins/coassembly_clean/SemiBin_370.fa")
contigs = {}
name, seq = None, []
with open(bin_fa) as fh:
    for line in fh:
        if line.startswith('>'):
            if name: contigs[name] = ''.join(seq)
            name = line[1:].split()[0].strip(); seq = []
        else:
            seq.append(line.strip())
if name: contigs[name] = ''.join(seq)
print(f"contigs: {len(contigs)}")
with open('contig_stats.tsv', 'w') as out:
    out.write("contig\tlength\tgc\n")
    for k, s in contigs.items():
        L = len(s)
        gc = (s.upper().count('G') + s.upper().count('C')) / max(1, L)
        out.write(f"{k}\t{L}\t{gc:.4f}\n")
PYEOF

# Per-sample coverage using samtools coverage (no --region-file)
BAMS=(/data/emptyquarter/ecology-paper-runs/csp_mag/map/V26Dr2.sorted.bam \
      /data/emptyquarter/ecology-paper-runs/csp_mag/map/V30Sr2.sorted.bam \
      /data/emptyquarter/ecology-paper-runs/csp_mag/map/V32PRr1.sorted.bam \
      /data/emptyquarter/ecology-paper-runs/csp_mag/map/V38PRr3.sorted.bam \
      /data/emptyquarter/ecology-paper-runs/csp_mag/map/V39Sr2.sorted.bam)

# Extract SemiBin_370 contig names
cut -f1 contig_stats.tsv | tail -n +2 | sort -u > bin_contigs.txt
echo "SemiBin_370 has $(wc -l < bin_contigs.txt) contigs"

echo "=== compute coverage per contig per sample (filter bin's contigs) ==="
for BAM in "${BAMS[@]}"; do
  SAMPLE=$(basename $BAM .sorted.bam)
  OUT=${SAMPLE}.cov.tsv
  if [ ! -s "$OUT" ]; then
    samtools coverage -H -o /tmp/${SAMPLE}.cov.full.tsv "$BAM"
    # Filter to bin's contigs
    awk -F'\t' 'NR==FNR {c[$1]=1; next} $1 in c' bin_contigs.txt /tmp/${SAMPLE}.cov.full.tsv > $OUT
    rm -f /tmp/${SAMPLE}.cov.full.tsv
  fi
  n=$(wc -l < $OUT)
  echo "  $SAMPLE: $n contig coverage rows"
done

python3 << 'PYEOF'
import pandas as pd, numpy as np
from pathlib import Path

samples = ["V26Dr2","V30Sr2","V32PRr1","V38PRr3","V39Sr2"]
stats = pd.read_csv("contig_stats.tsv", sep="\t")
merged = stats.copy()
# samtools coverage has columns: rname, startpos, endpos, numreads, covbases, coverage, meandepth, meanbaseq, meanmapq
for s in samples:
    f = f"{s}.cov.tsv"
    if not Path(f).exists(): continue
    df = pd.read_csv(f, sep="\t", header=None,
                     names=["contig","startpos","endpos","numreads","covbases","coverage","meandepth","meanbaseq","meanmapq"])
    merged = merged.merge(df[["contig","meandepth"]].rename(columns={"meandepth": s}), on="contig", how="left")

cov_cols = [c for c in samples if c in merged.columns]
merged["cov_mean"] = merged[cov_cols].mean(axis=1)
merged["cov_log"] = np.log10(merged.cov_mean.clip(lower=0.1))
merged.to_csv("merged_stats.tsv", sep="\t", index=False)
print(f"merged: {merged.shape}, total size: {merged.length.sum()/1e6:.2f} Mb")
print(f"cov_mean median: {merged.cov_mean.median():.2f}, p10: {merged.cov_mean.quantile(0.1):.2f}, p90: {merged.cov_mean.quantile(0.9):.2f}")
print(f"GC median: {merged.gc.median():.3f}, std: {merged.gc.std():.3f}")

# Outlier flagging
gc_med = merged.gc.median()
gc_mad = 1.4826 * np.median(np.abs(merged.gc - gc_med))  # MAD-based
gc_z = (merged.gc - gc_med) / max(gc_mad, 0.005)
outlier_gc = abs(gc_z) > 3  # more lenient than T-score

cov_log_med = merged.cov_log.median()
cov_log_mad = 1.4826 * np.median(np.abs(merged.cov_log - cov_log_med))
cov_z = (merged.cov_log - cov_log_med) / max(cov_log_mad, 0.1)
outlier_cov = abs(cov_z) > 3

# Very short contigs are suspect: remove <1500 bp
short = merged.length < 1500

merged["outlier"] = outlier_gc | outlier_cov | short
n_out = merged.outlier.sum()
kept_size = merged[~merged.outlier].length.sum()
print(f"\nOutlier contigs: {n_out} ({merged[merged.outlier].length.sum()/1e3:.0f} kb)")
print(f"  outlier_gc:    {outlier_gc.sum()} ({merged[outlier_gc].length.sum()/1e3:.0f} kb)")
print(f"  outlier_cov:   {outlier_cov.sum()} ({merged[outlier_cov].length.sum()/1e3:.0f} kb)")
print(f"  short (<1500): {short.sum()} ({merged[short].length.sum()/1e3:.0f} kb)")
print(f"Kept: {len(merged)-n_out} contigs, {kept_size/1e6:.2f} Mb")

merged[merged.outlier].to_csv("removed_contigs.tsv", sep="\t", index=False)
keep = merged[~merged.outlier].contig.tolist()
with open("keep_contigs.txt","w") as fh:
    for k in keep: fh.write(k+"\n")
PYEOF

# Write refined FASTA
python3 << 'PYEOF'
keep = set(open("keep_contigs.txt").read().splitlines())
w = False
with open("/data/emptyquarter/ecology-paper-runs/csp_mag/bins/coassembly_clean/SemiBin_370.fa") as fh, \
     open("SemiBin_370_refined.fa","w") as out:
    for line in fh:
        if line.startswith(">"):
            name = line[1:].split()[0].strip()
            w = name in keep
        if w: out.write(line)
PYEOF
echo "refined bin size: $(wc -c < SemiBin_370_refined.fa) bytes, contigs: $(grep -c '^>' SemiBin_370_refined.fa)"

echo
echo "=== CheckM2 on refined bin ==="
mkdir -p refined_input
cp SemiBin_370_refined.fa refined_input/SemiBin_370_refined.fa
conda activate checkm2-v3
CHECKM2=/storage/miniforge3/envs/checkm2-v3/bin/checkm2
CHECKM_DB=/data/emptyquarter/ecology-paper-runs/checkm2_db/uniref100.KO.1.dmnd
mkdir -p checkm2_refined
$CHECKM2 predict --input refined_input -x fa --output-directory checkm2_refined \
  --threads 8 --database_path $CHECKM_DB --force 2>&1 | tail -10

echo
echo "=== Before vs After ==="
echo "Before (SemiBin_370):"
awk '$1=="SemiBin_370" {print "C="$2"%  X="$3"%  size="$9"  N50="$7}' /data/emptyquarter/ecology-paper-runs/csp_mag/quality/checkm2_coassembly/quality_report.tsv
echo "After (refined):"
awk 'NR>1 {print "C="$2"%  X="$3"%  size="$9"  N50="$7}' checkm2_refined/quality_report.tsv

echo "=== DONE ==="
