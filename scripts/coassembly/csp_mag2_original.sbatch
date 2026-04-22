#!/bin/bash
#SBATCH --job-name=eq-csp2
#SBATCH --output=/data/emptyquarter/ecology-paper-runs/logs/csp_mag2.%j.out
#SBATCH --error=/data/emptyquarter/ecology-paper-runs/logs/csp_mag2.%j.err
#SBATCH --partition=debug
#SBATCH --cpus-per-task=16
#SBATCH --mem=120G
#SBATCH --time=96:00:00
#SBATCH --exclude=node001,node007

set -euo pipefail
export HOME=/tmp
CORES=$SLURM_CPUS_PER_TASK

ROOT=/data/emptyquarter/ecology-paper-runs/csp_mag
mkdir -p "$ROOT"/{reads_extracted,assembly,map,bins,quality,taxonomy,logs}

SAMPLES=(V26Dr2 V30Sr2 V39Sr2 V32PRr1 V38PRr3)
READS_BASE=/data/emptyquarter/sequencing-results/sites/trip5

source /storage/miniforge3/etc/profile.d/conda.sh
conda activate metagenomics

echo "=== [$(date)] 0. Extract paired reads from BAMs ==="
for s in "${SAMPLES[@]}"; do
  site_num=$(echo "$s" | sed -E "s/V([0-9]+).*/\1/")
  bam="${READS_BASE}/site${site_num}/${s}/wgs/reads/${s}.sorted.bam"
  r1="$ROOT/reads_extracted/${s}_R1.fastq.gz"
  r2="$ROOT/reads_extracted/${s}_R2.fastq.gz"
  if [ -f "$r1" ] && [ -f "$r2" ]; then
    echo "  $s: already extracted"; continue
  fi
  echo "  extracting $s from $bam"
  samtools collate -O -@ $CORES "$bam" | \
    samtools fastq -@ $CORES -1 "$r1" -2 "$r2" -0 /dev/null -s /dev/null -n -
  ls -lh "$r1" "$r2"
done

echo "=== [$(date)] 1. MEGAHIT co-assembly ==="
ASM="$ROOT/assembly/coassembly"
if [ ! -f "$ASM/final.contigs.fa" ]; then
  R1_CSV=$(printf "$ROOT/reads_extracted/%s_R1.fastq.gz," "${SAMPLES[@]}"); R1_CSV=${R1_CSV%,}
  R2_CSV=$(printf "$ROOT/reads_extracted/%s_R2.fastq.gz," "${SAMPLES[@]}"); R2_CSV=${R2_CSV%,}
  rm -rf "$ASM"
  megahit -1 "$R1_CSV" -2 "$R2_CSV" -o "$ASM" \
    -t $CORES --min-contig-len 1500 --presets meta-large \
    2>&1 | tail -10
fi
CONTIGS="$ASM/final.contigs.fa"
echo "contigs: $(grep -c '^>' $CONTIGS)"

echo "=== [$(date)] 2. Map reads back ==="
if [ ! -f "${CONTIGS}.bwt" ]; then bwa index "$CONTIGS"; fi
for s in "${SAMPLES[@]}"; do
  bam="$ROOT/map/${s}.sorted.bam"
  [ -f "$bam" ] && continue
  bwa mem -t $CORES "$CONTIGS" "$ROOT/reads_extracted/${s}_R1.fastq.gz" "$ROOT/reads_extracted/${s}_R2.fastq.gz" 2>/dev/null \
    | samtools sort -@ $CORES -o "$bam" -
  samtools index "$bam"
done

echo "=== [$(date)] 3. SemiBin2 multi-sample binning ==="
conda activate semibin2-env
BIN_OUT=$ROOT/bins
rm -rf "$BIN_OUT/semibin_run"
SemiBin2 single_easy_bin \
  --input-fasta "$CONTIGS" \
  --input-bam $ROOT/map/*.sorted.bam \
  --environment global \
  -o "$BIN_OUT/semibin_run" \
  --threads $CORES 2>&1 | tail -10
ls -lh "$BIN_OUT/semibin_run/output_bins/" | head -20

echo "=== [$(date)] 4. CheckM2 quality ==="
conda activate checkm2-v3
CHECKM2_DB=$(find /storage -maxdepth 5 -name "uniref100.KO.1.dmnd" -type f 2>/dev/null | head -1)
mkdir -p "$ROOT/bins/gz2fa"
for f in "$BIN_OUT/semibin_run/output_bins/"*.fa.gz; do
  bn=$(basename "$f" .fa.gz)
  zcat "$f" > "$ROOT/bins/gz2fa/${bn}.fa"
done
if [ -n "$CHECKM2_DB" ]; then
  checkm2 predict --input "$ROOT/bins/gz2fa/" \
    --output-directory "$ROOT/quality/checkm2" \
    --threads $CORES --extension fa \
    --database_path "$CHECKM2_DB" 2>&1 | tail -5
  head -20 "$ROOT/quality/checkm2/quality_report.tsv"
fi

echo "=== [$(date)] 5. GTDB-Tk ==="
conda activate metagenomics
export GTDBTK_DATA_PATH=/storage/software/databases/gtdbtk
gtdbtk classify_wf \
  --genome_dir "$ROOT/bins/gz2fa" \
  --out_dir "$ROOT/taxonomy/gtdbtk" \
  --extension fa --skip_ani_screen --force \
  --cpus $CORES --pplacer_cpus $CORES 2>&1 | tail -8

echo "=== [$(date)] 6. Dadabacteria/CSP1-2 hits ==="
grep -iE "Dadabacteri|CSP1-2" "$ROOT/taxonomy/gtdbtk/"*.summary.tsv 2>/dev/null || echo "(none)"
grep -iE "Dadabacteri|CSP1-2" "$ROOT/taxonomy/gtdbtk/classify/"*.summary.tsv 2>/dev/null || echo "(none in classify)"

echo "=== [$(date)] done ==="
