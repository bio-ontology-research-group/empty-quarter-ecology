#!/usr/bin/env bash
#SBATCH --job-name=eq-gurb-cd
#SBATCH --partition=debug
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=12:00:00
#SBATCH --exclude=node001,node004,node007
#SBATCH --output=/data/emptyquarter/ecology-paper-runs/crossdesert/logs/gurb_%j.log
#SBATCH --error=/data/emptyquarter/ecology-paper-runs/crossdesert/logs/gurb_%j.err

set -eu -o pipefail
WD=/data/emptyquarter/ecology-paper-runs/crossdesert
cd "$WD"
mkdir -p logs raw/Gurbantunggut processed/trimmed/Gurbantunggut processed/merged/Gurbantunggut processed/derep/Gurbantunggut processed/otus processed/tables

echo "[gurb] start: $(date)"

# -------- Stage 1: download (idempotent) --------
awk -F'\t' 'NR>1 {print $1"\t"$11}' gurb_manifest.tsv | while IFS=$'\t' read -r run ftp; do
  if [ -f raw/Gurbantunggut/${run}.done ]; then continue; fi
  IFS=';' read -r -a urls <<< "$ftp"
  ok=1
  for url in "${urls[@]}"; do
    [ -z "$url" ] && continue
    target=raw/Gurbantunggut/$(basename "$url")
    if [ -s "$target" ]; then continue; fi
    if ! [[ "$url" == http* ]]; then url="https://$url"; fi
    if ! wget -q --tries=3 --timeout=180 -O "$target" "$url"; then
      rm -f "$target"; ok=0; break
    fi
  done
  if [ $ok -eq 1 ]; then touch raw/Gurbantunggut/${run}.done; fi
done
DL=$(ls raw/Gurbantunggut/*.done 2>/dev/null | wc -l)
echo "[gurb] downloaded: $DL"

# -------- Stage 2: SKIP cutadapt (reads are pre-trimmed per sample_alias 'trim.1.fq') --------
# Just quality-filter and derep
source /storage/miniforge3/etc/profile.d/conda.sh
conda activate metagenomics

awk -F'\t' 'NR>1 {print $1}' gurb_manifest.tsv | while read -r run; do
  [ -f raw/Gurbantunggut/${run}.done ] || continue
  done_marker=processed/derep/Gurbantunggut/${run}.done
  [ -f "$done_marker" ] && continue

  se=raw/Gurbantunggut/${run}.fastq.gz
  [ -f "$se" ] || { echo "[gurb] SKIP $run no file"; continue; }

  outm=processed/merged/Gurbantunggut
  outd=processed/derep/Gurbantunggut

  # Permissive quality filter: maxee 2.0 for single-end pre-trimmed reads
  vsearch --fastq_filter "$se" --fastq_maxee 2.0 --fastq_minlen 200 --fastq_maxlen 300 \
          --fastaout $outm/${run}.fa --threads 4 2>/dev/null \
          || { echo "[gurb] FAIL filter $run"; continue; }

  n_filt=$(grep -c '^>' $outm/${run}.fa 2>/dev/null || echo 0)
  if [ "$n_filt" -lt 100 ]; then
    echo "[gurb] SKIP $run too few filtered reads ($n_filt)"
    continue
  fi

  vsearch --derep_fulllength $outm/${run}.fa --output $outd/${run}.derep.fa \
          --sizein --sizeout --minuniquesize 2 --relabel "${run}_" --threads 4 >/dev/null 2>&1 \
          || { echo "[gurb] FAIL derep $run"; continue; }
  touch "$done_marker"
done
DD=$(ls processed/derep/Gurbantunggut/*.done 2>/dev/null | wc -l)
echo "[gurb] derep done: $DD"

# -------- Per-desert OTU build --------
cat processed/derep/Gurbantunggut/*.derep.fa 2>/dev/null > processed/otus/Gurbantunggut_all.fa
N_ALL=$(grep -c '^>' processed/otus/Gurbantunggut_all.fa 2>/dev/null || echo 0)
echo "[gurb] concatenated derep: $N_ALL sequences"

if [ "$N_ALL" -lt 100 ]; then
  echo "[gurb] FATAL: not enough derep sequences; check filter settings"
  exit 1
fi

vsearch --derep_fulllength processed/otus/Gurbantunggut_all.fa \
        --output processed/otus/Gurbantunggut_global_derep.fa \
        --sizein --sizeout --minuniquesize 2 --threads 16 >/dev/null 2>&1
vsearch --uchime3_denovo processed/otus/Gurbantunggut_global_derep.fa \
        --nonchimeras processed/otus/Gurbantunggut_nochim.fa \
        --sizein --sizeout --threads 16 >/dev/null 2>&1
vsearch --cluster_size processed/otus/Gurbantunggut_nochim.fa --id 0.97 \
        --centroids processed/otus/Gurbantunggut_otus.fa --sizein --sizeout --threads 16 >/dev/null 2>&1
n_otus=$(grep -c '^>' processed/otus/Gurbantunggut_otus.fa || echo 0)
echo "[gurb] OTUs: $n_otus"

# -------- Per-sample OTU table --------
cat processed/derep/Gurbantunggut/*.derep.fa > processed/tables/Gurbantunggut_all_forotu.fa
vsearch --usearch_global processed/tables/Gurbantunggut_all_forotu.fa \
        --db processed/otus/Gurbantunggut_otus.fa --id 0.97 --strand plus \
        --otutabout processed/tables/Gurbantunggut_otutab.tsv --threads 16 2>/dev/null
echo "[gurb] OTU table rows (incl header): $(wc -l < processed/tables/Gurbantunggut_otutab.tsv)"

# -------- CSP1-2 search --------
for PID in 0.97 0.90 0.85 0.75; do
  vsearch --usearch_global processed/otus/Gurbantunggut_otus.fa \
          --db csp1-2_asvs.fasta --id $PID --strand both \
          --blast6out processed/tables/Gurbantunggut_csp_hits_id${PID}.tsv \
          --threads 16 2>/dev/null || true
  H=$(wc -l < processed/tables/Gurbantunggut_csp_hits_id${PID}.tsv 2>/dev/null || echo 0)
  echo "[gurb] CSP1-2 OTU hits at $PID: $H"
done

echo "[gurb] end: $(date)"
