#!/bin/bash
#SBATCH --job-name=gem_dl
#SBATCH --partition=debug
#SBATCH --exclude=node001,node004,node006,node007
#SBATCH --cpus-per-task=16
#SBATCH --mem=16G
#SBATCH --time=08:00:00
#SBATCH --output=/data/emptyquarter/ecology-paper-runs/public_metagenomes/logs/download_%j.out
#SBATCH --error=/data/emptyquarter/ecology-paper-runs/public_metagenomes/logs/download_%j.err
set -euo pipefail
source /storage/miniforge3/etc/profile.d/conda.sh
conda activate metagenomics
WORK=/data/emptyquarter/ecology-paper-runs/public_metagenomes
MAGS=$WORK/mags
mkdir -p $MAGS
cd $WORK
awk -F'\t' 'NR>1 {print $1}' metadata/selected_mags.tsv > metadata/_ids.txt
TOTAL=$(wc -l < metadata/_ids.txt)
echo "[start] $(date) Downloading $TOTAL MAGs to $MAGS"
awk '{print "https://portal.nersc.gov/GEM/genomes/fna/"$1".fna.gz\n  out="$1".fna.gz"}' metadata/_ids.txt > metadata/_urls.txt
aria2c -i metadata/_urls.txt -d $MAGS -j 16 -x 4 -c --retry-wait=5 --max-tries=3 --summary-interval=30 || true
echo "[mid] $(date) Downloaded so far:"
ls $MAGS | wc -l
while read id; do
  if [ ! -s $MAGS/$id.fna.gz ]; then
    curl -sSL -o $MAGS/$id.fna.gz "https://portal.nersc.gov/GEM/genomes/fna/$id.fna.gz" || echo "FAIL $id" >> logs/download_failed.txt
  fi
done < metadata/_ids.txt
echo "[end] $(date) Final count:"
ls $MAGS | wc -l
