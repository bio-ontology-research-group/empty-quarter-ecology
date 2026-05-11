#!/bin/bash
#SBATCH --job-name=gem_prodigal
#SBATCH --partition=debug
#SBATCH --exclude=node001,node004,node006,node007
#SBATCH --cpus-per-task=40
#SBATCH --mem=32G
#SBATCH --time=08:00:00
#SBATCH --output=/data/emptyquarter/ecology-paper-runs/public_metagenomes/logs/prodigal_%j.out
#SBATCH --error=/data/emptyquarter/ecology-paper-runs/public_metagenomes/logs/prodigal_%j.err
set -euo pipefail
source /storage/miniforge3/etc/profile.d/conda.sh
conda activate metagenomics
WORK=/data/emptyquarter/ecology-paper-runs/public_metagenomes
MAGS=$WORK/mags
ORFS=$WORK/orfs
TMP=$WORK/mags_decompressed
mkdir -p $ORFS $TMP
cd $WORK
echo "[start] $(date)  $(ls $MAGS | wc -l) MAGs to process"
run_one() {
  local gz=$1
  local id=$(basename $gz .fna.gz)
  if [ -s $ORFS/$id.faa ]; then return; fi
  local fna=$TMP/$id.fna
  zcat $gz > $fna
  prodigal -p meta -q -i $fna -a $ORFS/$id.faa -o /dev/null 2>/dev/null || echo "FAIL $id" >> logs/prodigal_failed.txt
  rm -f $fna
}
export -f run_one
export ORFS TMP MAGS
ls $MAGS/*.fna.gz | parallel -j 40 --no-notice run_one {}
echo "[end] $(date) ORF files:"
ls $ORFS/*.faa 2>/dev/null | wc -l
rm -rf $TMP
