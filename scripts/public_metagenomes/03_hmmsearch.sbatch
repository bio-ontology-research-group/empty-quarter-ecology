#!/bin/bash
#SBATCH --job-name=gem_hmm
#SBATCH --partition=debug
#SBATCH --exclude=node001,node004,node006,node007
#SBATCH --cpus-per-task=40
#SBATCH --mem=32G
#SBATCH --time=04:00:00
#SBATCH --output=/data/emptyquarter/ecology-paper-runs/public_metagenomes/logs/hmm_%j.out
#SBATCH --error=/data/emptyquarter/ecology-paper-runs/public_metagenomes/logs/hmm_%j.err
set -euo pipefail
source /storage/miniforge3/etc/profile.d/conda.sh
conda activate metagenomics
WORK=/data/emptyquarter/ecology-paper-runs/public_metagenomes
ORFS=$WORK/orfs
HMM=$WORK/hmm
HMM_FILE=/data/emptyquarter/ecology-paper-runs/t1_mechanism/kofam/target.hmm
KOLIST=/data/emptyquarter/ecology-paper-runs/t1_mechanism/kofam/ko_list.txt
KOTHR=/data/emptyquarter/ecology-paper-runs/t1_mechanism/kofam/ko_thresholds.tsv
mkdir -p $HMM $HMM/per_genome
cd $WORK
echo "[start] $(date)  $(ls $ORFS/*.faa 2>/dev/null | wc -l) ORF files"
run_one() {
  local faa=$1
  local id=$(basename $faa .faa)
  local out=$HMM/per_genome/$id.tbl
  if [ -s $out ]; then return; fi
  hmmsearch --cpu 1 --tblout $out /data/emptyquarter/ecology-paper-runs/t1_mechanism/kofam/target.hmm $faa > /dev/null 2>&1 || echo "FAIL $id" >> logs/hmm_failed.txt
}
export -f run_one
export HMM
ls $ORFS/*.faa | parallel -j 40 --no-notice run_one {}
echo "[end] $(date) tbl files:"
ls $HMM/per_genome/*.tbl 2>/dev/null | wc -l
