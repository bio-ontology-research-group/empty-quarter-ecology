#!/bin/bash
#SBATCH --job-name=pfam_gmc
#SBATCH --partition=debug
#SBATCH --exclude=node001,node004,node006,node007
#SBATCH --cpus-per-task=40
#SBATCH --mem=32G
#SBATCH --time=02:00:00
#SBATCH --output=/data/emptyquarter/ecology-paper-runs/public_metagenomes/logs/pfam_%j.out
#SBATCH --error=/data/emptyquarter/ecology-paper-runs/public_metagenomes/logs/pfam_%j.err
set -euo pipefail
source /storage/miniforge3/etc/profile.d/conda.sh
conda activate metagenomics
WORK=/data/emptyquarter/ecology-paper-runs/public_metagenomes
ORFS=$WORK/orfs
PFAM=$WORK/pfam
mkdir -p $PFAM/per_genome
cd $WORK
echo [start] $(date)  $(ls $ORFS/*.faa 2>/dev/null | wc -l) ORF files
export PFAM
run_one() {
  local faa=$1
  local id=$(basename $faa .faa)
  local out=$PFAM/per_genome/$id.tbl
  if [ -s $out ]; then return; fi
  hmmsearch --cpu 1 --cut_ga --tblout $out $PFAM/GMC_oxred.hmm $faa > /dev/null 2>&1 || echo FAIL $id >> logs/pfam_failed.txt
}
export -f run_one
ls $ORFS/*.faa | xargs -I{} -P 40 bash -c 'run_one "$@"' _ {}
echo [end] $(date) tbl files: $(ls $PFAM/per_genome/*.tbl 2>/dev/null | wc -l)
