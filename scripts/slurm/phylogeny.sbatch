#!/bin/bash
#SBATCH --job-name=eq-phylo
#SBATCH --output=/data/emptyquarter/ecology-paper-runs/logs/phylogeny.%j.out
#SBATCH --error=/data/emptyquarter/ecology-paper-runs/logs/phylogeny.%j.err
#SBATCH --partition=debug
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=24:00:00

set -euo pipefail

source /storage/miniforge3/etc/profile.d/conda.sh
conda activate metagenomics

IN=/data/emptyquarter/ecology-paper-runs/inputs/ASV_seqs_filtered.fasta
OUT=/data/emptyquarter/ecology-paper-runs/phylogeny
mkdir -p "$OUT"

echo "[$(date)] MAFFT alignment on $(grep -c '^>' "$IN") sequences with ${SLURM_CPUS_PER_TASK} threads"
# --parttree + --retree 1 --maxiterate 0 is the standard large-scale mode for >10k sequences
mafft --parttree --retree 1 --maxiterate 0 --nofft \
      --thread ${SLURM_CPUS_PER_TASK} \
      "$IN" > "$OUT/ASV_aligned.fasta"

echo "[$(date)] MAFFT done. Alignment size:"
ls -lh "$OUT/ASV_aligned.fasta"

echo "[$(date)] FastTree GTR+Gamma on alignment"
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}
FastTreeMP -nt -gtr -gamma -fastest \
    "$OUT/ASV_aligned.fasta" > "$OUT/ASV_tree.nwk" 2> "$OUT/fasttree.log"

echo "[$(date)] FastTree done."
ls -lh "$OUT/ASV_tree.nwk"
echo "[$(date)] Pipeline complete."
