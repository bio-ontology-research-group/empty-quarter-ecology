#!/bin/bash
#SBATCH --job-name=betA_tree
#SBATCH --partition=debug
#SBATCH --exclude=node001,node004,node006,node007
#SBATCH --cpus-per-task=16
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --output=/data/emptyquarter/ecology-paper-runs/public_metagenomes/logs/msatree_%j.out
#SBATCH --error=/data/emptyquarter/ecology-paper-runs/public_metagenomes/logs/msatree_%j.err
set -euo pipefail
source /storage/miniforge3/etc/profile.d/conda.sh
conda activate metagenomics
MSA=/data/emptyquarter/ecology-paper-runs/public_metagenomes/msa
cd $MSA
echo [MAFFT] $(date)
mafft --auto --thread 16 betA_combined.faa > betA_combined.aln.faa 2> logs_mafft.err
echo [trimAl] $(date)
trimal -in betA_combined.aln.faa -out betA_combined.trim.faa -automated1 > logs_trimal.out
grep -c '^>' betA_combined.trim.faa
awk 'BEGIN{l=0} /^>/ {next} {l+=length($0)} END {print "Alignment length:", l/NR"*nseq"}' betA_combined.trim.faa
echo [IQ-TREE] $(date)
iqtree -s betA_combined.trim.faa -m MFP -B 1000 -alrt 1000 -T 16 -pre betA_tree --seqtype AA --redo 2>&1 | tail -20
echo [DONE] $(date)
ls -la betA_tree.*
