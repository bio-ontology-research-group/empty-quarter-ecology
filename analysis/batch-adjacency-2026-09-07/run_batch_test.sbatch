#!/bin/bash
#SBATCH --job-name=eq-batch-test
#SBATCH --partition=batch
#SBATCH --time=06:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=96G
#SBATCH --output=/ibex/user/hohndor/eq-batch-test-2026-09-07/slurm-%j.out
set -euo pipefail
cd /ibex/user/hohndor/eq-batch-test-2026-09-07
module load python/3.9.16
export OMP_NUM_THREADS=8
python3 batch_adjacency_test.py feature-table-trips1-5.tsv.gz batch_meta.tsv .
