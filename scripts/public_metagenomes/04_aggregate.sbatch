#!/bin/bash
#SBATCH --job-name=gem_agg
#SBATCH --partition=debug
#SBATCH --nodelist=unimatrix
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=01:00:00
#SBATCH --output=/data/emptyquarter/ecology-paper-runs/public_metagenomes/logs/agg_%j.out
#SBATCH --error=/data/emptyquarter/ecology-paper-runs/public_metagenomes/logs/agg_%j.err
source /storage/miniforge3/etc/profile.d/conda.sh
conda activate metagenomics
cd /data/emptyquarter/ecology-paper-runs/public_metagenomes
python3 04_aggregate.py
