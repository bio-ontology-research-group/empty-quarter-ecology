#!/bin/bash
#SBATCH --job-name=eq-bnti
#SBATCH --output=/data/emptyquarter/ecology-paper-runs/logs/bnti.%j.out
#SBATCH --error=/data/emptyquarter/ecology-paper-runs/logs/bnti.%j.err
#SBATCH --partition=debug
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=48:00:00

set -euo pipefail

ROOT=/data/emptyquarter/ecology-paper-runs
FT_TSV=$ROOT/inputs/feature_table_bnti.tsv
META_TSV=$ROOT/inputs/metadata.tsv
TREE=${TREE:-$ROOT/phylogeny/ASV_tree.nwk}
TAG=${TAG:-fasttree}
OUT=$ROOT/bnti/${TAG}_${SLURM_JOB_ID:-manual}
export R_LIBS_USER=$ROOT/envs/R-libs
export N_NULL=${N_NULL:-999}
export N_TOP=${N_TOP:-2000}
mkdir -p "$R_LIBS_USER" "$OUT" "$ROOT/bnti"

echo "[$(date)] β-NTI run (tag=$TAG)"
echo "  FT:     $FT_TSV"
echo "  META:   $META_TSV"
echo "  TREE:   $TREE"
echo "  OUT:    $OUT"
echo "  N_NULL: $N_NULL   N_TOP: $N_TOP"
echo "  CPUs:   $SLURM_CPUS_PER_TASK"

source /storage/miniforge3/etc/profile.d/conda.sh
conda activate anvio-9
which Rscript
Rscript $ROOT/inputs/bnti.R "$FT_TSV" "$TREE" "$META_TSV" "$OUT"
echo "[$(date)] done"
