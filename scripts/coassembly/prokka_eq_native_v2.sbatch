#!/usr/bin/env bash
#SBATCH --job-name=eq-prokka2
#SBATCH --output=/data/emptyquarter/ecology-paper-runs/logs/eq_prokka2_%A_%a.log
#SBATCH --error=/data/emptyquarter/ecology-paper-runs/logs/eq_prokka2_%A_%a.err
#SBATCH --partition=debug
#SBATCH --cpus-per-task=6
#SBATCH --mem=16G
#SBATCH --time=2:00:00
#SBATCH --array=0-9
#SBATCH --exclude=node001,node004,node007

set -u -o pipefail
source /storage/miniforge3/etc/profile.d/conda.sh
conda activate prokka-env

# Shim ALL BLAST+ tools that Prokka checks — they're compiled BLAST 2.16 but Prokka's
# version comparison treats "2.16" < "2.2" lexicographically.
SHIM=/tmp/prokka_shim_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}
mkdir -p $SHIM
for tool in blastp blastn tblastn tblastx makeblastdb rpsblast psiblast; do
  REAL=$(which $tool 2>/dev/null) || continue
  cat > $SHIM/$tool << SHIM_EOF
#!/bin/bash
if [ "\${1:-}" = "-version" ] || [ "\${1:-}" = "-v" ]; then
  echo "$tool: 2.9.0+"
  echo " Package: blast 2.9.0, build Mar 30 2020"
  exit 0
fi
exec $REAL "\$@"
SHIM_EOF
  chmod +x $SHIM/$tool
done
export PATH=$SHIM:$PATH

ROOT=/data/emptyquarter/ecology-paper-runs
BINS=$ROOT/csp_mag/bins/coassembly_clean
OUT_BASE=$ROOT/csp_mag/prokka_out

declare -A MAG_LABEL=(
  ["SemiBin_370"]="coassembly_CSP12"
  ["SemiBin_120"]="Herpetosiphonaceae_BPHZ01_1"
  ["SemiBin_145"]="Herpetosiphonaceae_BPHZ01_2"
  ["SemiBin_66"]="Rubrobacter_1"
  ["SemiBin_234"]="Rubrobacter_JALYGY01"
  ["SemiBin_9"]="Rubrobacter_DATHTK01"
  ["SemiBin_236"]="Streptomyces"
  ["SemiBin_124"]="Myxococcus"
  ["SemiBin_191"]="CSP1-4_Limnocylindria"
  ["SemiBin_69"]="CSP1-6_Rokubacteriales"
)
MAG_IDS=("SemiBin_370" "SemiBin_120" "SemiBin_145" "SemiBin_66" "SemiBin_234" "SemiBin_9" "SemiBin_236" "SemiBin_124" "SemiBin_191" "SemiBin_69")

MAG=${MAG_IDS[$SLURM_ARRAY_TASK_ID]}
LABEL=${MAG_LABEL[$MAG]}
FNA=$BINS/${MAG}.fa
[ -s "$FNA" ] || { echo "Missing $FNA"; exit 1; }

OUT=$OUT_BASE/${LABEL}_${MAG}
rm -rf "$OUT"
mkdir -p $OUT_BASE
prokka --outdir "$OUT" --prefix "${LABEL}_${MAG}" \
       --cpus 6 --fast --force --addgenes "$FNA" 2>&1 | tail -20

N=$(grep -c 'CDS' $OUT/${LABEL}_${MAG}.gff 2>/dev/null || echo 0)
echo "[done] $LABEL ($MAG): $N CDS"
rm -rf $SHIM
