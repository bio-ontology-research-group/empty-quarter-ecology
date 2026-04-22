#!/usr/bin/env bash
#SBATCH --job-name=csp-qual3
#SBATCH --output=/data/emptyquarter/ecology-paper-runs/logs/csp_qual3.%j.out
#SBATCH --error=/data/emptyquarter/ecology-paper-runs/logs/csp_qual3.%j.err
#SBATCH --partition=debug
#SBATCH --cpus-per-task=16
#SBATCH --mem=120G
#SBATCH --time=1-00:00:00
#SBATCH --exclude=node001,node004,node007

set -u -o pipefail
source /storage/miniforge3/etc/profile.d/conda.sh

ROOT=/data/emptyquarter/ecology-paper-runs/csp_mag
CLEAN_BINS=$ROOT/bins/coassembly_clean
CORES=16
CHECKM2=/storage/miniforge3/envs/checkm2-v3/bin/checkm2
GTDBTK=/storage/miniforge3/envs/metagenomics/bin/gtdbtk
CHECKM_DB=/data/emptyquarter/ecology-paper-runs/checkm2_db/uniref100.KO.1.dmnd
GTDB_DB=/storage/software/databases/gtdbtk
export GTDBTK_DATA_PATH=$GTDB_DB

NFA=$(ls $CLEAN_BINS/*.fa 2>/dev/null | wc -l)
echo "input bins: $NFA"
[ $NFA -lt 1 ] && { echo "No bins in $CLEAN_BINS; aborting"; exit 1; }

echo
echo "=== [$(date)] 4. CheckM2 (using /data/emptyquarter/.../checkm2_db) ==="
conda activate checkm2-v3
rm -rf $ROOT/quality/checkm2_coassembly
$CHECKM2 predict --input $CLEAN_BINS -x fa \
  --output-directory $ROOT/quality/checkm2_coassembly \
  --threads $CORES --database_path "$CHECKM_DB" --force 2>&1 | tail -20

QREPORT=$ROOT/quality/checkm2_coassembly/quality_report.tsv
echo
if [ -s "$QREPORT" ]; then
  echo "=== Top 20 bins by Completeness ==="
  head -1 "$QREPORT"
  tail -n +2 "$QREPORT" | sort -t$'\t' -k2,2 -rn | head -20
  echo
  echo "=== Quality summary ==="
  awk -F'\t' 'NR>1 {c=$2; x=$3;
    if (c>=90 && x<=5) hi++;
    else if (c>=70 && x<=10) mq++;
    else if (c>=50 && x<=10) lq++;
    else frag++}
    END {print "HQ (C>=90, X<=5): " (hi+0); print "MQ (C>=70, X<=10): " (mq+0); print "Med (C>=50, X<=10): " (lq+0); print "Frag: " (frag+0)}' "$QREPORT"
else
  echo "CheckM2 report missing"; exit 2
fi

mkdir -p $ROOT/taxonomy/coassembly_gtdbtk_in
rm -f $ROOT/taxonomy/coassembly_gtdbtk_in/*.fa
awk -F'\t' 'NR>1 && $2>=70 && $3<=10 {print $1}' "$QREPORT" | while read bin; do
  src=$CLEAN_BINS/${bin}.fa
  [ -s "$src" ] && ln -sf "$src" $ROOT/taxonomy/coassembly_gtdbtk_in/${bin}.fa
done
NHQ=$(ls $ROOT/taxonomy/coassembly_gtdbtk_in/*.fa 2>/dev/null | wc -l)
echo
echo "bins for GTDB-Tk (>=70% C, <=10% X): $NHQ"

echo
echo "=== [$(date)] 5. GTDB-Tk classify_wf on $NHQ bins ==="
conda activate metagenomics
mkdir -p $ROOT/taxonomy/gtdbtk_coassembly
$GTDBTK classify_wf --genome_dir $ROOT/taxonomy/coassembly_gtdbtk_in -x fa \
  --out_dir $ROOT/taxonomy/gtdbtk_coassembly \
  --cpus $CORES --skip_ani_screen --pplacer_cpus $CORES 2>&1 | tail -25

echo
echo "=== CSP1-2 / Dadabacteria / Desulfobacterota_D hits ==="
grep -rhiE 'CSP1-2|Dadabacteria|Desulfobacterota_D|UBA2774' $ROOT/taxonomy/gtdbtk_coassembly/ 2>/dev/null | head || echo "  no matches"

echo "=== [$(date)] DONE ==="
