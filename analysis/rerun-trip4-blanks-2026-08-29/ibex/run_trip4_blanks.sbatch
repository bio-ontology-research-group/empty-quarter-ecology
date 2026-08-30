#!/bin/bash
#SBATCH --job-name=eq_trip4_blanks
#SBATCH --output=/ibex/user/hohndor/eq-trip4-blanks-2026-08-29/slurm-%j.out
#SBATCH --error=/ibex/user/hohndor/eq-trip4-blanks-2026-08-29/slurm-%j.out
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=24:00:00
# Head job for nf-core/ampliseq 2.14.0 on the six Trip 4 extraction blanks.
# Same settings as the canonical Trips 1-4 run (novaseq_14_07_25/raw_reads/ampliseq,
# Aug 2025) and the Trip 5 run (Trip5/analysis/ampliseq/results_trip5_fixed, Mar 2026):
# -profile kaust (slurm executor + singularity), 341F/806R primers, --illumina_novaseq,
# SILVA 138.2, defaults otherwise. truncLen 240/238 is what those runs derived
# automatically (dada2/args/filterAndTrim.args.txt); it is fixed here so six
# low-depth blanks do not re-derive a different value.  --skip_qiime is the only
# other deviation: QIIME 2 steps do not touch dada2/ASV_table.tsv.
set -euo pipefail
D=/ibex/user/hohndor/eq-trip4-blanks-2026-08-29
cd "$D"
module load nextflow/24.10.6
module load singularity/3.9.7
export NXF_SINGULARITY_CACHEDIR=/ibex/user/hohndor/.singularity/nf_images
echo "start $(date -u +%FT%TZ) host=$(hostname) jobid=${SLURM_JOB_ID:-none}"
nextflow -version | head -5
nextflow run nf-core/ampliseq -r 2.14.0 \
  -profile kaust \
  -c "$D/trip4_blanks.config" \
  -params-file "$D/params.json" \
  -w "$D/work" \
  -resume
echo "nextflow finished $(date -u +%FT%TZ)"
python3 "$D/make_md5_table.py" "$D/results/dada2/ASV_table.tsv" "$D/results/dada2/ASV_seqs.fasta" "$D/trip4_blanks_md5.tsv"
touch "$D/DONE"
echo "end $(date -u +%FT%TZ)"
