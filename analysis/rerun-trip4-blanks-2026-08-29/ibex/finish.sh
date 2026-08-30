#!/usr/bin/env bash
# Finish the Trip 4 blank re-run once Ibex job $(cat job.id) has completed
# (needs the KAUST network/VPN for ssh to ibex):
#   1. rsync the blank ASV table (md5 ids), sequences, args and pipeline_info from Ibex
#   2. run the per-trip contaminant screen with the Trip 4 blanks (run_rerun.sh)
#   3. compare.py is run by run_rerun.sh -> outputs/comparison.md
# Usage: bash ibex/finish.sh            (from anywhere)
#        MODE=pooled bash ibex/finish.sh  (single pooled screen instead of per-trip)
set -euo pipefail
HERE=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
TOP=$(cd "$HERE/.." && pwd)
REMOTE=/ibex/user/hohndor/eq-trip4-blanks-2026-08-29
JOB=$(cat "$HERE/job.id")
mkdir -p "$HERE/results"
state=$(ssh -o BatchMode=yes ibex "sacct -j $JOB -X -n -o State 2>/dev/null | head -1 | tr -d ' '; test -e $REMOTE/DONE && echo DONE" | tr '\n' ' ')
echo "job $JOB: $state"
if ! grep -q DONE <<<"$state"; then
  echo "Ibex job has not written $REMOTE/DONE yet (state: $state). Check: ssh ibex 'tail -30 $REMOTE/slurm-$JOB.out'"
  exit 1
fi
rsync -a --info=stats1 \
  "ibex:$REMOTE/trip4_blanks_md5.tsv" \
  "ibex:$REMOTE/slurm-$JOB.out" \
  "ibex:$REMOTE/results/dada2/ASV_table.tsv" \
  "ibex:$REMOTE/results/dada2/ASV_seqs.fasta" \
  "ibex:$REMOTE/results/dada2/ASV_tax.silva_138.2.tsv" \
  "ibex:$REMOTE/results/dada2/DADA2_stats.tsv" \
  "$HERE/results/"
rsync -a "ibex:$REMOTE/results/dada2/args" "ibex:$REMOTE/results/pipeline_info" "$HERE/results/"
# Rebuild the md5 table locally as a check (must match what Ibex wrote).
python3 "$HERE/make_md5_table.py" "$HERE/results/ASV_table.tsv" "$HERE/results/ASV_seqs.fasta" "$HERE/results/trip4_blanks_md5.local.tsv"
cmp "$HERE/results/trip4_blanks_md5.tsv" "$HERE/results/trip4_blanks_md5.local.tsv" && echo "md5 table verified"
# Keep the validation-run outputs (Trip 5 only) before they are overwritten.
if [[ -d "$TOP/outputs" && ! -d "$TOP/outputs_validation_trip5_only" ]]; then
  cp -a "$TOP/outputs" "$TOP/outputs_validation_trip5_only"
fi
cd "$TOP"
MODE=${MODE:-per-trip} TRIP4_BLANK_TABLE="$HERE/results/trip4_blanks_md5.tsv" bash run_rerun.sh 2>&1 | tee "$HERE/finish.log"
echo "done: see $TOP/outputs/comparison.md"
