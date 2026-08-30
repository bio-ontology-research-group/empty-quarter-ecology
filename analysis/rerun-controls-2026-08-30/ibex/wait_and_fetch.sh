#!/bin/bash
# Poll Ibex until eq-controls-2026-08-30 writes DONE, then rsync the outputs here.
HERE=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
RD=/ibex/user/hohndor/eq-controls-2026-08-30
JOB=$(cat "$HERE/job.id")
for i in $(seq 1 72); do
  if ssh -o BatchMode=yes -o ConnectTimeout=20 ibex "test -f $RD/DONE" 2>/dev/null; then
    echo "DONE seen at $(date)"
    mkdir -p "$HERE/results"
    rsync -a "ibex:$RD/controls_md5.tsv" "ibex:$RD/slurm-$JOB.out" "ibex:$RD/results/dada2/ASV_table.tsv" "ibex:$RD/results/dada2/ASV_seqs.fasta" "ibex:$RD/results/dada2/ASV_tax.silva_138.2.tsv" "ibex:$RD/results/dada2/DADA2_stats.tsv" "$HERE/results/"
    rsync -a "ibex:$RD/results/dada2/args" "ibex:$RD/results/pipeline_info" "$HERE/results/"
    python3 "$HERE/make_md5_table.py" "$HERE/results/ASV_table.tsv" "$HERE/results/ASV_seqs.fasta" "$HERE/results/controls_md5.local.tsv"
    cmp "$HERE/results/controls_md5.tsv" "$HERE/results/controls_md5.local.tsv" && echo "md5 table verified"
    exit 0
  fi
  st=$(ssh -o BatchMode=yes -o ConnectTimeout=20 ibex "squeue -j $JOB -h -o %T" 2>/dev/null | tail -1)
  echo "$(date +%H:%M) job $JOB state: ${st:-not-in-queue}"
  if [ -z "$st" ]; then
    sleep 60
    if ! ssh -o BatchMode=yes ibex "test -f $RD/DONE" 2>/dev/null; then echo "job left the queue without DONE; check slurm-$JOB.out"; exit 2; fi
  fi
  sleep 120
done
echo "timeout"; exit 3
