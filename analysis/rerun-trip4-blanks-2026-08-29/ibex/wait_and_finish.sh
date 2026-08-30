#!/bin/bash
# Wait for the Ibex denoising job to write DONE, then run finish.sh (rsync + per-trip screen + comparison).
DONE=/ibex/user/hohndor/eq-trip4-blanks-2026-08-29/DONE
for i in $(seq 1 60); do
  if ssh -o BatchMode=yes -o ConnectTimeout=20 ibex "test -f $DONE" 2>/dev/null; then
    echo "DONE seen at $(date)"; bash "$(dirname "$0")/finish.sh"; exit $?
  fi
  st=$(ssh -o BatchMode=yes -o ConnectTimeout=20 ibex "squeue -j 50996981 -h -o %T" 2>/dev/null)
  echo "$(date +%H:%M) job state: ${st:-not-in-queue}"
  if [ -z "$st" ] && ! ssh -o BatchMode=yes ibex "test -f $DONE" 2>/dev/null; then echo "job left the queue without DONE; check slurm-50996981.out"; exit 2; fi
  sleep 300
done
echo "timeout after 5h"; exit 3
