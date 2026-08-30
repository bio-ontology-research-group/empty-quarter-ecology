#!/usr/bin/env bash
# Re-run of the extraction-blank contaminant screen with the Trip 4 blanks.
# Usage: bash run_rerun.sh                      (validation: Trip 5 screen only, Trip 4 not run)
#        TRIP4_BLANK_TABLE=/path/to/trip4_blanks.tsv bash run_rerun.sh   (full re-run, per-trip)
#        MODE=pooled TRIP4_BLANK_TABLE=... bash run_rerun.sh             (single pooled screen)
set -euo pipefail
HERE=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
ROOT=$(cd "$HERE/../.." && pwd)
PY=${PYTHON:-$ROOT/../data-paper/.venv/bin/python}
MODE=${MODE:-per-trip}
T4=${TRIP4_BLANK_TABLE:-}
cd "$HERE"
mkdir -p logs outputs
echo "start $(date -u +%FT%TZ) host=$(hostname) python=$PY mode=$MODE trip4_blank_table=${T4:-none}"
"$PY" --version
"$PY" build_trip4_batch_map.py > logs/build_trip4_batch_map.log
echo "batch map built: $(grep -c . inputs/extraction_batch_map_extended.tsv) rows"
if [[ -n "$T4" ]]; then
  "$PY" screen_extended.py --mode "$MODE" --trip4-blank-table "$T4" --require-trip4 --output-dir outputs 2>&1 | tee logs/screen_extended.log
else
  "$PY" screen_extended.py --mode "$MODE" --output-dir outputs 2>&1 | tee logs/screen_extended.log
fi
"$PY" compare.py > outputs/comparison.md 2> logs/compare.log
echo "comparison written to outputs/comparison.md"
echo "end $(date -u +%FT%TZ)"
