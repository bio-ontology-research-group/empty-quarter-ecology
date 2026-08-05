#!/usr/bin/env bash
set -euo pipefail

root=$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)
data_repo=${1:-"$root/../empty-quarter-data-paper"}
data_repo=$(cd "$data_repo" && pwd)
lock="$root/DATA_REPOSITORY.lock"

lock_value() {
  awk -F '\t' -v key="$1" '$1 == key { print $2 }' "$lock"
}

require_digest() {
  local field=$1
  local relative=$2
  local expected
  local observed
  expected=$(lock_value "$field")
  observed=$(sha256sum "$data_repo/$relative" | cut -d ' ' -f 1)
  if [[ -z "$expected" || "$observed" != "$expected" ]]; then
    printf 'data dependency digest mismatch: %s\n' "$relative" >&2
    printf 'expected %s; observed %s\n' "$expected" "$observed" >&2
    exit 65
  fi
}

expected_commit=$(lock_value commit)
observed_commit=$(git -C "$data_repo" rev-parse HEAD)
if [[ "$observed_commit" != "$expected_commit" ]]; then
  printf 'data dependency must be at %s; observed %s\n' \
    "$expected_commit" "$observed_commit" >&2
  exit 65
fi

require_digest file_manifest_sha256 FILE_MANIFEST.tsv
require_digest bulk_manifest_sha256 BULK_ARTIFACTS.tsv
require_digest release_manifest_sha256 PRE_RELEASE_MANIFEST.tsv
require_digest environment_lock_sha256 environment/requirements.lock.txt
require_digest workflow_sha256 workflow/main.nf

bash "$data_repo/scripts/release/bootstrap_package_layout.sh" "$data_repo"
python3 "$data_repo/scripts/release/verify_repository.py" "$data_repo"

link_path() {
  local destination=$1
  local target_abs=$2
  mkdir -p "$(dirname "$destination")"
  local parent_abs
  local target
  parent_abs=$(cd "$(dirname "$destination")" && pwd -P)
  target=$(realpath --relative-to="$parent_abs" "$target_abs")
  if [[ -L "$destination" ]]; then
    if [[ $(readlink "$destination") == "$target" ]]; then
      return
    fi
    printf 'existing link has an unexpected target: %s\n' "$destination" >&2
    exit 73
  fi
  if [[ -e "$destination" ]]; then
    printf 'refusing to replace existing path: %s\n' "$destination" >&2
    exit 73
  fi
  ln -s "$target" "$destination"
}

link_path "$root/data" "$data_repo/data"
link_path "$root/data-paper" "$data_repo/paper"
link_path "$root/relic-dna/ASV_table.tsv" \
  "$data_repo/metadata/relic-dna/PMA_ASV_table.tsv"
link_path "$root/analysis/xrf_audit" "$data_repo/evidence/xrf_audit"

printf 'PASS: ecology dependency verified at %s\n' "$observed_commit"
