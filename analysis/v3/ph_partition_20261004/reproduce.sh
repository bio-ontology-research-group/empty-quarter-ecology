#!/usr/bin/env bash
set -euo pipefail
package_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
if [ "$#" -ne 1 ]; then echo 'Usage: bash reproduce.sh NEW_OUTPUT_DIRECTORY' >&2; exit 2; fi
if [ -e "$1" ]; then echo 'Output directory must not already exist.' >&2; exit 2; fi
ph_env_dir="$(mktemp -d -t eq-ph-env.XXXXXXXX)"
trap 'rm -rf -- "$ph_env_dir"' EXIT
uv venv --python 3.11.14 "$ph_env_dir/venv"
uv pip install --python "$ph_env_dir/venv/bin/python" -r "$package_dir/requirements.lock.txt"
export OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 OMP_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1
"$ph_env_dir/venv/bin/python" -m unittest discover -s "$package_dir/tests" -v
"$ph_env_dir/venv/bin/python" "$package_dir/analyze.py" --output "$1" --replicates 1999
