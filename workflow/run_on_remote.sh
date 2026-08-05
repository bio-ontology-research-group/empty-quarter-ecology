#!/usr/bin/env bash
set -euo pipefail

host=$(hostname -s)
case "$host" in
  cbontsr01|leechuck-office) ;;
  *)
    printf 'Refusing real workflow on %s; use ws or Ontolinator.\n' "$host" >&2
    exit 69
    ;;
esac

root=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
data_repo=${1:-"$root/../empty-quarter-data-paper"}
output_dir=${2:-"$root/results/remote-validation-$(date -u +%Y%m%dT%H%M%SZ)"}
data_repo=$(cd "$data_repo" && pwd)

if [[ -e "$output_dir" ]]; then
  printf 'refusing to overwrite existing output: %s\n' "$output_dir" >&2
  exit 73
fi

bash "$root/scripts/release/bootstrap_data_dependency.sh" "$data_repo"
python3 "$root/scripts/release/verify_repository.py" "$root"
EQ_DATA_REPO="$data_repo" python3 -m pytest -q "$root/tests"
python3 "$root/scripts/release/render_figures.py" "$root"

python3 "$data_repo/scripts/manuscript/test_manuscript_consistency.py"
python3 -m pytest -q "$data_repo/tests" "$data_repo/workflow/tests"

"$data_repo/workflow/bin/bootstrap_nextflow.sh" run \
  "$data_repo/workflow/main.nf" \
  -profile bare \
  --project_root "$data_repo" \
  --ecology_paper "$root/empty-quarter-amplicon" \
  --stage full \
  --run_advanced true \
  --run_kg true \
  --build_papers true \
  --pma_asv_table "$data_repo/metadata/relic-dna/PMA_ASV_table.tsv" \
  --coverm_dir "$data_repo/metadata/metagenome/coverm_profiles.tar.gz" \
  --eggnog_annotations "$data_repo/metadata/metagenome/eq.emapper.annotations.gz" \
  --measured_function_inputs \
    "$data_repo/metadata/metagenome/measured_function_inputs.tar.gz" \
  --taxonomy_source_taxonomy \
    "$data_repo/metadata/taxonomy/taxonomy-trips1-5.tsv" \
  --taxonomy_feature_table \
    "$data_repo/metadata/taxonomy/feature-table-trips1-5.tsv" \
  --taxonomy_canonical_mapping "$data_repo/ontology/mapped_taxonomy.tsv" \
  --taxonomy_ncbi_owl "$data_repo/data/ontologies/ncbitaxon.owl" \
  --taxonomy_sra_sheet \
    "$data_repo/metadata/sra-submissions/submission-sheet.tsv" \
  --outdir "$output_dir"

printf 'PASS: complete remote workflow written to %s\n' "$output_dir"
