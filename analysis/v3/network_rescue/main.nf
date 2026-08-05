nextflow.enable.dsl = 2

params.project_root = params.project_root ?:
    file("${projectDir}/../../..").toAbsolutePath()
params.input_table = params.input_table ?:
    "analysis/v2/review/cache/genus_counts.tsv"
params.outdir = params.outdir ?:
    file("${params.project_root}/analysis/v3/network_rescue/results")

process NETWORK_RESCUE {
    tag "matched compositional networks"
    cpus 4
    memory "16 GB"
    time "4h"
    conda "${projectDir}/../../../workflow/environment.yml"
    publishDir params.outdir, mode: "copy", overwrite: true

    input:
    val project_root
    val input_table

    output:
    path "network_rescue/*"

    script:
    """
    python3 '${project_root}/analysis/v3/network_rescue/run_network_rescue.py' \
      --project-root '${project_root}' \
      --input-table '${input_table}' \
      --output-dir network_rescue
    """
}

workflow {
    NETWORK_RESCUE(
        Channel.value(params.project_root.toString()),
        Channel.value(params.input_table.toString())
    )
}
