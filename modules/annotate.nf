process ANNOTATE {
    publishDir "${params.outdir}/annotate", mode: 'copy'
    container 'local/annotate:latest'

    input:
    path(adata_markers)

    output:
    path("adata_annotated.h5ad"),        emit: adata
    path("*.png"),                        emit: figures
    path("cluster_cell_type_mapping.csv"), emit: mapping

    script:
    """
    annotate.py \\
        --adata "${adata_markers}" \\
        --model "${params.celltypist_model}" \\
        --seed  ${params.seed}
    """
}
