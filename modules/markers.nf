process MARKERS {
    publishDir "${params.outdir}/markers", mode: 'copy'
    container 'local/markers:latest'

    input:
    path(adata_clustered)

    output:
    path("adata_markers.h5ad"),    emit: adata
    path("*.png"),                  emit: figures
    path("markers_cluster_*.csv"), emit: marker_tables

    script:
    """
    markers.py \\
        --adata         "${adata_clustered}" \\
        --n_top_markers ${params.n_top_markers} \\
        --seed          ${params.seed}
    """
}
