process CLUSTER {
    publishDir "${params.outdir}/cluster", mode: 'copy'
    container 'local/cluster:latest'

    input:
    path(adata_integrated)

    output:
    path("adata_clustered.h5ad"), emit: adata
    path("*.png"),                 emit: figures

    script:
    """
    cluster.py \\
        --adata             "${adata_integrated}" \\
        --leiden_resolution ${params.leiden_resolution} \\
        --seed              ${params.seed}
    """
}
