process NORMALIZE {
    tag "${meta.sample}"
    publishDir "${params.outdir}/normalize/${meta.sample}", mode: 'copy'
    container 'local/normalize:latest'

    input:
    tuple val(meta), path(adata_qc)

    output:
    tuple val(meta), path("adata_norm.h5ad"), emit: adata
    path("*.png"),                             emit: figures

    script:
    """
    normalize.py \\
        --sample  "${meta.sample}" \\
        --adata   "${adata_qc}" \\
        --n_hvgs  ${params.n_hvgs} \\
        --seed    ${params.seed}
    """
}
