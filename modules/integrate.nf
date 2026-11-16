process INTEGRATE {
    publishDir "${params.outdir}/integrate", mode: 'copy'
    container 'local/integrate:latest'

    input:
    path(adata_merged)

    output:
    path("adata_integrated.h5ad"), emit: adata
    path("*.png"),                  emit: figures

    script:
    """
    integrate.py \\
        --adata  "${adata_merged}" \\
        --n_hvgs ${params.n_hvgs} \\
        --seed   ${params.seed}
    """
}
