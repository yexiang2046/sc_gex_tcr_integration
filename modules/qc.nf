process QC {
    tag "${meta.sample}"
    publishDir "${params.outdir}/qc/${meta.sample}", mode: 'copy'
    container 'local/qc:latest'

    input:
    tuple val(meta), path(cellranger_dir)

    output:
    tuple val(meta), path("adata_qc.h5ad"), emit: adata
    path("*.png"),                           emit: figures
    path("qc_metrics.csv"),                  emit: metrics

    script:
    """
    qc.py \\
        --sample      "${meta.sample}" \\
        --cellranger_dir "${cellranger_dir}" \\
        --chemistry   "${meta.chemistry}" \\
        --min_genes   ${params.min_genes} \\
        --max_genes   ${params.max_genes} \\
        --mt_threshold ${params.mt_threshold} \\
        --min_cells   ${params.min_cells} \\
        --seed        ${params.seed}
    """
}
