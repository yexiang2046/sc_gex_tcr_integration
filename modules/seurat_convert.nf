process SEURAT_CONVERT {
    publishDir "${params.outdir}/seurat_convert", mode: 'copy'
    container 'xiang2019/seurat_convert:v1.0.0'

    input:
    path(adata_clustered)

    output:
    path("adata_clustered.rds"), emit: seurat_object

    script:
    def all_genes_flag = params.seurat_all_genes ? "--all_genes" : ""
    """
    seurat_convert.R \\
        --adata     "${adata_clustered}" \\
        --seed      ${params.seed} \\
        --max_cells ${params.seurat_max_cells} \\
        ${all_genes_flag}
    """
}
