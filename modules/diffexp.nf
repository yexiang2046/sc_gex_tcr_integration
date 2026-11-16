process DIFFEXP {
    publishDir "${params.outdir}/diffexp", mode: 'copy'
    container 'xiang2019/seurat_convert:v1.0.0'

    input:
    path(seurat_rds)
    path(adata_with_condition)

    output:
    path("adata_de.h5ad"),              emit: adata
    path("de_celltype_*.csv"),          emit: de_celltype, optional: true

    script:
    """
    diffexp_seurat.R \\
        --seurat "${seurat_rds}" \\
        --adata  "${adata_with_condition}" \\
        --seed   ${params.seed}
    """
}
