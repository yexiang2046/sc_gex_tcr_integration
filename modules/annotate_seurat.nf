process ANNOTATE_SEURAT {
    publishDir "${params.outdir}/annotate_seurat", mode: 'copy'
    container 'xiang2019/seurat_convert:v1.0.0'

    input:
    path(seurat_rds)

    output:
    path("seurat_annotated.rds"),      emit: seurat_object
    path("mapquery_labels.csv"),       emit: mapquery_labels, optional: true
    path("*.png"),                     emit: figures

    script:
    """
    annotate_seurat.R \\
        --seurat    "${seurat_rds}" \\
        --reference /references/pbmc_multimodal.h5seurat \\
        --seed      ${params.seed}
    """
}
