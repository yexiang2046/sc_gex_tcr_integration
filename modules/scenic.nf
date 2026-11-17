process SCENIC {
    publishDir "${params.outdir}/scenic", mode: 'copy'
    container 'local/scenic:latest'

    // Bind the database directory into the container at its absolute host path.
    // file().toAbsolutePath() resolves relative params before the container starts.
    containerOptions "-v ${file(params.scenic_db_dir).toAbsolutePath()}:${file(params.scenic_db_dir).toAbsolutePath()}:ro"

    input:
    path(adata_annotated)
    path(mapquery_labels)

    output:
    path("adata_scenic.h5ad"),  emit: adata
    path("auc_matrix.csv"),     emit: auc_matrix
    path("regulons.pkl"),       emit: regulons
    path("adj.tsv"),            emit: adjacencies
    path("rss_scores.csv"),     emit: rss,     optional: true
    path("*.png"),              emit: figures, optional: true

    script:
    def abs_db_dir = file(params.scenic_db_dir).toAbsolutePath()
    def abs_tfs    = file(params.scenic_tfs).toAbsolutePath()
    def abs_motifs = file(params.scenic_motifs).toAbsolutePath()
    """
    # Run pySCENIC with MapQuery cell type labels (predicted.celltype.l2)
    # from ANNOTATE_SEURAT's mapquery_labels.csv — no R needed in this container.
    scenic.py \\
        --adata       "${adata_annotated}" \\
        --cell_types_csv "${mapquery_labels}" \\
        --db_dir      "${abs_db_dir}" \\
        --tfs         "${abs_tfs}" \\
        --motifs      "${abs_motifs}" \\
        --n_cells          ${params.scenic_n_cells} \\
        --num_workers      ${params.scenic_num_workers} \\
        --ctx_num_workers  ${params.scenic_ctx_num_workers} \\
        --n_hvgs           ${params.scenic_n_hvgs} \\
        --top_targets      ${params.scenic_top_targets} \\
        --module_chunksize ${params.scenic_module_chunksize} \\
        --seed             ${params.seed}
    """
}
