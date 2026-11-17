process PATHWAY_ANALYSIS {
    publishDir "${params.outdir}/pathway_analysis", mode: 'copy'
    container 'local/pathway_analysis:latest'

    input:
    path(de_csvs)

    output:
    path("gsea_results.csv"),              emit: results, optional: true
    path("gsea_*.png"),                    emit: figures, optional: true

    script:
    """
    pathway_analysis.R \\
        --de_csvs . \\
        --outdir . \\
        --categories "${params.pathway_categories}" \\
        --species "${params.pathway_species}" \\
        --pval_thresh ${params.pathway_pval_thresh} \\
        --top_pathways ${params.pathway_top_pathways} \\
        --seed ${params.seed}
    """
}
