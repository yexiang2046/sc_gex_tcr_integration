process TCR_INTEGRATE {
    publishDir "${params.outdir}/tcr", mode: 'copy'
    container 'local/tcr_integrate:latest'

    // Bind the TCR root so absolute host paths from the samplesheet resolve in-container.
    // Mirrors the SCENIC db bind-mount pattern for large inputs kept outside the repo.
    containerOptions {
        def tcr_dir = file(params.tcr_dir)
        // Quote the bind so paths with spaces (e.g. data/TCR 12113) work in Docker.
        return tcr_dir.exists()
            ? "-v \"${tcr_dir.toAbsolutePath()}:${tcr_dir.toAbsolutePath()}:ro\""
            : ""
    }

    input:
    path(adata_annotated)
    path(mapquery_labels)
    val(tcr_map)

    output:
    path("adata_gex_tcr.h5ad"),       emit: adata
    path("mudata_gex_tcr.h5mu"),      emit: mudata
    path("tcr_summary.csv"),          emit: summary
    path("clonotypes_by_sample.csv"), emit: clonotypes
    path("tcr_input_status.csv"),     emit: input_status
    path("*.png"),                    emit: figures, optional: true

    script:
    def tcr_json = groovy.json.JsonOutput.toJson(tcr_map)
    """
    tcr_integrate.py \\
        --adata "${adata_annotated}" \\
        --mapquery_labels "${mapquery_labels}" \\
        --tcr_map '${tcr_json}' \\
        --seed ${params.seed}
    """
}
