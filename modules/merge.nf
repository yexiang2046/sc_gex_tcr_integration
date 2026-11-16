process MERGE {
    publishDir "${params.outdir}/merge", mode: 'copy'
    container 'local/merge:latest'

    input:
    val(metas)
    path(adatas, stageAs: "inputs/adata_norm_?.h5ad")

    output:
    path("adata_merged.h5ad"), emit: adata
    path("*.png"),              emit: figures

    script:
    def metas_json = groovy.json.JsonOutput.toJson(metas)
    """
    merge.py \\
        --metas  '${metas_json}' \\
        --adatas ${adatas instanceof List ? adatas.join(' ') : adatas}
    """
}
