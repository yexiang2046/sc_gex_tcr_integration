#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include { QC        } from './modules/qc'
include { NORMALIZE } from './modules/normalize'
include { MERGE     } from './modules/merge'
include { INTEGRATE } from './modules/integrate'
include { CLUSTER   } from './modules/cluster'
include { MARKERS   } from './modules/markers'
include { ANNOTATE  } from './modules/annotate'
include { DIFFEXP         } from './modules/diffexp'
include { SEURAT_CONVERT  } from './modules/seurat_convert'
include { SCENIC          } from './modules/scenic'
include { ANNOTATE_SEURAT } from './modules/annotate_seurat'

workflow {
    ch_samples = Channel
        .fromPath(params.input)
        .splitCsv(header: true)
        .map { row ->
            def meta = [
                sample     : row.Sample,
                patient_id : row['Patient ID'],
                description: row.Description,
                chemistry  : row.Chemistry
            ]
            tuple(meta, file("${projectDir}/${row.cellranger_outs}"))
        }

    QC(ch_samples)
    NORMALIZE(QC.out.adata)

    // Separate metas and paths before collecting — .collect() on a tuple channel
    // produces a flat list of tuples that cannot be destructured as tuple val, path.
    ch_norm_metas  = NORMALIZE.out.adata.map { meta, adata -> meta  }.collect()
    ch_norm_adatas = NORMALIZE.out.adata.map { meta, adata -> adata }.collect()
    MERGE(ch_norm_metas, ch_norm_adatas)

    INTEGRATE(MERGE.out.adata)
    CLUSTER(INTEGRATE.out.adata)
    MARKERS(CLUSTER.out.adata)
    SEURAT_CONVERT(CLUSTER.out.adata)
    ANNOTATE(MARKERS.out.adata)
    ANNOTATE_SEURAT(SEURAT_CONVERT.out.seurat_object)
    SCENIC(ANNOTATE.out.adata, ANNOTATE_SEURAT.out.mapquery_labels)
    DIFFEXP(ANNOTATE_SEURAT.out.seurat_object, ANNOTATE.out.adata)

    // Note: PATHWAY_ANALYSIS is run separately via Docker on published DIFFEXP results
    // because Nextflow channel emission of glob patterns has issues with optional outputs
}
