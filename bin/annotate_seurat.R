#!/usr/bin/env Rscript
"Annotate Seurat object: Seurat MapQuery label transfer against pbmc_multimodal reference (l2 by default)."

suppressPackageStartupMessages({
    library(optparse)
    library(Seurat)
    library(SeuratDisk)
    library(ggplot2)
    library(patchwork)
})

option_list <- list(
    make_option("--seurat",    type="character", help="Path to adata_clustered.rds (from SEURAT_CONVERT)"),
    make_option("--reference", type="character", default="/references/pbmc_multimodal.h5seurat"),
    make_option("--seed",      type="integer",   default=42L)
)
args <- parse_args(OptionParser(option_list=option_list))
set.seed(args$seed)

# ── Load Seurat object ────────────────────────────────────────────────────────
cat("[ANNOTATE_SEURAT] Loading Seurat object:", args$seurat, "\n")
seurat_obj <- readRDS(args$seurat)
cat(sprintf("[ANNOTATE_SEURAT] %d cells × %d genes\n", ncol(seurat_obj), nrow(seurat_obj)))

seurat_obj <- NormalizeData(seurat_obj, normalization.method="LogNormalize",
                            scale.factor=1e4, verbose=FALSE)
gc()

# ── Seurat label transfer against pbmc_multimodal reference ──────────────────
cat("[ANNOTATE_SEURAT] Loading reference (counts slot only, will free immediately) …\n")
reference <- LoadH5Seurat(args$reference,
                          assays     = list(SCT = "counts"),
                          reductions = c("spca", "wnn.umap"),
                          graphs     = NULL,
                          neighbors  = NULL,
                          images     = NULL)

# Replace SCT with an all-zero RNA assay carrying the SPCA gene names so that
# FindTransferAnchors can intersect reference genes with the query.
# Expression values are never read when reference.reduction="spca".
spca_genes <- rownames(reference@reductions[["spca"]]@feature.loadings)
cat(sprintf("[ANNOTATE_SEURAT] SPCA has %d genes; query has %d genes\n",
            length(spca_genes), nrow(seurat_obj)))
dummy <- Matrix::sparseMatrix(i=integer(0), j=integer(0), x=numeric(0),
                              dims=c(length(spca_genes), ncol(reference)),
                              dimnames=list(spca_genes, colnames(reference)))
reference[["RNA"]] <- CreateAssayObject(counts=dummy)
DefaultAssay(reference) <- "RNA"
reference[["SCT"]] <- NULL
gc()

# tryCatch returns the annotated Seurat object on success, NULL on failure.
# Assigning to seurat_obj happens outside tryCatch with plain <- to avoid the
# <<- compound-assignment scoping issue inside tryCatch bodies.
annotated <- tryCatch({
    features_use <- intersect(rownames(reference), rownames(seurat_obj))
    cat(sprintf("[ANNOTATE_SEURAT] FindTransferAnchors on %d cells, %d common features …\n",
                ncol(seurat_obj), length(features_use)))
    if (length(features_use) == 0L)
        stop("No common features between reference SPCA genes and query — gene name mismatch?")

    anchors <- FindTransferAnchors(
        reference            = reference,
        query                = seurat_obj,
        normalization.method = "LogNormalize",
        reference.reduction  = "spca",
        features             = features_use,
        dims                 = 1:50
    )
    cat("[ANNOTATE_SEURAT] TransferData …\n")
    result <- TransferData(
        anchorset = anchors,
        refdata   = list(
            celltype.l1 = "celltype.l1",
            celltype.l2 = "celltype.l2"
        ),
        reference = reference,
        query     = seurat_obj
    )
    rm(anchors, reference); gc()
    result
}, error = function(e) {
    cat("[ANNOTATE_SEURAT] WARNING: label transfer failed —", conditionMessage(e),
        "\nSaving Seurat object without MapQuery labels.\n")
    NULL
})

mapquery_ok <- !is.null(annotated)
if (mapquery_ok) {
    seurat_obj <- annotated
    cat("[ANNOTATE_SEURAT] pbmc_multimodal l1:\n")
    print(sort(table(seurat_obj$predicted.celltype.l1), decreasing=TRUE))
    cat("[ANNOTATE_SEURAT] pbmc_multimodal l2:\n")
    print(sort(table(seurat_obj$predicted.celltype.l2), decreasing=TRUE))
}
rm(annotated); gc()

# ── Plots ─────────────────────────────────────────────────────────────────────
umap_name <- if ("umap" %in% names(seurat_obj@reductions)) "umap" else
             names(seurat_obj@reductions)[1]

plot_dimplot <- function(obj, group_by, file, title, label=TRUE) {
    tryCatch({
        p <- DimPlot(obj, reduction=umap_name, group.by=group_by,
                     label=label, label.size=3, pt.size=0.1, repel=TRUE) +
             ggtitle(title) +
             theme(legend.text=element_text(size=7),
                   legend.key.size=unit(0.4, "cm"))
        ggsave(file, p, width=10, height=8, dpi=150)
        cat("[ANNOTATE_SEURAT] Saved", file, "\n")
    }, error=function(e)
        cat("[ANNOTATE_SEURAT] WARNING: could not save", file, ":", conditionMessage(e), "\n"))
}

if (mapquery_ok) {
    plot_dimplot(seurat_obj, "predicted.celltype.l2", "umap_mapquery_l2.png",
                 "MapQuery l2 (pbmc_multimodal)")
    plot_dimplot(seurat_obj, "predicted.celltype.l1", "umap_mapquery_l1.png",
                 "MapQuery l1 (pbmc_multimodal)", label=FALSE)
    tryCatch({
        p <- FeaturePlot(seurat_obj, features="predicted.celltype.l2.score",
                         reduction=umap_name, pt.size=0.1) +
             scale_color_viridis_c() + ggtitle("MapQuery l2 mapping score")
        ggsave("umap_mapquery_score.png", p, width=8, height=7, dpi=150)
    }, error=function(e) NULL)
}

# ── Save Seurat ────────────────────────────────────────────────────────────────
cat("[ANNOTATE_SEURAT] Saving seurat_annotated.rds …\n")
saveRDS(seurat_obj, "seurat_annotated.rds")

# ── Export MapQuery labels for DIFFEXP ──────────────────────────────────────────
if (mapquery_ok) {
    ct_export <- data.frame(
        barcode = colnames(seurat_obj),
        predicted.celltype.l2 = seurat_obj$predicted.celltype.l2,
        predicted.celltype.l2.score = seurat_obj$predicted.celltype.l2.score,
        stringsAsFactors = FALSE
    )
    write.csv(ct_export, "mapquery_labels.csv", row.names=FALSE)
    cat("[ANNOTATE_SEURAT] Exported MapQuery cell types to mapquery_labels.csv\n")
}

cat(sprintf("[ANNOTATE_SEURAT] Done. %d cells. Meta.data columns: %s\n",
            ncol(seurat_obj),
            paste(colnames(seurat_obj@meta.data), collapse=", ")))
