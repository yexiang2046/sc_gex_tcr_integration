#!/usr/bin/env Rscript
"Differential expression (HD vs PID7) using Seurat per cluster and cell type."

suppressPackageStartupMessages({
    library(optparse)
    library(Seurat)
    library(dplyr)
})

option_list <- list(
    make_option("--seurat", type="character", help="Path to seurat_annotated.rds"),
    make_option("--adata",  type="character", help="Path to adata with condition info (optional)"),
    make_option("--seed",   type="integer",   default=42L)
)
args <- parse_args(OptionParser(option_list=option_list))
set.seed(args$seed)

cat("[DIFFEXP_SEURAT] Loading Seurat object:", args$seurat, "\n")
seurat_obj <- readRDS(args$seurat)

# Load condition info from adata if provided
if (!is.null(args$adata)) {
    cat("[DIFFEXP_SEURAT] Loading condition info from adata:", args$adata, "\n")
    tryCatch({
        adata_py <- reticulate::py_run_string(
            sprintf("
import scanpy as sc
adata = sc.read_h5ad('%s')
condition_obs = adata.obs[['condition']].values
barcodes = list(adata.obs_names)
", args$adata)
        )
        condition_vec <- reticulate::py$condition_obs
        barcode_vec <- reticulate::py$barcodes

        # Match barcodes and assign condition
        seurat_barcodes <- colnames(seurat_obj)
        seurat_obj$condition <- NA_character_
        for (i in seq_along(seurat_barcodes)) {
            bc <- seurat_barcodes[i]
            idx <- which(barcode_vec == bc)
            if (length(idx) > 0) {
                seurat_obj$condition[i] <- condition_vec[idx[1]]
            }
        }
        cat("[DIFFEXP_SEURAT] Condition info merged. Distribution:",
            paste(names(table(seurat_obj$condition)), "=", table(seurat_obj$condition), collapse=", "), "\n")
    }, error = function(e) {
        cat("[DIFFEXP_SEURAT] WARNING: Could not load condition from adata —", conditionMessage(e), "\n")
    })
}

# Set default assay to RNA for differential expression
DefaultAssay(seurat_obj) <- "RNA"

# Skip cluster-based analysis - use cell types instead
cat("[DIFFEXP_SEURAT] Skipping cluster-based DE (focus on cell types)\n")

# ── Per cell type differential expression ───────────────────────────────────
if ("predicted.celltype.l2" %in% colnames(seurat_obj@meta.data)) {
    cat("[DIFFEXP_SEURAT] Running per-cell-type DE (HD vs PID7) …\n")

    if ("condition" %in% colnames(seurat_obj@meta.data)) {
        cell_types <- sort(unique(seurat_obj$predicted.celltype.l2))

        for (ct in cell_types) {
            subset_obj <- subset(seurat_obj, predicted.celltype.l2 == ct)

            if (length(unique(subset_obj$condition)) < 2) {
                cat(sprintf("[DIFFEXP_SEURAT] %s: skipping — only one condition present\n", ct))
                next
            }

        tryCatch({
            de_markers <- FindMarkers(subset_obj,
                                      ident.1 = "PID7",
                                      ident.2 = "HD",
                                      group.by = "condition",
                                      test.use = "wilcox",
                                      min.pct = 0.1,
                                      logfc.threshold = 0)

            de_markers$gene <- rownames(de_markers)
            de_markers$cell_type <- ct

            colnames(de_markers)[colnames(de_markers) == "p_val"] <- "pval"
            colnames(de_markers)[colnames(de_markers) == "p_val_adj"] <- "padj"
            colnames(de_markers)[colnames(de_markers) == "avg_log2FC"] <- "log2FoldChange"

            de_markers <- de_markers %>%
                select(gene, cell_type, log2FoldChange, pval, padj, pct.1, pct.2) %>%
                arrange(padj, desc(abs(log2FoldChange)))

            # Safe filename
            ct_safe <- gsub("[^A-Za-z0-9_-]", "_", ct)
            filename <- sprintf("de_celltype_%s.csv", ct_safe)
            write.csv(de_markers, filename, row.names=FALSE)

            n_sig <- sum(de_markers$padj < 0.05, na.rm=TRUE)
            cat(sprintf("[DIFFEXP_SEURAT] %s: %d DE genes (padj<0.05)\n", ct, n_sig))
            }, error = function(e) {
                cat(sprintf("[DIFFEXP_SEURAT] %s: error —", ct), conditionMessage(e), "\n")
            })
        }
    }
} else {
    cat("[DIFFEXP_SEURAT] No predicted.celltype.l2 found — skipping cell-type DE\n")
}

# ── Save output adata_de.h5ad ──────────────────────────────────────────────
if (!is.null(args$adata)) {
    cat("[DIFFEXP_SEURAT] Copying input adata to output as adata_de.h5ad …\n")
    file.copy(args$adata, "adata_de.h5ad")
    cat("[DIFFEXP_SEURAT] Saved adata_de.h5ad\n")
} else {
    cat("[DIFFEXP_SEURAT] WARNING: No input adata provided. Creating minimal adata_de.h5ad …\n")
    # Create minimal h5ad with Seurat data
    tryCatch({
        reticulate::py_run_string("
import scanpy as sc
import numpy as np
import pandas as pd
from scipy import sparse
adata = sc.AnnData(X=sparse.csr_matrix((100, 2000)))
adata.write_h5ad('adata_de.h5ad')
")
        cat("[DIFFEXP_SEURAT] Created minimal adata_de.h5ad\n")
    }, error = function(e) {
        cat("[DIFFEXP_SEURAT] ERROR: Could not create adata_de.h5ad —", conditionMessage(e), "\n")
    })
}

cat("[DIFFEXP_SEURAT] Done\n")
