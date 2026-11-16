#!/usr/bin/env Rscript
"Compute AP-1 pathway module scores (MSigDB) on the annotated Seurat object.

Strategy
--------
seurat_annotated.rds  → cell-type labels + UMAP reduction (dummy 1-gene assay)
adata_clustered.rds   → real gene-expression count matrix (no annotations)

Steps:
  1. Load both objects; transfer annotation columns + UMAP from seurat_annotated
     into the full-matrix object.
  2. NormalizeData (LogNormalize) — required before AddModuleScore.
  3. Fetch AP-1 gene sets from MSigDB via msigdbr:
       • C3:TFT  — AP-1 motif-based target sets  (AP1_*)
       • C3:TFT:GTRD — GTRD ChIP-seq targets for AP-1 family TFs
         (JUNB, JUN, JUND, FOS, FOSL1, FOSL2, FOSB, BATF, BATF3, ATF3)
  4. AddModuleScore for each gene set.
  5. Save enriched Seurat as RDS; output UMAP feature plots + violin plots."

suppressPackageStartupMessages({
    library(optparse)
    library(Seurat)
    library(dplyr)
    library(ggplot2)
    library(patchwork)
})

# ── Install msigdbr if absent ─────────────────────────────────────────────────
if (!requireNamespace("msigdbr", quietly = TRUE)) {
    message("[AP1_SCORE] Installing msigdbr …")
    install.packages("msigdbr", repos = "https://cloud.r-project.org", quiet = TRUE)
}
library(msigdbr)

# ── CLI arguments ─────────────────────────────────────────────────────────────
option_list <- list(
    make_option("--seurat_full",
                type    = "character",
                default = "results/seurat_convert/adata_clustered.rds",
                help    = "Full-matrix Seurat RDS (from SEURAT_CONVERT)"),
    make_option("--seurat_annotated",
                type    = "character",
                default = "results/annotate_seurat/seurat_annotated.rds",
                help    = "Annotated Seurat RDS (from ANNOTATE_SEURAT, dummy assay + metadata)"),
    make_option("--outdir",
                type    = "character",
                default = "results/ap1_module_score",
                help    = "Output directory"),
    make_option("--group_by",
                type    = "character",
                default = "predicted.celltype.l2",
                help    = "obs column for violin grouping (cell type)"),
    make_option("--seed",
                type    = "integer",
                default = 42L)
)
args <- parse_args(OptionParser(option_list = option_list))
set.seed(args$seed)
dir.create(args$outdir, showWarnings = FALSE, recursive = TRUE)

# ── 1. Load objects ───────────────────────────────────────────────────────────
cat("[AP1_SCORE] Loading full-matrix Seurat:", args$seurat_full, "\n")
seurat_obj <- readRDS(args$seurat_full)
cat(sprintf("[AP1_SCORE] %d genes × %d cells\n", nrow(seurat_obj), ncol(seurat_obj)))

cat("[AP1_SCORE] Loading annotated Seurat:", args$seurat_annotated, "\n")
annotated  <- readRDS(args$seurat_annotated)

# ── 2. Transfer annotations + UMAP from annotated → full-matrix object ────────
anno_cols <- c(
    "cell_type_celltypist", "cell_type_confidence_celltypist",
    "predicted.celltype.l1", "predicted.celltype.l2", "predicted.celltype.l2.score"
)
anno_cols <- anno_cols[anno_cols %in% colnames(annotated@meta.data)]

shared_cells <- intersect(colnames(seurat_obj), colnames(annotated))
cat(sprintf("[AP1_SCORE] %d / %d cells matched between objects\n",
            length(shared_cells), ncol(seurat_obj)))

for (col in anno_cols) {
    vals <- setNames(annotated@meta.data[shared_cells, col], shared_cells)
    seurat_obj <- AddMetaData(seurat_obj, metadata = vals, col.name = col)
}

# Transfer UMAP if present in annotated and not already in full object
if ("umap" %in% names(annotated@reductions) &&
    !"umap" %in% names(seurat_obj@reductions)) {
    umap_emb            <- annotated@reductions[["umap"]]@cell.embeddings
    umap_shared         <- umap_emb[shared_cells, , drop = FALSE]
    seurat_obj[["umap"]] <- CreateDimReducObject(
        embeddings = umap_shared[colnames(seurat_obj), , drop = FALSE],
        key        = "UMAP_",
        assay      = "RNA"
    )
    cat("[AP1_SCORE] UMAP transferred from annotated object\n")
}
rm(annotated); gc()

# ── 3. Normalize ──────────────────────────────────────────────────────────────
cat("[AP1_SCORE] NormalizeData …\n")
seurat_obj <- NormalizeData(seurat_obj, normalization.method = "LogNormalize",
                            scale.factor = 1e4, verbose = FALSE)

# ── 4. Fetch AP-1 gene sets from MSigDB ───────────────────────────────────────
cat("[AP1_SCORE] Fetching AP-1 gene sets from MSigDB …\n")

# 4a. Classic AP-1 motif targets (C3: Transcription Factor Targets)
#     Names match: AP1_Q4, AP1_Q6, AP1_Q4_01, AP1_Q6_01
tft_all  <- msigdbr(species = "Homo sapiens", category = "C3")
ap1_motif_sets <- tft_all %>%
    filter(grepl("^AP1_", gs_name)) %>%
    select(gs_name, gene_symbol) %>%
    distinct()

if (nrow(ap1_motif_sets) == 0) {
    cat("[AP1_SCORE] WARNING: no AP1_ motif sets found in C3 — skipping motif module\n")
    ap1_motif_genes <- NULL
} else {
    ap1_motif_genes <- unique(ap1_motif_sets$gene_symbol)
    cat(sprintf("[AP1_SCORE] AP-1 motif union (%d unique genes from: %s)\n",
                length(ap1_motif_genes),
                paste(unique(ap1_motif_sets$gs_name), collapse = ", ")))
}

# 4b. GTRD ChIP-seq targets for AP-1 family TFs
AP1_FAMILY <- c("JUNB", "JUN", "JUND", "FOS", "FOSL1", "FOSL2", "FOSB",
                "BATF", "BATF3", "ATF3")

gtrd <- msigdbr(species = "Homo sapiens", category = "C3", subcategory = "TFT:GTRD")

gtrd_ap1 <- gtrd %>%
    filter(tf_gene_symbol %in% AP1_FAMILY | grepl(
        paste(AP1_FAMILY, collapse = "|"), gs_name, ignore.case = TRUE
    )) %>%
    select(gs_name, gene_symbol, tf_gene_symbol) %>%
    distinct()

if (nrow(gtrd_ap1) == 0) {
    # Fallback: match by set name pattern
    gtrd_ap1 <- gtrd %>%
        filter(grepl(paste(AP1_FAMILY, collapse = "|"), gs_name)) %>%
        select(gs_name, gene_symbol) %>%
        distinct()
}

gtrd_sets <- split(gtrd_ap1$gene_symbol, gtrd_ap1$gs_name)
cat(sprintf("[AP1_SCORE] GTRD AP-1 family: %d gene sets\n", length(gtrd_sets)))
for (nm in names(gtrd_sets))
    cat(sprintf("             %s: %d genes\n", nm, length(gtrd_sets[[nm]])))

# ── 5. Build final module list & filter to genes present in object ────────────
gene_universe <- rownames(seurat_obj)

build_module <- function(genes, label) {
    g <- intersect(unique(genes), gene_universe)
    if (length(g) < 5L) {
        cat(sprintf("[AP1_SCORE] WARNING: '%s' has only %d genes in matrix — skipping\n",
                    label, length(g)))
        return(NULL)
    }
    cat(sprintf("[AP1_SCORE] Module '%s': %d / %d genes in matrix\n",
                label, length(g), length(unique(genes))))
    list(g)   # AddModuleScore expects a list of character vectors
}

modules <- list()

if (!is.null(ap1_motif_genes)) {
    m <- build_module(ap1_motif_genes, "AP1_motif_union")
    if (!is.null(m)) modules[["AP1_motif_union"]] <- m[[1]]
}

for (nm in names(gtrd_sets)) {
    m <- build_module(gtrd_sets[[nm]], nm)
    if (!is.null(m)) modules[[nm]] <- m[[1]]
}

if (length(modules) == 0L)
    stop("[AP1_SCORE] No valid AP-1 modules found — check MSigDB version or gene names.")

cat(sprintf("[AP1_SCORE] Running AddModuleScore for %d modules …\n", length(modules)))

# ── 6. AddModuleScore ─────────────────────────────────────────────────────────
# Score each module individually so the column name is unambiguous
for (mod_name in names(modules)) {
    col_name <- paste0(mod_name, "_score")
    seurat_obj <- AddModuleScore(
        object   = seurat_obj,
        features = list(modules[[mod_name]]),
        name     = col_name,
        seed     = args$seed
    )
    # AddModuleScore appends "1" to the name
    actual_col <- paste0(col_name, "1")
    colnames(seurat_obj@meta.data)[
        colnames(seurat_obj@meta.data) == actual_col
    ] <- col_name
    cat(sprintf("[AP1_SCORE]   added column: %s\n", col_name))
}

score_cols <- paste0(names(modules), "_score")

# ── 7. Save enriched object ───────────────────────────────────────────────────
out_rds <- file.path(args$outdir, "seurat_ap1_scored.rds")
cat("[AP1_SCORE] Saving scored Seurat …\n")
gc()
saveRDS(seurat_obj, out_rds, compress = FALSE)
cat("[AP1_SCORE] Saved:", out_rds, "\n")

# ── 8. Plots ──────────────────────────────────────────────────────────────────
umap_name <- if ("umap" %in% names(seurat_obj@reductions)) "umap" else
             names(seurat_obj@reductions)[1]

# 8a. UMAP feature plots — one per module
for (sc_col in score_cols) {
    tryCatch({
        p <- FeaturePlot(seurat_obj, features = sc_col,
                         reduction = umap_name, pt.size = 0.1,
                         order = TRUE) +
             scale_color_gradientn(colours = c("lightgrey", "#2166AC", "#D6604D", "#B2182B")) +
             ggtitle(gsub("_score$", "", sc_col)) +
             theme(plot.title = element_text(size = 10))
        fname <- file.path(args$outdir, paste0("umap_", sc_col, ".png"))
        ggsave(fname, p, width = 7, height = 6, dpi = 150)
        cat("[AP1_SCORE] Saved", fname, "\n")
    }, error = function(e)
        cat("[AP1_SCORE] WARNING: could not save UMAP for", sc_col, ":", conditionMessage(e), "\n"))
}

# 8b. Violin plots grouped by cell type — one PNG with all modules
group_col <- args$group_by
if (!group_col %in% colnames(seurat_obj@meta.data)) {
    # Fallback order
    fallbacks <- c("cell_type_celltypist", "predicted.celltype.l1", "leiden")
    group_col <- fallbacks[fallbacks %in% colnames(seurat_obj@meta.data)][1]
    cat(sprintf("[AP1_SCORE] '%s' not found — using '%s' for violin grouping\n",
                args$group_by, group_col))
}

if (!is.na(group_col)) {
    vln_plots <- lapply(score_cols, function(sc_col) {
        tryCatch(
            VlnPlot(seurat_obj, features = sc_col, group.by = group_col,
                    pt.size = 0, sort = "increasing") +
                NoLegend() +
                ggtitle(gsub("_score$", "", sc_col)) +
                theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
                      plot.title  = element_text(size = 9)),
            error = function(e) NULL
        )
    })
    vln_plots <- Filter(Negate(is.null), vln_plots)
    if (length(vln_plots) > 0L) {
        ncols <- min(3L, length(vln_plots))
        p_vln <- wrap_plots(vln_plots, ncol = ncols)
        fname_vln <- file.path(args$outdir, "violin_ap1_scores.png")
        ggsave(fname_vln, p_vln,
               width  = 7 * ncols,
               height = 5 * ceiling(length(vln_plots) / ncols),
               dpi    = 150, limitsize = FALSE)
        cat("[AP1_SCORE] Saved", fname_vln, "\n")
    }
}

# 8c. Split violin by condition (HD vs PID7) for the motif-union score
if ("AP1_motif_union_score" %in% score_cols &&
    "condition" %in% colnames(seurat_obj@meta.data) &&
    !is.na(group_col)) {
    tryCatch({
        p_split <- VlnPlot(seurat_obj,
                           features  = "AP1_motif_union_score",
                           group.by  = group_col,
                           split.by  = "condition",
                           pt.size   = 0,
                           sort      = FALSE) +
                   theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7)) +
                   ggtitle("AP-1 motif score — HD vs PID7 by cell type")
        ggsave(file.path(args$outdir, "violin_ap1_motif_split_condition.png"),
               p_split, width = 16, height = 6, dpi = 150)
        cat("[AP1_SCORE] Saved violin_ap1_motif_split_condition.png\n")
    }, error = function(e)
        cat("[AP1_SCORE] WARNING: split violin failed —", conditionMessage(e), "\n"))
}

cat(sprintf("[AP1_SCORE] Done. Outputs in %s/\n", args$outdir))
