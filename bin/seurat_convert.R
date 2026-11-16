#!/usr/bin/env Rscript
"Convert adata_clustered.h5ad to a Seurat RDS object.
 Reads the CSR count matrix in batches to avoid loading the full ~10 GB array.
 Optionally subsets to HVGs (var/highly_variable) and/or downsamples cells."

suppressPackageStartupMessages({
    library(optparse)
    library(Matrix)
    library(hdf5r)
    library(Seurat)
})

option_list <- list(
    make_option("--adata",     type="character", help="Path to adata_clustered.h5ad"),
    make_option("--seed",      type="integer",   default=42L,    help="Random seed"),
    make_option("--max_cells", type="integer",   default=50000L, help="Downsample to N cells (0=disabled)"),
    make_option("--all_genes", action="store_true", default=FALSE, help="Keep all genes instead of HVGs only")
)
args <- parse_args(OptionParser(option_list=option_list))
set.seed(args$seed)

cat("[SEURAT_CONVERT] Opening", args$adata, "\n")
f <- H5File$new(args$adata, mode="r")
on.exit(f$close_all(), add=TRUE)

# ── Gene names and gene subset mask ──────────────────────────────────────────
var_names <- f[["var/_index"]]$read()
n_genes   <- length(var_names)

if (args$all_genes) {
    cat(sprintf("[SEURAT_CONVERT] Keeping all %d genes\n", n_genes))
    hvg_flag  <- rep(TRUE, n_genes)
} else {
    hvg_flag  <- tryCatch(as.logical(f[["var/highly_variable"]]$read()), error=function(e) NULL)
    if (is.null(hvg_flag) || sum(hvg_flag) == 0) {
        cat("[SEURAT_CONVERT] WARNING: var/highly_variable not found — using all genes\n")
        hvg_flag  <- rep(TRUE, n_genes)
    }
}
hvg_idx   <- which(hvg_flag)
hvg_names <- var_names[hvg_idx]
n_hvg     <- length(hvg_idx)
cat(sprintf("[SEURAT_CONVERT] %d / %d genes selected\n", n_hvg, n_genes))

# fast lookup: gene_index (1-based) → output row (1-based, 0=excluded)
hvg_lookup            <- integer(n_genes)
hvg_lookup[hvg_idx]   <- seq_len(n_hvg)

# ── Cell barcodes ─────────────────────────────────────────────────────────────
obs_raw   <- f[["obs/_index"]]$read()
obs_names <- make.unique(obs_raw)
n_cells   <- length(obs_names)

# Optional downsampling
max_cells <- args$max_cells
if (max_cells > 0L && n_cells > max_cells) {
    keep_cells <- sort(sample(n_cells, max_cells))
    obs_names  <- obs_names[keep_cells]
    n_cells_sub <- max_cells
    cat(sprintf("[SEURAT_CONVERT] Downsampled: %d / %d cells kept\n", n_cells_sub, n_cells))
} else {
    keep_cells  <- seq_len(n_cells)
    n_cells_sub <- n_cells
}
cat(sprintf("[SEURAT_CONVERT] Building matrix: %d genes × %d cells\n", n_hvg, n_cells_sub))

# ── Batch-read count matrix (CSR, rows=cells) ─────────────────────────────────
# Read 500 original cells at a time so peak memory per batch ≈ 30 MB.
# For each batch, filter to selected cells AND HVG genes before accumulation.
indptr    <- f[["layers/counts/indptr"]]$read()   # length n_cells+1, cheap
data_ds   <- f[["layers/counts/data"]]
idx_ds    <- f[["layers/counts/indices"]]

batch_size_orig <- 500L
n_orig_batches  <- ceiling(n_cells / batch_size_orig)
results <- vector("list", n_orig_batches)

for (b in seq_len(n_orig_batches)) {
    orig_start <- (b - 1L) * batch_size_orig + 1L
    orig_end   <- min(b * batch_size_orig, n_cells)

    # Which selected cells fall in this batch of original cells?
    sel_in_batch <- which(keep_cells >= orig_start & keep_cells <= orig_end)
    if (length(sel_in_batch) == 0L) next

    # Contiguous read for ALL original cells in [orig_start, orig_end]
    ptr_s <- indptr[orig_start]
    ptr_e <- indptr[orig_end + 1L]
    n_ent <- ptr_e - ptr_s
    if (n_ent == 0L) next

    gene_idx_raw <- idx_ds[(ptr_s + 1L):(ptr_s + n_ent)] + 1L
    vals_raw     <- data_ds[(ptr_s + 1L):(ptr_s + n_ent)]

    # Assign each entry to its original cell (1-based)
    cell_entry_counts <- diff(indptr[orig_start:(orig_end + 1L)])
    entry_cell        <- rep(orig_start:orig_end, cell_entry_counts)

    # Mask: entry belongs to a SELECTED cell AND is an HVG gene
    hvg_row    <- hvg_lookup[gene_idx_raw]         # 1-based HVG row, 0=not HVG
    sel_sub_j  <- match(entry_cell, keep_cells)     # 1-based subsampled col, NA=not selected
    keep       <- (hvg_row > 0L) & !is.na(sel_sub_j)

    if (any(keep)) {
        results[[b]] <- list(
            rows = hvg_row[keep],
            cols = sel_sub_j[keep],
            vals = as.numeric(vals_raw[keep])
        )
    }

    if (b %% 50L == 0L)
        cat(sprintf("[SEURAT_CONVERT]   batch %d / %d\n", b, n_orig_batches))
}

cat("[SEURAT_CONVERT] Assembling sparse matrix …\n")
all_rows <- unlist(lapply(results, `[[`, "rows"))
all_cols <- unlist(lapply(results, `[[`, "cols"))
all_vals <- unlist(lapply(results, `[[`, "vals"))

counts_mat <- sparseMatrix(
    i    = all_rows,
    j    = all_cols,
    x    = all_vals,
    dims = c(n_hvg, n_cells_sub),
    dimnames = list(hvg_names, obs_names)
)
rm(all_rows, all_cols, all_vals, results); gc()

# ── Observation metadata ──────────────────────────────────────────────────────
cat("[SEURAT_CONVERT] Reading obs metadata …\n")
obs_grp  <- f[["obs"]]
obs_keys <- obs_grp$names

read_obs_col <- function(key) {
    tryCatch({
        vals <- obs_grp[[key]]$read()
        cats_path <- paste0(key, "/__categories")
        if (obs_grp$exists(cats_path))
            vals <- obs_grp[[cats_path]]$read()[vals + 1L]
        vals
    }, error = function(e) NULL)
}

meta_keys <- Filter(
    function(k) !k %in% c("_index") && !endsWith(k, "/__categories"),
    obs_keys
)
obs_list <- setNames(lapply(meta_keys, read_obs_col), meta_keys)
obs_list <- Filter(Negate(is.null), obs_list)
# Subset to selected cells
obs_list <- lapply(obs_list, `[`, keep_cells)
obs_df   <- as.data.frame(obs_list, row.names=obs_names, stringsAsFactors=FALSE)

# ── Build Seurat object ───────────────────────────────────────────────────────
cat("[SEURAT_CONVERT] Creating Seurat object …\n")
seurat_obj <- CreateSeuratObject(counts=counts_mat, meta.data=obs_df)
rm(counts_mat); gc()

# ── Reductions (obsm) ─────────────────────────────────────────────────────────
cat("[SEURAT_CONVERT] Reading obsm reductions …\n")
obsm_grp  <- f[["obsm"]]
obsm_keys <- obsm_grp$names

for (rd_key in obsm_keys) {
    emb <- obsm_grp[[rd_key]]$read()
    if (is.matrix(emb) && nrow(emb) != n_cells) emb <- t(emb)
    if (!is.matrix(emb)) emb <- matrix(emb, nrow=n_cells)
    emb <- emb[keep_cells, , drop=FALSE]          # subset to selected cells
    rownames(emb) <- obs_names
    rd_name <- sub("^X_", "", rd_key)
    key_str <- paste0(rd_name, "_")
    colnames(emb) <- paste0(key_str, seq_len(ncol(emb)))
    seurat_obj[[rd_name]] <- CreateDimReducObject(
        embeddings = emb, key = key_str, assay = "RNA"
    )
    cat(sprintf("[SEURAT_CONVERT]   added reduction '%s' (%d dims)\n", rd_name, ncol(emb)))
}

# ── Summary & save ────────────────────────────────────────────────────────────
cat("[SEURAT_CONVERT] Seurat object summary:\n")
print(seurat_obj)
cat("[SEURAT_CONVERT] Meta.data columns:", paste(colnames(seurat_obj@meta.data), collapse=", "), "\n")

saveRDS(seurat_obj, "adata_clustered.rds")
cat(sprintf("[SEURAT_CONVERT] Saved adata_clustered.rds (%d genes × %d cells)\n",
            n_hvg, n_cells_sub))
cat("[SEURAT_CONVERT] NOTE: Call NormalizeData() before downstream Seurat workflows.\n")
