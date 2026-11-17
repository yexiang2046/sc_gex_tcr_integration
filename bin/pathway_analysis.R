#!/usr/bin/env Rscript
"Pathway analysis using fgsea on DE results from CSV files - separate outputs per cell type."

suppressPackageStartupMessages({
    library(optparse)
    library(fgsea)
    library(data.table)
    library(dplyr)
    library(ggplot2)
})

option_list <- list(
    make_option("--de_csvs",      type="character", help="Directory containing de_celltype_*.csv files"),
    make_option("--outdir",       type="character", default=".", help="Output directory"),
    make_option("--categories",   type="character", default="H,C2,C5", help="MSigDB categories (comma-separated)"),
    make_option("--species",      type="character", default="Homo sapiens", help="Species name"),
    make_option("--pval_thresh",  type="double", default=0.05, help="P-value threshold"),
    make_option("--top_pathways", type="integer", default=15, help="Number of top pathways to plot"),
    make_option("--seed",         type="integer", default=42L, help="Random seed")
)
args <- parse_args(OptionParser(option_list=option_list))
set.seed(args$seed)

cat("[PATHWAY] Pathway analysis from DE CSV files - separate per cell type\n")
cat("[PATHWAY] Working directory:", args$de_csvs, "\n")

# Find all DE CSV files
if (dir.exists(args$de_csvs)) {
    de_files <- list.files(args$de_csvs, pattern="de_celltype_.*\\.csv$", full.names=TRUE)
} else {
    de_files <- list.files(".", pattern="de_celltype_.*\\.csv$", full.names=TRUE)
}

cat("[PATHWAY] Found", length(de_files), "DE result files\n")

if (length(de_files) == 0) {
    cat("[PATHWAY] No DE CSV files found. Exiting.\n")
    quit(status=0)
}

# Parse categories
categories <- trimws(unlist(strsplit(args$categories, ",")))
cat("[PATHWAY] Loading MSigDB categories:", paste(categories, collapse=", "), "\n")

# Load gene sets
pathways <- list()
for (cat in categories) {
    pathways_cat <- msigdbr::msigdbr(species=args$species, category=cat) %>%
        split(., .$gs_name)
    for (i in seq_along(pathways_cat)) {
        pathways[[names(pathways_cat)[i]]] <- pathways_cat[[i]]$gene_symbol
    }
}
cat("[PATHWAY] Total pathways:", length(pathways), "\n\n")

# Process each DE file separately
for (de_file in de_files) {
    cell_type <- gsub("de_celltype_|\\.csv$", "", basename(de_file))
    cat("[PATHWAY] Processing:", cell_type, "\n")

    tryCatch({
        de_data <- read.csv(de_file, row.names="gene")

        # Create ranked list by log2FoldChange
        if ("log2FoldChange" %in% colnames(de_data)) {
            ranks <- setNames(de_data$log2FoldChange, rownames(de_data))
            ranks <- ranks[!is.na(ranks)]

            # Run fgsea
            gsea_res <- fgsea(pathways, ranks, eps=0, minSize=5, maxSize=500)
            gsea_res <- gsea_res[order(gsea_res$pval), ]

            # Remove list columns (leadingEdge) for CSV output
            list_cols <- sapply(gsea_res, is.list)
            gsea_res_csv <- as.data.frame(gsea_res)[, !list_cols, drop=FALSE]

            # Save individual CSV per cell type
            csv_file <- file.path(args$outdir, paste0("gsea_results_", cell_type, ".csv"))
            write.csv(gsea_res_csv, csv_file, row.names=FALSE)

            n_sig <- sum(gsea_res$padj < args$pval_thresh, na.rm=TRUE)
            cat(sprintf("  %d significant pathways (padj<%g)\n", n_sig, args$pval_thresh))

            # Plot top pathways if significant results exist
            top_pathways <- gsea_res %>%
                filter(padj < args$pval_thresh) %>%
                arrange(pval) %>%
                slice_head(n = args$top_pathways)

            if (nrow(top_pathways) > 0) {
                # Create plot with reordered pathways
                top_pathways <- top_pathways %>%
                    arrange(NES) %>%
                    mutate(pathway = factor(pathway, levels = pathway))

                p <- ggplot(top_pathways, aes(x=NES, y=pathway, fill=ifelse(padj < args$pval_thresh, "Sig", "NS"))) +
                    geom_col() +
                    scale_fill_manual(values=c("Sig"="steelblue", "NS"="grey80")) +
                    theme_minimal() +
                    theme(
                        axis.text.y = element_text(size=8),
                        plot.title = element_text(size=10, face="bold"),
                        legend.position = "none"
                    ) +
                    labs(
                        title = paste("Top Enriched Pathways -", cell_type),
                        x = "Normalized Enrichment Score",
                        y = "Pathway"
                    )

                plot_file <- file.path(args$outdir, paste0("gsea_top_pathways_", cell_type, ".png"))
                ggsave(plot_file, p, width=10, height=6, dpi=150)
                cat(sprintf("  Saved: %s\n", basename(plot_file)))
            }
        } else {
            cat("[PATHWAY]", cell_type, ": No log2FoldChange column found\n")
        }
    }, error = function(e) {
        cat("[PATHWAY]", cell_type, ": Error -", conditionMessage(e), "\n")
    })
}

cat("\n[PATHWAY] Done\n")
