# Pathway Analysis Module

Fast gene set enrichment analysis (GSEA) using **fgsea** and **msigdbr** R packages.

## Overview

The pathway analysis module performs rapid GSEA on differential expression results to identify enriched biological pathways and gene sets. It integrates with the Nextflow pipeline to process DE results and generate publication-ready visualizations. **Supports cell-type-specific pathway analysis**: run GSEA separately for each cell type to discover cell-type-restricted pathway programs.

## Input

- `adata_de.h5ad`: AnnData object from differential expression analysis with results stored in `.var`:
  - `logFC` (or specified ranking column): log2 fold-change or t-statistic
  - `adj_pval` (or specified p-value column): adjusted p-values

## Output

**Results tables:**
- `gsea_results_by_group.csv`: Combined results across all groups with columns:
  - `pathway`: Pathway name
  - `ES` / `NES`: Enrichment and normalized enrichment scores
  - `pval` / `padj`: p-value and adjusted p-value
  - `group`: Cell type or grouping variable
  - `leadingEdge`: Top contributing genes
- `gsea_results_[celltype].csv`: Per-group results (one file per cell type)

**Visualizations:**
- `gsea_ranked_[celltype].png`: Bar plots of top pathways per cell type
- `gsea_volcano_[celltype].png`: NES vs significance scatter plots per cell type
- `gsea_pathway_heatmap.png`: Comparative heatmap of pathway NES across cell types (if multiple groups)

## Configuration

Parameters in `nextflow.config` (already included):

```nextflow
params {
    pathway_group_by         = "predicted.celltype.l2"  // obs column for cell-type grouping (empty for no grouping)
    pathway_categories       = "H,C2,C5"  // MSigDB categories: H, C2, C5, C3, etc.
    pathway_species          = "Homo sapiens"
    pathway_pval_thresh      = 0.05       // adjusted p-value threshold
    pathway_top_pathways     = 15         // number of pathways to plot per group
}
```

### Grouping

- **With grouping** (default): Set `pathway_group_by = "predicted.celltype.l2"` to run pathway analysis separately for each cell type
  - Outputs separate enrichment files and plots per cell type
  - Generates a comparative heatmap across cell types
  
- **Without grouping**: Set `pathway_group_by = ""` to analyze all cells together
  - Single results table and visualization per analysis

### Supported MSigDB Categories

| Code | Name | Description |
|------|------|-------------|
| H | Hallmark | Hallmark gene sets |
| C1 | Positional | Positional gene sets |
| C2 | Curated | Curated gene sets (KEGG, Reactome, etc.) |
| C3 | Regulatory | Transcription factor targets |
| C4 | Computational | Computationally defined sets |
| C5 | GO | Gene Ontology gene sets |
| C6 | Oncogenic | Oncogenic signatures |
| C7 | Immunologic | Immunologic signatures |
| C8 | Cell Type | Cell type signatures |

## Usage

### Basic Usage (in Nextflow pipeline)

```nextflow
// Include the module
include { PATHWAY_ANALYSIS } from './modules/pathway_analysis'

workflow {
    // Run differential expression
    DIFFEXP(adata_clustered)
    
    // Run pathway analysis on DE results
    PATHWAY_ANALYSIS(DIFFEXP.out.adata)
}
```

Uncomment the DIFFEXP and PATHWAY_ANALYSIS calls in `main.nf`:

```bash
# Enable both modules in main.nf
sed -i 's|// DIFFEXP|DIFFEXP|; s|// PATHWAY_ANALYSIS|PATHWAY_ANALYSIS|' main.nf

# Build the Docker image
docker build -f containers/Dockerfile.pathway_analysis -t local/pathway_analysis:latest .

# Run the pipeline
nextflow run main.nf -profile docker --input samplesheet.csv --outdir results -resume
```

### Standalone Usage

```bash
# Build image
docker build -f containers/Dockerfile.pathway_analysis -t local/pathway_analysis:latest .

# Run with cell-type grouping
docker run --rm -v $(pwd):/data local/pathway_analysis:latest \
    pathway_analysis.R \
    --adata_de /data/adata_de.h5ad \
    --outdir /data/pathway_results \
    --group_by predicted.celltype.l2 \
    --categories "H,C2,C5" \
    --top_pathways 20

# Run on all cells (no grouping)
docker run --rm -v $(pwd):/data local/pathway_analysis:latest \
    pathway_analysis.R \
    --adata_de /data/adata_de.h5ad \
    --outdir /data/pathway_results \
    --group_by "" \
    --categories "H,C2,C5" \
    --top_pathways 20
```

## Command-line Arguments

```
--adata_de [PATH]
    Input AnnData object with DE results in .var

--outdir [PATH]
    Output directory (default: results/pathway_analysis)

--group_by [STRING]
    obs column to group by (e.g., predicted.celltype.l2). If empty, analyze all cells.
    (default: NULL)

--gene_col [STRING]
    Column name for gene identifiers (default: gene_names)

--ranking_col [STRING]
    Column name for ranking metric (default: logFC)

--pval_col [STRING]
    Column name for adjusted p-value (default: adj_pval)

--categories [STRING]
    Comma-separated MSigDB categories (default: H,C2,C5)

--species [STRING]
    Species for msigdbr (default: Homo sapiens)

--pval_thresh [DOUBLE]
    Adjusted p-value threshold for significant DE genes (default: 0.05)

--top_pathways [INT]
    Number of top pathways to plot per group (default: 15)

--seed [INT]
    Random seed for reproducibility (default: 42)
```

## Docker Image

The module uses `rocker/tidyverse:latest` with:
- R 4.x (with tidyverse/ggplot2/dplyr)
- fgsea: Fast gene set enrichment analysis
- msigdbr: MSigDB R interface
- anndata: Read H5AD files
- patchwork: Plot composition

Build locally:

```bash
cd /path/to/JunB_scRNAseq
docker build -f containers/Dockerfile.pathway_analysis -t local/pathway_analysis:latest .
```

## Example Workflow

### Overall Analysis
1. Differential Expression: Identify DE genes (all cells pooled) with log2FC and adjusted p-values
2. Pathway Enrichment: Rank genes and test against MSigDB sets
3. Visualization: Generate ranked pathway plots and heatmap
4. Interpretation: Identify key biological processes

### Cell-Type-Specific Analysis
1. Differential Expression: Compute DE genes per cell type (patient vs. healthy)
2. Pathway Enrichment (grouped): Run GSEA separately for each cell type
3. Visualization: 
   - Per-type pathway plots (ranked and volcano)
   - Cross-type comparative heatmap showing pathway specialization
4. Interpretation:
   - Which pathways are shared across cell types?
   - Which are cell-type-restricted?
   - Are immune cells selectively activated in JunB patient?

## Cell-Type-Specific Analysis Details

When `pathway_group_by` is set to a cell type column (e.g., `predicted.celltype.l2`):

- **Output**: Separate enrichment results for each cell type
  - `gsea_results_CD4_T_cells.csv`
  - `gsea_results_B_cells.csv`
  - etc.

- **Plots**: Per-type ranked pathway plots
  - `gsea_ranked_CD4_T_cells.png`
  - `gsea_volcano_CD4_T_cells.png`
  
- **Comparative Heatmap**: `gsea_pathway_heatmap.png`
  - Rows = top pathways across all cell types
  - Columns = cell types
  - Values = normalized enrichment scores (NES)
  - Identifies cell-type-specific and shared pathway programs

## Notes

- **Species**: Currently configured for "Homo sapiens"; modify `pathway_species` for other organisms
- **Gene Names**: Script expects HGNC gene symbols; ensure gene identifiers match MSigDB
- **Computational Cost**: fgsea is optimized for speed; scales well to thousands of gene sets
- **P-value Adjustment**: Uses Benjamini-Hochberg FDR correction by default

## References

- **fgsea**: Korotkevich et al. (2019). doi:10.1038/s41598-019-45872-z
- **msigdbr**: Liberzon et al. (2015). doi:10.1093/nar/gku1227
- **MSigDB**: Subramanian et al. (2005). doi:10.1073/pnas.0506580102
