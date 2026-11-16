# Quick Start Guide

Complete scRNA-seq analysis pipeline for JunB patient vs. healthy control PBMCs.

## Prerequisites

1. **Docker** — container runtime for reproducible environments
2. **Nextflow** — workflow orchestration (≥25.10)
3. **Git** — for version control
4. **macOS/Linux** — tested on Darwin 25.5.0 with zsh

### Install Nextflow

```bash
# Download and install (or update PATH if already installed)
curl -s https://get.nextflow.io | bash
mv nextflow ~/.local/bin/  # or any directory in your PATH
chmod +x ~/.local/bin/nextflow
```

### Verify Installation

```bash
docker --version
nextflow --version
```

## Quick Setup

### 1. Build Docker Images (first time only)

All images use pinned versions from [REPRODUCIBILITY.md](REPRODUCIBILITY.md) for reproducibility.

```bash
# Build all images (takes ~20-30 min depending on Internet speed)
./run_ht_pipeline.sh --build

# OR build individual images:
docker build -f containers/Dockerfile.qc -t local/qc:latest .
docker build -f containers/Dockerfile.normalize -t local/normalize:latest .
# ... (see REPRODUCIBILITY.md for full list)
```

### 2. Run Analysis

#### HT Chemistry (Samples 01–03, 09–11, 16)

```bash
# First run (builds from scratch)
./run_ht_pipeline.sh --build

# Subsequent runs (resume from last checkpoint)
./run_ht_pipeline.sh

# Or manually:
unset JAVA_CMD JAVA_HOME && nextflow run main.nf \
  -profile docker \
  --input samplesheet_HT.csv \
  --outdir results_HVG6000_HT \
  -resume
```

#### GEMX Chemistry (Samples 17–20)

```bash
# Assumes images already built (by run_ht_pipeline.sh --build)
./run_gemx_pipeline.sh
```

## Pipeline Outputs

### HT Results: `results_HVG6000_HT/`

```
├── qc/                    # QC metrics & filtering plots per sample
├── normalize/             # Normalized, HVG-selected counts
├── merge/                 # Merged across all HT samples
├── integrate/             # Batch-corrected (BBKNN)
├── cluster/               # Clustering, UMAP, Leiden communities
├── markers/               # Cluster marker genes
├── annotate/              # Manual cluster annotation
├── seurat_convert/        # Seurat object (clustered)
├── annotate_seurat/       # MapQuery cell-type labels (predicted.celltype.l2)
├── scenic/                # Regulon activity (AUC matrix)
├── diffexp/               # Differential expression per cell type
│   ├── adata_de.h5ad                    # h5ad with DE results
│   ├── de_celltype_CD4_TCM.csv          # DE genes for each cell type
│   ├── de_celltype_Plasma_cells.csv
│   ├── de_celltype_*.csv                # ... 18 more cell types
│   ├── gsea_results_CD4_TCM.csv         # Pathway enrichment (GSEA)
│   ├── gsea_results_Plasma_cells.csv
│   ├── gsea_results_*.csv               # ... per cell type
│   ├── gsea_top_pathways_CD4_TCM.png    # Visualization
│   ├── gsea_top_pathways_*.png          # ... per cell type
│   └── ...
└── tcr/                   # GEX + TCR integration (scirpy)
    ├── adata_gex_tcr.h5ad
    ├── mudata_gex_tcr.h5mu
    ├── tcr_summary.csv
    ├── clonotypes_by_sample.csv
    └── *.png
```

### TCR inputs

Primary V(D)J outs live under `data/TCR 12113/` (space in the name): `01_TCR` … `16_TCR`. HT samplesheet rows quote those paths, e.g. `"data/TCR 12113/02_TCR"`. `params.tcr_dir` defaults to that folder for the Docker bind-mount. Build the TCR image once:

```bash
docker build -f containers/Dockerfile.tcr_integrate -t local/tcr_integrate:latest .
```

Libraries without TCR (including GEMX 17–20) leave `tcr_outs` blank; the module still runs and marks those cells as `no_tcr_input`.

### GEMX Results: `results_HVG6000_GEMX/`

Same structure as HT results (separate batch).

## Interpreting Results

### Cell Type Annotations
See `results_HVG6000_HT/annotate_seurat/` for cell type labels assigned via MapQuery.

### Differential Expression
Each `de_celltype_<name>.csv` contains:
- `gene` — gene name
- `p_val` — Wilcoxon rank-sum p-value
- `avg_log2FC` — log2 fold change (PID7 vs HD)
- `pct.1` — % cells expressing in PID7
- `pct.2` — % cells expressing in HD
- `p_val_adj` — adjusted p-value (Bonferroni)

### Pathway Enrichment (GSEA)
Each `gsea_results_<celltype>.csv` contains:
- `pathway` — gene set name (from MSigDB)
- `pval`, `padj` — significance
- `NES` — Normalized Enrichment Score
  - **Positive NES** = pathway upregulated in PID7 (JunB patient)
  - **Negative NES** = pathway downregulated in PID7
- `size` — number of genes in pathway present in dataset

## Common Commands

### Resume from Last Checkpoint
```bash
./run_ht_pipeline.sh
```

### Clean and Re-run from Scratch
```bash
./run_ht_pipeline.sh --clean --no-resume --skip-build
```

### Run Pathway Analysis on Existing DE Results
```bash
docker run -v $(pwd)/results_HVG6000_HT/diffexp:/data local/pathway_analysis:latest \
  pathway_analysis.R \
  --de_csvs /data \
  --outdir /data \
  --categories "H,C2,C5"
```

### View Pipeline Logs
```bash
# Live monitoring (if pipeline is running)
tail -f .nextflow.log

# Detailed logs in work directory:
cat work/<task-hash>/command.log
```

### Check Docker Image Versions
```bash
# Verify reproducible build date
docker image inspect local/qc:latest | grep -i reproducible

# Check all pinned versions
cat REPRODUCIBILITY.md
```

## Troubleshooting

### Java Path Error
If you see: `ERROR: Cannot find Java or it's a wrong version`

```bash
# Solution: clear Java environment variables
unset JAVA_CMD JAVA_HOME
./run_ht_pipeline.sh
```

### Docker Image Not Found
```bash
# Rebuild specific image
docker build -f containers/Dockerfile.<module> -t local/<module>:latest .

# Or rebuild all
./run_ht_pipeline.sh --build
```

### Out of Memory
The pipeline uses aggressive caching to minimize memory. If OOM errors occur:

1. Ensure Docker has ≥16 GB allocated
2. Check `nextflow.config` memory limits
3. Consider processing samples in smaller batches

### Incomplete Results
Check output directory exists and is writable:
```bash
ls -lh results_HVG6000_HT/
```

## Performance

**Typical runtime:**
- Full pipeline (HT, 7 samples): **4–6 hours** (first run)
- Resumed run (cached): **30 min** (re-runs only changed steps)
- Pathway analysis: **15–30 min** (per chemistry)

**Peak memory:** ~16 GB (SCENIC + INTEGRATE steps)

**Disk space:** ~80 GB per chemistry (temp work directory + results)

## Documentation

- **[CLAUDE.md](CLAUDE.md)** — architecture, conventions, module descriptions
- **[REPRODUCIBILITY.md](REPRODUCIBILITY.md)** — version pinning, build instructions
- **[nextflow.config](nextflow.config)** — pipeline parameters
- **[main.nf](main.nf)** — workflow DAG
- **[modules/](modules/)** — individual process definitions

## Support

For issues, check:
1. `.nextflow.log` for Nextflow errors
2. `work/<task-hash>/command.log` for process-specific errors
3. `docker logs` for container issues
4. [CLAUDE.md](CLAUDE.md) for pipeline architecture

## Next Steps

After running the pipeline:

1. **Review QC metrics** — sample quality, cell counts
2. **Inspect UMAP plots** — cluster separation, cell type overlap
3. **Check DE results** — top genes per cell type, log2FC distribution
4. **Explore pathway enrichment** — biological themes in patient vs. controls
5. **Cross-reference with literature** — validate findings
