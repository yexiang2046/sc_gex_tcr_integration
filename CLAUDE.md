# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

Single-cell RNA-seq (scRNA-seq) analysis of **JunB patient vs. healthy control** PBMCs. Analysis starts from pre-computed 10x Genomics Cell Ranger outputs — no FASTQ/BCL processing required.

The project uses a **custom Nextflow DSL2 pipeline** orchestrating Scanpy/scverse steps in containerized, modular processes.

**Experimental design**: 7 healthy controls (HD) + 1 JunB patient (PID7), two Cell Ranger chemistries:
- **HT**: samples 01–03, 09–11, 16 — outputs at `data/count_12113-XZ-<id>/`
- **GEMX**: samples 17–20 — outputs at `data/GEMX/count_12113-XZ-<id>_GEMX/per_sample_outs/count_12113-XZ-<id>_GEMX/count/`

HT and GEMX chemistries use **different h5 filenames**:
- HT: `filtered_feature_bc_matrix.h5`
- GEMX: `sample_filtered_feature_bc_matrix.h5`

## Directory Structure

```
project-root/
├── data/                             # Cell Ranger outputs (not committed)
├── samplesheet.csv                   # sample → path mapping
├── main.nf                           # pipeline entrypoint
├── nextflow.config                   # parameters & profiles
├── bin/                              # scripts copied into Docker PATH by Nextflow
│   ├── qc.py                         # QC module script
│   └── sync_cellranger.sh            # rsync helper to pull outs from HPC
├── containers/                       # one Dockerfile per module (Dockerfile.<module>)
│   └── Dockerfile.qc
├── modules/                          # one .nf file per Nextflow process
│   └── qc.nf
└── results/                          # pipeline outputs (publishDir target)
```

## Running the Pipeline

**Before first run**, build Docker images for all modules. See [REPRODUCIBILITY.md](REPRODUCIBILITY.md) for pinned version details and full build instructions.

Quick build (all modules at once):

```bash
for module in qc normalize merge integrate cluster markers annotate seurat_convert scenic pathway_analysis; do
  docker build -f containers/Dockerfile.$module -t local/$module:latest .
done
# Special: annotate_seurat uses external registry
docker build -f containers/Dockerfile.annotate_seurat -t xiang2019/seurat_convert:v1.0.0 .
```

**Run the pipeline:**

```bash
nextflow run main.nf -profile docker --input samplesheet.csv --outdir results -resume
```

**Key params** (override on CLI with `--param value`):

| Param | Default | Notes |
|---|---|---|
| `min_genes` | 200 | Lower bound for genes/cell |
| `max_genes` | 6000 | Upper bound for genes/cell |
| `mt_threshold` | 15 | Max % mitochondrial reads |
| `min_cells` | 3 | Min cells expressing a gene |
| `seed` | 42 | Fixed for reproducibility |

## Pipeline Architecture

`main.nf` reads `samplesheet.csv`, constructs a channel of `(meta, cellranger_dir)` tuples, and passes them to module processes. `meta` carries `{sample, patient_id, description, chemistry}`.

**Current modules** (each in `modules/<name>.nf` + `bin/<name>.py/R` + `containers/Dockerfile.<name>`):

- **QC** (`modules/qc.nf`, Python): loads h5 → calculates QC metrics → filters cells/genes → Scrublet doublet detection → outputs `adata_qc.h5ad` + PNG QC plots + `qc_metrics.csv`
- **NORMALIZE** (`modules/normalize.nf`, Python): log-normalization, HVG selection (6000 genes by default)
- **MERGE** (`modules/merge.nf`, Python): concatenate normalized samples across all batches
- **INTEGRATE** (`modules/integrate.nf`, Python): BBKNN batch correction + re-compute neighbors
- **CLUSTER** (`modules/cluster.nf`, Python): PCA → UMAP → Leiden clustering
- **MARKERS** (`modules/markers.nf`, Python): cell-type marker detection per cluster (scanpy.tl.rank_genes_groups)
- **ANNOTATE** (`modules/annotate.nf`, Python): manual annotation of Leiden clusters → per-sample metadata mapping
- **SEURAT_CONVERT** (`modules/seurat_convert.nf`, R): convert clustered AnnData to Seurat object for downstream R analysis
- **ANNOTATE_SEURAT** (`modules/annotate_seurat.nf`, R): MapQuery cell-type label transfer from pbmc_multimodal reference; adds `predicted.celltype.l1`/`predicted.celltype.l2` and scores to Seurat metadata. **`predicted.celltype.l2` (level 2) is the definitive cell-type annotation** used downstream by SCENIC and DIFFEXP
- **SCENIC** (`modules/scenic.nf`, Python): regulon activity inference (GRNBoost → AUC scoring). Reads MapQuery `predicted.celltype.l2` labels from ANNOTATE_SEURAT's `mapquery_labels.csv` (no R in the container — the CSV is consumed directly); outputs AUC matrix + regulons PKL + regulon specificity scores (RSS) per cell type
- **DIFFEXP** (`modules/diffexp.nf`, R): Seurat FindMarkers per cell type (PID7 vs HD); outputs `de_celltype_<name>.csv` per cell type + combined `adata_de.h5ad`
- **PATHWAY_ANALYSIS** (`modules/pathway_analysis.nf`, R): GSEA on DE results; outputs `gsea_results_<celltype>.csv` and `gsea_top_pathways_<celltype>.png` per cell type (run separately post-pipeline)

## Conventions for Adding Modules

1. Add the Python script to `bin/<module>.py` (Nextflow auto-adds `bin/` to PATH inside containers)
2. Add the Nextflow process to `modules/<module>.nf` with `container 'local/<module>:latest'`
3. Add `containers/Dockerfile.<module>` with pinned versions
4. Build the image: `docker build -f containers/Dockerfile.<module> -t local/<module>:latest .`
5. `include` the process in `main.nf` and wire it into the workflow
6. Document any new `params` with rationale in `nextflow.config`

The `Chemistry` column must drive branching logic for HT vs GEMX where behavior differs.

## Samplesheet Format

`samplesheet.csv` columns: `Sample`, `cellranger_outs`, `Patient ID`, `Description`, `Chemistry`

```csv
Sample,cellranger_outs,Patient ID,Description,Chemistry
sample1,data/count_12113-XZ-01,HD31,Healthy control,HT
sample16,data/count_12113-XZ-16,PID7,JunB patient,HT
sample17,data/GEMX/count_12113-XZ-17_GEMX/per_sample_outs/count_12113-XZ-17_GEMX/count,HD67,Healthy control,GEMX
```

## Standard Workflow Order

1. **QC** — per-sample QC, filtering, doublet detection → `adata_qc.h5ad`
2. **NORMALIZE** — log-normalization, HVG selection (6000 genes) → `adata_norm.h5ad`
3. **MERGE** — concatenate all normalized samples → `adata_merged.h5ad`
4. **INTEGRATE** — BBKNN batch correction, recompute neighbors → `adata_integrated.h5ad`
5. **CLUSTER** — PCA → UMAP → Leiden clustering → `adata_clustered.h5ad`
6. **MARKERS** — scanpy marker detection → `markers.csv`
7. **ANNOTATE** — manual cluster annotation + metadata → `adata_annotated.h5ad`
8. **SEURAT_CONVERT** — convert to Seurat object → `seurat_clustered.rds`
9. **ANNOTATE_SEURAT** — MapQuery reference-based cell type labeling → `seurat_annotated.rds` (with `predicted.celltype.l2`)
10. **SCENIC** — regulon inference (GRNBoost) + AUC scoring per MapQuery cell type → `auc_matrix.csv`, `regulons.pkl`, `rss_scores.csv` (regulon specificity per cell type), `rss_heatmap.png`
11. **DIFFEXP** — Seurat FindMarkers per cell type (PID7 vs HD) → `de_celltype_<name>.csv` (per cell type), `adata_de.h5ad`
12. **PATHWAY_ANALYSIS** (post-pipeline) — GSEA on DE results → `gsea_results_<celltype>.csv`, `gsea_top_pathways_<celltype>.png` (per cell type)

## Code Generation Rules

- Modular Nextflow DSL2 — never monolithic scripts
- Prefer Scanpy/scverse for Python; use Seurat for R-based downstream analysis
- Use relative paths + Nextflow channels; never hard-code absolute paths
- **Pin all container dependencies** — all Dockerfiles use pinned versions (see [REPRODUCIBILITY.md](REPRODUCIBILITY.md) for complete version manifest). Never use unpinned package names or `latest` tags in RUN commands.
- **Never call `sc.pp.scale()` in per-sample steps** (normalize.py). Scaling converts sparse matrices to dense float64, inflating each h5ad ~10× and causing OOM when MERGE loads all samples simultaneously. Scale only in integrate.py on the HVG subset immediately before PCA.
- All `bin/` scripts must be executable (`chmod +x`) — Nextflow copies them into containers via PATH but does not set permissions.
- All Dockerfiles need `procps` installed (for Nextflow's `ps`-based task metrics) and `ENV NUMBA_CACHE_DIR=/tmp MPLCONFIGDIR=/tmp` (containers run as the host user who cannot write to site-packages or `/.config`).
- Pin `pandas==2.2.0` (not `<3`) — anndata 0.10.x cannot write `pandas.arrays.StringArray` to HDF5; use pandas 2.2.0 which is compatible with scipy 1.13.0 and numpy 1.24.3.
- When collecting per-sample outputs into an aggregation step (MERGE), use `stageAs: "inputs/adata_norm_?.h5ad"` on the `path` input to avoid filename collisions.
- **R scripts** — start with `#!/usr/bin/env Rscript`, use `optparse` for CLI args, use `suppressPackageStartupMessages()` for cleaner output, use `remotes::install_version()` to pin exact CRAN versions in Dockerfiles.
- **Pathway analysis as post-processing** — PATHWAY_ANALYSIS runs separately via Docker on published DIFFEXP results to avoid Nextflow channel/glob pattern issues with optional file outputs
- **System libraries** — pin Debian package versions in Dockerfiles (e.g., `libhdf5-dev=1.10.6+repack-5+deb11u1`) to ensure reproducibility across rebuilds

## Seurat-Based Differential Expression & Pathway Analysis

**DIFFEXP module** (`bin/diffexp_seurat.R`):
- Input: `seurat_annotated.rds` (from ANNOTATE_SEURAT) + `adata_with_condition` (for condition metadata)
- Extracts `patient_id` from AnnData metadata (maps PID7 → "PID7", others → "HD")
- Iterates per cell type (uses `predicted.celltype.l2` from MapQuery)
- Runs Seurat `FindMarkers()` with Wilcoxon rank-sum test, comparing PID7 vs HD
- Outputs:
  - `de_celltype_<celltype>.csv` — one CSV per cell type with columns: `gene`, `p_val`, `avg_log2FC`, `pct.1`, `pct.2`, `p_val_adj`
  - `adata_de.h5ad` — combined h5ad with DE results embedded

**PATHWAY_ANALYSIS module** (`bin/pathway_analysis.R`):
- **Post-pipeline step** — run separately on published DE CSV files (not integrated into main Nextflow workflow due to Nextflow channel/glob pattern limitations with optional outputs)
- Input: directory of `de_celltype_*.csv` files
- Loads MSigDB gene sets (default: H, C2, C5 categories); customizable via `--categories` flag
- Runs `fgsea()` independently per cell type
- Outputs per cell type:
  - `gsea_results_<celltype>.csv` — pathway enrichment scores with columns: `pathway`, `pval`, `padj`, `log2err`, `ES`, `NES`, `size`
  - `gsea_top_pathways_<celltype>.png` — bar plot of top 15 NES-ranked pathways
- **NES interpretation**: positive = upregulated in PID7; negative = downregulated in PID7

**To run pathway analysis post-pipeline:**

```bash
docker run -v $(pwd):/data local/pathway_analysis:latest bin/pathway_analysis.R --de_csvs /data/results_HVG6000_HT/diffexp/ --outdir /data/results_HVG6000_HT/diffexp/ --categories H,C2,C5
```

## Helper Script

`bin/sync_cellranger.sh` — rsync Cell Ranger `outs/` from HPC, skipping BAM, raw matrices, and `molecule_info.h5`:

```bash
./bin/sync_cellranger.sh -n user@cluster:/data/run/outs/ data/count_12113-XZ-01/  # dry-run
./bin/sync_cellranger.sh user@cluster:/data/run/outs/ data/count_12113-XZ-01/     # real sync
```

## Environment

- **Nextflow**: 24.10+
- **Container engine**: Docker (preferred) or Singularity/Apptainer (HPC)
- **Conda** (for local script dev, not pipeline runs): `sc-rna-nf` with Scanpy ≥ 1.10, AnnData ≥ 0.10
- **R Docker image for Seurat**: `xiang2019/seurat_convert:v1.0.0` (ANNOTATE_SEURAT, DIFFEXP modules)
- **R Docker image for pathway analysis**: `local/pathway_analysis:latest` (built from `containers/Dockerfile.pathway_analysis`)
