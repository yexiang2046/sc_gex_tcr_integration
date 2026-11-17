# Reproducibility Manifest

**Date**: 2026-07-01  
**Commitment**: All container images built from these Dockerfiles use pinned dependency versions for complete reproducibility.

## Build Instructions

Rebuild all images with the same versions using:

```bash
# Python-based modules (Scanpy/scverse)
docker build -f containers/Dockerfile.qc -t local/qc:latest .
docker build -f containers/Dockerfile.normalize -t local/normalize:latest .
docker build -f containers/Dockerfile.merge -t local/merge:latest .
docker build -f containers/Dockerfile.integrate -t local/integrate:latest .
docker build -f containers/Dockerfile.cluster -t local/cluster:latest .
docker build -f containers/Dockerfile.markers -t local/markers:latest .
docker build -f containers/Dockerfile.annotate -t local/annotate:latest .
docker build -f containers/Dockerfile.diffexp -t local/diffexp:latest .
docker build -f containers/Dockerfile.tcr_integrate -t local/tcr_integrate:latest .

# Python-based SCENIC
docker build -f containers/Dockerfile.scenic -t local/scenic:latest .

# R-based modules (Seurat)
docker build -f containers/Dockerfile.seurat_convert -t local/seurat_convert:latest .
docker build -f containers/Dockerfile.annotate_seurat -t xiang2019/seurat_convert:v1.0.0 .
docker build -f containers/Dockerfile.pathway_analysis -t local/pathway_analysis:latest .
```

## Version Manifest

### Python Modules (Python 3.11.9-slim)

| Module | Python | scanpy | anndata | pandas | scipy | numpy | matplotlib | seaborn | Other |
|--------|--------|--------|---------|--------|-------|-------|-----------|---------|-------|
| QC | 3.11.9 | 1.10.3 | 0.10.9 | 2.2.0 | 1.13.0 | 1.24.3 | 3.8.4 | 0.13.2 | scrublet 0.2.3 |
| NORMALIZE | 3.11.9 | 1.10.3 | 0.10.9 | 2.2.0 | 1.13.0 | 1.24.3 | 3.8.4 | — | — |
| MERGE | 3.11.9 | 1.10.3 | 0.10.9 | 2.2.0 | 1.13.0 | 1.24.3 | 3.8.4 | 0.13.2 | — |
| INTEGRATE | 3.11.9 | 1.10.3 | 0.10.9 | 2.2.0 | 1.13.0 | 1.24.3 | 3.8.4 | — | harmonypy 0.0.9 |
| CLUSTER | 3.11.9 | 1.10.3 | 0.10.9 | 2.2.0 | 1.13.0 | 1.24.3 | 3.8.4 | 0.13.2 | leidenalg 0.10.2, igraph 0.11.6 |
| MARKERS | 3.11.9 | 1.10.3 | 0.10.9 | 2.2.0 | 1.13.0 | 1.24.3 | 3.8.4 | 0.13.2 | leidenalg 0.10.2, igraph 0.11.6 |
| ANNOTATE | 3.11.9 | 1.10.3 | 0.10.9 | 2.2.0 | 1.13.0 | 1.24.3 | 3.8.4 | — | celltypist 1.6.3 |
| DIFFEXP | 3.11.9 | 1.10.3 | 0.10.9 | 2.2.0 | 1.13.0 | 1.24.3 | 3.8.4 | — | pydeseq2 0.4.12 |
| TCR_INTEGRATE | 3.11.9 | 1.10.3 | 0.10.9 | 2.2.0 | 1.13.0 | 1.24.3 | 3.8.4 | 0.13.2 | scirpy 0.17.2, mudata 0.3.1, awkward 2.6.5, airr 1.5.1, setuptools 70.3.0 |

### SCENIC (Python 3.10.14-slim)

| Package | Version |
|---------|---------|
| Python | 3.10.14 |
| dask | 2023.12.1 |
| distributed | 2023.12.1 |
| pyscenic | 0.12.1 |
| loompy | 3.1.0 |
| scanpy | 1.10.3 |
| anndata | 0.10.9 |
| pandas | 2.2.0 |
| scipy | 1.13.0 |
| matplotlib | 3.8.4 |
| seaborn | 0.13.2 |

### R Modules (rocker/r-ver:4.4.0)

| Module | R | Seurat | ggplot2 | dplyr | Bioc Version | Key Packages |
|--------|---|--------|---------|-------|--------------|--------------|
| SEURAT_CONVERT | 4.4.0 | 5.1.0 | — | 1.1.4 | 3.18 | Matrix 1.6-5, hdf5r 1.3.8, SeuratDisk |
| ANNOTATE_SEURAT | 4.4.0 | 5.1.0 | 3.5.0 | 1.1.4 | 3.18 | Matrix 1.6-5, hdf5r 1.3.8, presto, SeuratDisk, glmGamPoi |
| PATHWAY_ANALYSIS | 4.4.0 | — | 3.5.0 | 1.1.4 | 3.18 | optparse 1.7.3, data.table 1.14.10, fgsea, msigdbr, rhdf5 |

### System Libraries (Debian 11)

| Library | Version | Used By |
|---------|---------|---------|
| libhdf5-dev | 1.10.6+repack-5+deb11u1 | All R modules |
| libssl-dev | 1.1.1n-0+deb11u6 | All R modules |
| libcurl4-openssl-dev | 7.74.0-1.3+deb11u11 | All R modules |
| libxml2-dev | 2.9.10+dfsg-6.7+deb11u7 | All R modules |
| libpng-dev | 1.6.37-3 | Seurat modules |
| libboost-all-dev | 1.74.0-9 | Seurat modules |
| procps | 3.3.17-5 | All modules (Nextflow task metrics) |

## Reproducibility Notes

1. **Python ecosystem**: scipy 1.13.0 + numpy 1.24.3 + pandas 2.2.0 are verified compatible with scanpy 1.10.3 and anndata 0.10.9
2. **R ecosystem**: Seurat 5.1.0 requires Matrix 1.6-5 for stability; remotes::install_version() pins exact CRAN versions
3. **Bioconductor**: Version 3.18 matches R 4.4.0; used for fgsea, msigdbr, presto, SeuratDisk, glmGamPoi, rhdf5
4. **System libraries**: Pinned to Debian 11 stable repository versions; ensure `apt-get update` is called before install
5. **Base images**: 
   - `python:3.11.9-slim` — Debian 11 base with Python 3.11.9
   - `python:3.10.14-slim` — Debian 11 base with Python 3.10.14 (SCENIC compatibility)
   - `rocker/r-ver:4.4.0` — Built on Debian 11 with R 4.4.0 from CRAN

## Cache Management

To ensure reproducible builds **without cache artifacts**:

```bash
docker build --no-cache -f containers/Dockerfile.qc -t local/qc:latest .
```

This forces re-download of all packages but guarantees the pinned versions.

## Version Upgrade Path

To upgrade a single module's dependencies in the future:

1. Test new versions locally (e.g., `pip install scanpy==1.11.0`)
2. Update the Dockerfile with new pinned versions
3. Rebuild with `--no-cache`
4. Verify pipeline runs end-to-end with new versions
5. Update this manifest with new version numbers
6. Document breaking changes in git commit message

## Validation

Each Dockerfile includes a `LABEL reproducible="2026-07-01"` tag to identify images built with this manifest:

```bash
docker image inspect local/qc:latest | grep -i reproducible
```
