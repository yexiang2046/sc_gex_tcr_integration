# sc_gex_tcr_integration

Nextflow pipeline for JunB PBMC **gene expression (GEX)** analysis with optional **TCR / V(D)J** integration.

Starts from 10x Cell Ranger count outputs (and, when available, Cell Ranger V(D)J outs). GEX steps use Scanpy/Seurat; TCR joining uses [scirpy](https://scirpy.scverse.org/).

## Quick start

See [QUICKSTART.md](QUICKSTART.md) for image builds and run commands. TCR integration uses `data/TCR 12113/{01..16}_TCR` for HT libraries (quoted paths; space in the folder name). GEMX samples 17–20 leave `tcr_outs` blank.

```bash
docker build -f containers/Dockerfile.tcr_integrate -t local/tcr_integrate:latest .
nextflow run main.nf -profile docker --input samplesheet_HT.csv --outdir results_HVG6000_HT -resume
```

Outputs land in `results*/tcr/`: `adata_gex_tcr.h5ad`, `mudata_gex_tcr.h5mu`, match/clonotype tables, and QC plots.
