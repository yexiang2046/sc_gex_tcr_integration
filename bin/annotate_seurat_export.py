#!/usr/bin/env python3
"""Merge Seurat MapQuery labels into adata."""

import argparse
import pandas as pd
import scanpy as sc

def main():
    p = argparse.ArgumentParser()
    p.add_argument("--adata", required=True, help="Original adata_annotated.h5ad")
    p.add_argument("--mapquery_labels", required=True, help="mapquery_labels.csv from Seurat")
    args = p.parse_args()

    adata = sc.read_h5ad(args.adata)
    labels_df = pd.read_csv(args.mapquery_labels)
    labels_df.set_index("barcode", inplace=True)

    # Transfer MapQuery cell type labels
    adata.obs["predicted.celltype.l2"] = adata.obs_names.map(
        lambda bc: labels_df.loc[bc, "predicted.celltype.l2"]
        if bc in labels_df.index else None
    )
    adata.obs["predicted.celltype.l2.score"] = adata.obs_names.map(
        lambda bc: labels_df.loc[bc, "predicted.celltype.l2.score"]
        if bc in labels_df.index else None
    )

    adata.write_h5ad("adata_annotated_seurat.h5ad", compression="gzip")
    print(f"[ANNOTATE_SEURAT_EXPORT] Saved adata_annotated_seurat.h5ad with MapQuery labels")
    print(adata.obs["predicted.celltype.l2"].value_counts())

if __name__ == "__main__":
    main()
