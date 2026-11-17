#!/usr/bin/env python3
"""Cell type annotation using CellTypist with fine-grained internal over-clustering."""

import argparse
import warnings

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import celltypist
from celltypist import models

warnings.filterwarnings("ignore")
sc.settings.verbosity = 1


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--adata",              required=True)
    p.add_argument("--model",              default="Immune_All_Low.pkl")
    p.add_argument("--cells_per_cluster",  type=int, default=200)
    p.add_argument("--seed",               type=int, default=42)
    return p.parse_args()


def save_umap(adata, color, filename, **kwargs):
    sc.pl.umap(adata, color=color, show=False, frameon=False, **kwargs)
    plt.tight_layout()
    plt.savefig(filename, dpi=150, bbox_inches="tight")
    plt.close("all")


def main():
    args = parse_args()

    adata = sc.read_h5ad(args.adata)

    # CellTypist expects log1p-normalised counts (target_sum=1e4) in adata.X.
    # adata.X carries log-normalised values from normalize.py (scale was not
    # applied to the full object, only the HVG subset used for PCA).

    # Subsample per Leiden cluster to keep memory manageable.
    rng = np.random.default_rng(args.seed)
    idx = []
    for cluster in adata.obs["leiden"].cat.categories:
        positions = np.where(adata.obs["leiden"] == cluster)[0]
        n = min(args.cells_per_cluster, len(positions))
        idx.extend(rng.choice(positions, n, replace=False).tolist())
    adata_sub = adata[idx].copy()
    # Reset obs_names to unique integers — sc.concat produces duplicate barcodes
    # across samples which causes CellTypist's to_adata() to fail on reindexing.
    adata_sub.obs_names = [str(i) for i in range(adata_sub.n_obs)]
    print(f"[ANNOTATE] Subsampled {len(idx)} cells ({args.cells_per_cluster}/cluster) "
          f"from {adata.n_obs} total for CellTypist prediction")

    models.download_models(model=args.model, force_update=False)
    model = models.Model.load(model=args.model)

    # Do NOT pass over_clustering="leiden" — let CellTypist build its own
    # internal over-clustering at finer resolution, unlocking the full subtype
    # repertoire of Immune_All_Low rather than being capped at ~17 Leiden clusters.
    predictions = celltypist.annotate(
        adata_sub,
        model=model,
        majority_voting=True,
    )
    adata_sub = predictions.to_adata()

    # Transfer predictions back to full adata via original positional indices.
    # majority_voting reflects CellTypist's own fine-grained clusters.
    cell_types  = adata_sub.obs["majority_voting"].values
    conf_scores = adata_sub.obs["conf_score"].values

    adata.obs["cell_type"]            = pd.NA
    adata.obs["cell_type_confidence"] = np.nan
    ct_loc   = adata.obs.columns.get_loc("cell_type")
    conf_loc = adata.obs.columns.get_loc("cell_type_confidence")
    for sub_i, orig_i in enumerate(idx):
        adata.obs.iloc[orig_i, ct_loc]   = cell_types[sub_i]
        adata.obs.iloc[orig_i, conf_loc] = conf_scores[sub_i]

    # Fill unsampled cells using the dominant type among sampled cells in the
    # same Leiden cluster (preserves fine-grained labels rather than collapsing
    # all clusters to a single type).
    sampled_mask = adata.obs["cell_type"].notna()
    for cluster in adata.obs["leiden"].cat.categories:
        cluster_mask        = adata.obs["leiden"] == cluster
        sampled_in_cluster  = cluster_mask & sampled_mask
        unsampled_in_cluster = cluster_mask & ~sampled_mask
        if unsampled_in_cluster.sum() == 0:
            continue
        if sampled_in_cluster.sum() == 0:
            adata.obs.loc[unsampled_in_cluster, "cell_type"] = "Unassigned"
            adata.obs.loc[unsampled_in_cluster, "cell_type_confidence"] = 0.0
            continue
        dominant_type = adata.obs.loc[sampled_in_cluster, "cell_type"].value_counts().index[0]
        dominant_conf = adata.obs.loc[sampled_in_cluster, "cell_type_confidence"].median()
        adata.obs.loc[unsampled_in_cluster, "cell_type"] = dominant_type
        adata.obs.loc[unsampled_in_cluster, "cell_type_confidence"] = dominant_conf

    save_umap(adata, "cell_type",            "umap_cell_type.png", legend_loc="on data", legend_fontsize=7)
    save_umap(adata, "cell_type_confidence", "umap_cell_type_confidence.png")
    save_umap(adata, "condition",            "umap_condition.png")

    # Leiden cluster → dominant cell type mapping table
    mapping = (
        adata.obs.groupby("leiden")["cell_type"]
        .agg(lambda x: x.value_counts().index[0])
        .reset_index()
    )
    mapping.columns = ["leiden_cluster", "cell_type"]
    mapping.to_csv("cluster_cell_type_mapping.csv", index=False)

    adata.write_h5ad("adata_annotated.h5ad", compression="gzip")

    n_types = adata.obs["cell_type"].nunique()
    print(f"[ANNOTATE] {n_types} cell types assigned across {adata.obs['leiden'].nunique()} clusters")
    print(adata.obs["cell_type"].value_counts().to_string())


if __name__ == "__main__":
    main()
