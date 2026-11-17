#!/usr/bin/env python3
"""Neighbor graph, UMAP, and Leiden clustering on the Harmony-integrated object."""

import argparse
import warnings

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import scanpy as sc

warnings.filterwarnings("ignore")
sc.settings.verbosity = 1

UMAP_COLORS = ["sample", "condition", "chemistry", "leiden"]


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--adata",             required=True)
    p.add_argument("--leiden_resolution", type=float, default=0.5)
    p.add_argument("--seed",              type=int,   default=42)
    return p.parse_args()


def save_umap(adata, color, filename):
    fig, ax = plt.subplots(figsize=(7, 6))
    sc.pl.umap(adata, color=color, ax=ax, show=False, frameon=False)
    fig.tight_layout()
    fig.savefig(filename, dpi=150)
    plt.close(fig)


def main():
    args = parse_args()

    adata = sc.read_h5ad(args.adata)

    sc.pp.neighbors(adata, use_rep="X_pca_harmony", random_state=args.seed)
    sc.tl.umap(adata, random_state=args.seed)
    sc.tl.leiden(
        adata,
        resolution=args.leiden_resolution,
        random_state=args.seed,
        key_added="leiden",
    )

    for color in UMAP_COLORS:
        save_umap(adata, color, f"umap_{color}.png")

    n_clusters = adata.obs["leiden"].nunique()
    adata.write_h5ad("adata_clustered.h5ad", compression="gzip")
    print(f"[CLUSTER] {n_clusters} Leiden clusters at resolution "
          f"{args.leiden_resolution} | {adata.n_obs} cells")


if __name__ == "__main__":
    main()
