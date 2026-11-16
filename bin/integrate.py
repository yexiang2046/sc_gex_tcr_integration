#!/usr/bin/env python3
"""PCA + Harmony batch integration on the merged scRNA-seq object."""

import argparse
import warnings

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import scanpy as sc
import harmonypy as hm

warnings.filterwarnings("ignore")
sc.settings.verbosity = 1


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--adata",  required=True)
    p.add_argument("--n_hvgs", type=int, default=2000)
    p.add_argument("--seed",   type=int, default=42)
    return p.parse_args()


def plot_convergence(ho):
    obj = np.array(ho.objective_harmony)
    fig, ax = plt.subplots(figsize=(7, 4))
    ax.plot(obj, marker="o", markersize=3, linewidth=1.5)
    ax.set_xlabel("Iteration")
    ax.set_ylabel("Harmony objective")
    ax.set_title("Harmony convergence")
    fig.tight_layout()
    fig.savefig("harmony_convergence.png", dpi=150)
    plt.close(fig)


def main():
    args = parse_args()

    adata = sc.read_h5ad(args.adata)

    # PCA on HVG subset — scale here, not per-sample, to keep upstream files sparse
    adata_hvg = adata[:, adata.var["highly_variable"]].copy()
    sc.pp.scale(adata_hvg, max_value=10)
    sc.tl.pca(adata_hvg, n_comps=50, random_state=args.seed, use_highly_variable=False)

    adata.obsm["X_pca"] = adata_hvg.obsm["X_pca"]
    adata.uns["pca"]    = adata_hvg.uns["pca"]

    # Harmony integration — sample captures both patient identity and chemistry
    ho = hm.run_harmony(
        adata.obsm["X_pca"],
        adata.obs,
        vars_use=["sample"],
        random_state=args.seed,
        max_iter_harmony=20,
    )
    adata.obsm["X_pca_harmony"] = ho.Z_corr.T

    plot_convergence(ho)

    adata.write_h5ad("adata_integrated.h5ad", compression="gzip")
    print(f"[INTEGRATE] Done: {adata.n_obs} cells | "
          f"X_pca_harmony shape: {adata.obsm['X_pca_harmony'].shape}")


if __name__ == "__main__":
    main()
