#!/usr/bin/env python3
"""Per-sample normalization, HVG selection, and scaling for scRNA-seq."""

import argparse
import warnings

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import scanpy as sc

warnings.filterwarnings("ignore")
sc.settings.verbosity = 1


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--sample", required=True)
    p.add_argument("--adata",  required=True)
    p.add_argument("--n_hvgs", type=int, default=2000)
    p.add_argument("--seed",   type=int, default=42)
    return p.parse_args()


def plot_hvg(adata, sample):
    sc.pl.highly_variable_genes(adata, show=False)
    plt.title(f"{sample} — Highly Variable Genes")
    plt.tight_layout()
    plt.savefig(f"{sample}_hvg.png", dpi=150)
    plt.close("all")


def main():
    args = parse_args()

    adata = sc.read_h5ad(args.adata)

    # Preserve raw counts before normalisation (needed for pseudobulk DE later)
    adata.layers["counts"] = adata.X.copy()

    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)

    sc.pp.highly_variable_genes(
        adata,
        n_top_genes=args.n_hvgs,
        flavor="seurat",
        batch_key=None,
    )

    plot_hvg(adata, args.sample)

    # Scaling is deferred to integrate.py (applied after merge on HVG subset)
    # to avoid converting sparse matrices to dense here and inflating file size.
    adata.write_h5ad("adata_norm.h5ad", compression="gzip")
    n_hvg = adata.var["highly_variable"].sum()
    print(f"[{args.sample}] Normalised: {adata.n_obs} cells | {n_hvg} HVGs selected")


if __name__ == "__main__":
    main()
