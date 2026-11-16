#!/usr/bin/env python3
"""Per-sample QC, filtering, and doublet detection for scRNA-seq Cell Ranger outputs."""

import argparse
import os
import warnings

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd
import scanpy as sc
import scrublet as scr

warnings.filterwarnings("ignore")
sc.settings.verbosity = 1


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--sample",        required=True)
    p.add_argument("--cellranger_dir", required=True)
    p.add_argument("--chemistry",     required=True, choices=["HT", "GEMX"])
    p.add_argument("--min_genes",     type=int, default=200)
    p.add_argument("--max_genes",     type=int, default=6000)
    p.add_argument("--mt_threshold",  type=float, default=15.0)
    p.add_argument("--min_cells",     type=int, default=3)
    p.add_argument("--seed",          type=int, default=42)
    return p.parse_args()


def load_h5(cellranger_dir, chemistry):
    h5_name = (
        "sample_filtered_feature_bc_matrix.h5"
        if chemistry == "GEMX"
        else "filtered_feature_bc_matrix.h5"
    )
    h5_path = os.path.join(cellranger_dir, h5_name)
    if not os.path.exists(h5_path):
        raise FileNotFoundError(f"Expected h5 file not found: {h5_path}")
    return sc.read_10x_h5(h5_path)


def plot_violin_pre(adata, sample):
    fig, axes = plt.subplots(1, 3, figsize=(12, 4))
    metrics = ["n_genes_by_counts", "total_counts", "pct_counts_mt"]
    labels  = ["Genes per cell", "Total counts", "% Mitochondrial"]
    for ax, metric, label in zip(axes, metrics, labels):
        ax.violinplot(adata.obs[metric], showmedians=True)
        ax.set_title(label)
        ax.set_xticks([])
    fig.suptitle(f"{sample} — pre-filter QC")
    fig.tight_layout()
    fig.savefig(f"{sample}_qc_prefilter_violin.png", dpi=150)
    plt.close(fig)


def plot_scatter(adata, sample):
    fig, axes = plt.subplots(1, 2, figsize=(10, 4))
    axes[0].scatter(adata.obs["total_counts"], adata.obs["n_genes_by_counts"],
                    s=1, alpha=0.3)
    axes[0].set_xlabel("Total counts")
    axes[0].set_ylabel("Genes per cell")
    axes[1].scatter(adata.obs["total_counts"], adata.obs["pct_counts_mt"],
                    s=1, alpha=0.3)
    axes[1].set_xlabel("Total counts")
    axes[1].set_ylabel("% Mitochondrial")
    fig.suptitle(f"{sample} — QC scatter")
    fig.tight_layout()
    fig.savefig(f"{sample}_qc_scatter.png", dpi=150)
    plt.close(fig)


def plot_violin_post(adata, sample):
    fig, axes = plt.subplots(1, 3, figsize=(12, 4))
    metrics = ["n_genes_by_counts", "total_counts", "pct_counts_mt"]
    labels  = ["Genes per cell", "Total counts", "% Mitochondrial"]
    for ax, metric, label in zip(axes, metrics, labels):
        ax.violinplot(adata.obs[metric], showmedians=True)
        ax.set_title(label)
        ax.set_xticks([])
    fig.suptitle(f"{sample} — post-filter QC")
    fig.tight_layout()
    fig.savefig(f"{sample}_qc_postfilter_violin.png", dpi=150)
    plt.close(fig)


def run_scrublet(adata, sample, seed):
    import scipy.sparse as sp
    X = adata.X.toarray() if sp.issparse(adata.X) else adata.X
    scrub = scr.Scrublet(X, random_state=seed)
    scores, predicted = scrub.scrub_doublets(verbose=False)
    adata.obs["doublet_score"]     = scores
    adata.obs["predicted_doublet"] = predicted

    fig, ax = plt.subplots(figsize=(6, 4))
    ax.hist(scores, bins=50, color="steelblue", edgecolor="white")
    ax.axvline(scrub.threshold_, color="red", linestyle="--",
               label=f"threshold={scrub.threshold_:.3f}")
    ax.set_xlabel("Doublet score")
    ax.set_ylabel("Cells")
    ax.set_title(f"{sample} — Scrublet doublet scores")
    ax.legend()
    fig.tight_layout()
    fig.savefig(f"{sample}_scrublet.png", dpi=150)
    plt.close(fig)

    return predicted.sum()


def main():
    args = parse_args()

    adata = load_h5(args.cellranger_dir, args.chemistry)
    adata.var_names_make_unique()
    n_cells_raw = adata.n_obs

    adata.obs["sample"]    = args.sample
    adata.obs["chemistry"] = args.chemistry

    adata.var["mt"] = adata.var_names.str.startswith("MT-")
    sc.pp.calculate_qc_metrics(
        adata, qc_vars=["mt"], percent_top=None, log1p=False, inplace=True
    )

    plot_violin_pre(adata, args.sample)
    plot_scatter(adata, args.sample)

    # Cell filtering
    adata = adata[
        (adata.obs["n_genes_by_counts"] >= args.min_genes) &
        (adata.obs["n_genes_by_counts"] <= args.max_genes) &
        (adata.obs["pct_counts_mt"]     <= args.mt_threshold)
    ].copy()

    # Gene filtering
    sc.pp.filter_genes(adata, min_cells=args.min_cells)

    n_cells_filtered = adata.n_obs

    # Doublet detection
    n_doublets = run_scrublet(adata, args.sample, args.seed)
    adata = adata[~adata.obs["predicted_doublet"]].copy()

    plot_violin_post(adata, args.sample)

    # Summary metrics
    metrics = pd.DataFrame([{
        "sample":             args.sample,
        "chemistry":          args.chemistry,
        "n_cells_raw":        n_cells_raw,
        "n_cells_filtered":   n_cells_filtered,
        "n_doublets_removed": int(n_doublets),
        "n_cells_final":      adata.n_obs,
        "n_genes":            adata.n_vars,
        "median_genes":       adata.obs["n_genes_by_counts"].median(),
        "median_counts":      adata.obs["total_counts"].median(),
        "pct_mt_median":      adata.obs["pct_counts_mt"].median(),
    }])
    metrics.to_csv("qc_metrics.csv", index=False)

    adata.write_h5ad("adata_qc.h5ad", compression="gzip")
    print(f"[{args.sample}] Raw: {n_cells_raw} → Filtered: {n_cells_filtered} "
          f"→ After doublet removal: {adata.n_obs} cells | {adata.n_vars} genes")


if __name__ == "__main__":
    main()
