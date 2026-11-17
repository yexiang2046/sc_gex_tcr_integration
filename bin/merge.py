#!/usr/bin/env python3
"""Merge all per-sample normalised AnnData objects into a single object."""

import argparse
import json
import warnings

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import scanpy as sc

warnings.filterwarnings("ignore")
sc.settings.verbosity = 1

PID_PATIENT = "PID7"


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--metas",  required=True,
                   help="JSON string: list of meta dicts (one per sample)")
    p.add_argument("--adatas", required=True, nargs="+",
                   help="Paths to adata_norm.h5ad files (same order as metas)")
    return p.parse_args()


def plot_sample_composition(adata):
    counts = adata.obs["sample"].value_counts().sort_index()
    fig, ax = plt.subplots(figsize=(max(6, len(counts) * 0.8), 5))
    ax.bar(counts.index, counts.values, color="steelblue", edgecolor="white")
    ax.set_xlabel("Sample")
    ax.set_ylabel("Number of cells")
    ax.set_title("Sample composition after merge")
    ax.tick_params(axis="x", rotation=45)
    fig.tight_layout()
    fig.savefig("sample_composition.png", dpi=150)
    plt.close(fig)


def main():
    args = parse_args()
    metas = json.loads(args.metas)

    adatas = []
    for meta, path in zip(metas, args.adatas):
        ad = sc.read_h5ad(path)
        ad.obs["sample"]     = meta["sample"]
        ad.obs["patient_id"] = meta["patient_id"]
        ad.obs["chemistry"]  = meta["chemistry"]
        ad.obs["condition"]  = "PID7" if meta["patient_id"] == PID_PATIENT else "HD"
        adatas.append(ad)

    adata = sc.concat(
        adatas,
        label="sample",
        keys=[m["sample"] for m in metas],
        join="outer",
        merge="unique",
    )

    # Batch-aware HVG selection on the merged object
    sc.pp.highly_variable_genes(
        adata,
        n_top_genes=2000,
        flavor="seurat",
        batch_key="sample",
        subset=False,
    )

    plot_sample_composition(adata)

    adata.write_h5ad("adata_merged.h5ad", compression="gzip")
    print(f"[MERGE] {adata.n_obs} cells × {adata.n_vars} genes")
    print(adata.obs.groupby(["sample", "condition"])["chemistry"].first().to_string())


if __name__ == "__main__":
    main()
