#!/usr/bin/env python3
"""Rank genes per Leiden cluster and produce dotplot and heatmap."""

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
    p.add_argument("--adata",         required=True)
    p.add_argument("--n_top_markers", type=int, default=20)
    p.add_argument("--seed",          type=int, default=42)
    return p.parse_args()


def main():
    args = parse_args()

    adata = sc.read_h5ad(args.adata)

    sc.tl.rank_genes_groups(
        adata,
        groupby="leiden",
        method="wilcoxon",
        key_added="rank_genes_groups",
        use_raw=False,
        n_genes=args.n_top_markers,
    )

    clusters = adata.obs["leiden"].cat.categories.tolist()
    all_top_genes = []

    for cluster in clusters:
        df = sc.get.rank_genes_groups_df(
            adata,
            group=cluster,
            key="rank_genes_groups",
        ).head(args.n_top_markers)
        df.insert(0, "cluster", cluster)
        df.to_csv(f"markers_cluster_{cluster}.csv", index=False)
        all_top_genes.extend(df["names"].tolist())

    # Deduplicate while preserving order
    seen = set()
    top_genes = [g for g in all_top_genes if not (g in seen or seen.add(g))]

    sc.pl.dotplot(adata, var_names=top_genes, groupby="leiden", show=False)
    plt.savefig(f"markers_dotplot_top{args.n_top_markers}.png",
                dpi=150, bbox_inches="tight")
    plt.close("all")

    sc.pl.heatmap(adata, var_names=top_genes, groupby="leiden", show=False)
    plt.savefig(f"markers_heatmap_top{args.n_top_markers}.png",
                dpi=150, bbox_inches="tight")
    plt.close("all")

    adata.write_h5ad("adata_markers.h5ad", compression="gzip")
    print(f"[MARKERS] {len(clusters)} clusters | top {args.n_top_markers} markers each")


if __name__ == "__main__":
    main()
