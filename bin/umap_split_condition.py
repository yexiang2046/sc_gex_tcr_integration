#!/usr/bin/env python3
"""Plot UMAP split into one panel per condition value.

Each panel shows cells of that condition in the foreground (coloured)
and all other cells as a grey background, preserving spatial context.
"""

import argparse

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import scanpy as sc


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--adata",         required=True,            help="Path to adata_clustered.h5ad")
    p.add_argument("--condition_col", default="condition",       help="obs column with condition labels")
    p.add_argument("--out",           default="umap_condition_split.png", help="Output PNG path")
    return p.parse_args()


def main():
    args = parse_args()

    adata = sc.read_h5ad(args.adata)
    if args.condition_col not in adata.obs.columns:
        raise ValueError(f"Column '{args.condition_col}' not found in obs. "
                         f"Available: {list(adata.obs.columns)}")

    conditions = sorted(adata.obs[args.condition_col].astype(str).unique())
    palette    = sc.pl.palettes.default_20
    xy         = adata.obsm["X_umap"]

    fig, axes = plt.subplots(1, len(conditions), figsize=(7 * len(conditions), 6))
    if len(conditions) == 1:
        axes = [axes]

    for ax, (i, cond) in zip(axes, enumerate(conditions)):
        mask = adata.obs[args.condition_col].astype(str) == cond

        # grey background: all other cells
        ax.scatter(
            xy[~mask, 0], xy[~mask, 1],
            c="#CCCCCC", s=1, alpha=0.3, linewidths=0, rasterized=True,
        )
        # foreground: this condition
        ax.scatter(
            xy[mask, 0], xy[mask, 1],
            c=palette[i % len(palette)], s=1, alpha=0.6, linewidths=0, rasterized=True,
        )

        ax.set_title(f"{cond}  (n={mask.sum():,})", fontsize=12)
        ax.set_xlabel("UMAP 1")
        ax.set_ylabel("UMAP 2")
        ax.set_frame_on(False)

    fig.tight_layout()
    fig.savefig(args.out, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"[UMAP_SPLIT] Saved {args.out}  ({len(conditions)} panels: {', '.join(conditions)})")


if __name__ == "__main__":
    main()
