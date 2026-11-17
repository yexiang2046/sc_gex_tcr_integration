#!/usr/bin/env python3
"""Pseudobulk differential expression (HD vs PID7) per Leiden cluster and cell type."""

import argparse
import re
import warnings

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats

warnings.filterwarnings("ignore")
sc.settings.verbosity = 1


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--adata", required=True)
    p.add_argument("--seed",  type=int, default=42)
    return p.parse_args()


def safe_name(s):
    """Convert a string to a filesystem-safe token."""
    return re.sub(r"[^A-Za-z0-9_-]", "_", str(s))


def pseudobulk(adata, group_col, group_val):
    """Sum raw counts per sample for cells in a given group.

    Slices per sample before summing to avoid densifying the full group matrix.
    """
    import scipy.sparse as sp
    sub = adata[adata.obs[group_col] == group_val]

    samples = sorted(sub.obs["sample"].unique())
    bulk_rows, meta_rows = [], []
    for s in samples:
        mask = (sub.obs["sample"] == s).values
        X_s = sub.layers["counts"][mask]
        # Sum while keeping sparse — result is dense 1-D after squeeze
        row_sum = np.asarray(X_s.sum(axis=0)).flatten() if sp.issparse(X_s) \
                  else X_s.sum(axis=0)
        bulk_rows.append(row_sum)
        cond = sub.obs.loc[sub.obs["sample"] == s, "condition"].iloc[0]
        meta_rows.append({"sample": s, "condition": cond})

    counts_df = pd.DataFrame(
        np.vstack(bulk_rows),
        index=[r["sample"] for r in meta_rows],
        columns=sub.var_names,
    ).astype(int)
    meta_df = pd.DataFrame(meta_rows).set_index("sample")
    return counts_df, meta_df


def run_deseq2(counts_df, meta_df):
    counts_df = counts_df.loc[:, counts_df.sum(axis=0) >= 10]
    dds = DeseqDataSet(
        counts=counts_df,
        metadata=meta_df,
        design_factors="condition",
        ref_level=["condition", "HD"],
        refit_cooks=True,
        quiet=True,
    )
    dds.deseq2()
    stat_res = DeseqStats(dds, quiet=True)
    stat_res.summary()
    results = stat_res.results_df
    results.index.name = "gene"
    results.reset_index(inplace=True)
    return results


def volcano_plot(results_df, title, filename):
    df = results_df.dropna(subset=["log2FoldChange", "padj"])
    sig = df["padj"] < 0.05
    fig, ax = plt.subplots(figsize=(8, 6))
    ax.scatter(df.loc[~sig, "log2FoldChange"], -np.log10(df.loc[~sig, "padj"]),
               s=4, alpha=0.4, color="grey", label="n.s.")
    ax.scatter(df.loc[sig, "log2FoldChange"],  -np.log10(df.loc[sig, "padj"]),
               s=6, alpha=0.8, color="steelblue", label="padj<0.05")
    ax.axvline(0, color="black", linewidth=0.5)
    ax.axhline(-np.log10(0.05), color="red", linewidth=0.5, linestyle="--")
    ax.set_xlabel("log2 Fold Change (PID7 / HD)")
    ax.set_ylabel("-log10(padj)")
    ax.set_title(title)
    ax.legend(markerscale=2, fontsize=8)
    fig.tight_layout()
    fig.savefig(filename, dpi=150)
    plt.close(fig)


def run_group(adata, group_col, group_val, label_prefix, extra_cols=None):
    counts_df, meta_df = pseudobulk(adata, group_col, group_val)

    if meta_df["condition"].nunique() < 2:
        print(f"[DIFFEXP] {group_val}: skipping — only one condition present.")
        return None

    results = run_deseq2(counts_df, meta_df)
    if extra_cols:
        for i, (col, val) in enumerate(extra_cols.items()):
            results.insert(i, col, val)

    tag = safe_name(group_val)
    results.to_csv(f"de_{label_prefix}_{tag}.csv", index=False)
    volcano_plot(
        results,
        title=f"Volcano — {group_val} (PID7 / HD)",
        filename=f"volcano_{label_prefix}_{tag}.png",
    )
    n_sig = (results["padj"] < 0.05).sum()
    print(f"[DIFFEXP] {group_val}: {n_sig} DE genes (padj<0.05)")
    return results


def main():
    args = parse_args()
    adata = sc.read_h5ad(args.adata)

    # Use MapQuery cell types if available, otherwise fall back to CellTypist
    if "predicted.celltype.l2" in adata.obs.columns:
        ct_col = "predicted.celltype.l2"
        print("[DIFFEXP] Using MapQuery cell type annotations (predicted.celltype.l2)")
    elif "cell_type" in adata.obs.columns:
        ct_col = "cell_type"
        print("[DIFFEXP] Using CellTypist cell type annotations")
    else:
        ct_col = None
        print("[DIFFEXP] No cell type annotations found")

    # ── Per Leiden cluster ──────────────────────────────────────────────────
    print("[DIFFEXP] Running per-cluster DE …")
    for cluster in adata.obs["leiden"].cat.categories.tolist():
        extra = {"cluster": cluster}
        if "cell_type" in adata.obs.columns:
            ct = adata[adata.obs["leiden"] == cluster].obs["cell_type"].iloc[0]
            extra["cell_type"] = ct
        run_group(adata, "leiden", cluster, "cluster", extra_cols=extra)

    # ── Per cell type ────────────────────────────────────────────────────────
    if ct_col is None:
        print("[DIFFEXP] No cell type column found — skipping cell-type DE.")
    else:
        print(f"[DIFFEXP] Running per-cell-type DE using {ct_col} …")
        for cell_type in sorted(adata.obs[ct_col].unique()):
            run_group(adata, ct_col, cell_type, "celltype",
                      extra_cols={ct_col: cell_type})

    adata.write_h5ad("adata_de.h5ad")


if __name__ == "__main__":
    main()
