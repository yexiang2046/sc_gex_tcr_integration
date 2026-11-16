#!/usr/bin/env python3
import argparse
import gc
import glob
import os
import pickle
import subprocess
import warnings

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import seaborn as sns
from ctxcore.rnkdb import FeatherRankingDatabase
from functools import partial as _partial
from pyscenic.aucell import aucell
from pyscenic.prune import df2regulons, load_motif_annotations, module2features_auc1st_impl, modules2df
from pyscenic.rss import regulon_specificity_scores
from pyscenic.utils import modules_from_adjacencies

warnings.filterwarnings("ignore")


def parse_args():
    p = argparse.ArgumentParser(description="Run pySCENIC regulon inference pipeline with MapQuery cell type labels")
    p.add_argument("--adata",           required=True,  help="Path to adata_annotated.h5ad")
    p.add_argument("--cell_types_csv",  required=False, help="Path to mapquery_labels.csv (barcode, predicted.celltype.l2) from ANNOTATE_SEURAT")
    p.add_argument("--db_dir",          required=True,  help="Directory containing cisTarget .feather databases")
    p.add_argument("--tfs",             required=True,  help="Path to TF list file (allTFs_hg38.txt)")
    p.add_argument("--motifs",          required=True,  help="Path to motif annotations .tbl file")
    p.add_argument("--n_cells",         type=int, default=5000,  help="Cells to subsample for GRN step (0=all)")
    p.add_argument("--num_workers",     type=int, default=4,  help="Parallel workers for GRNBoost2 and AUCell")
    p.add_argument("--ctx_num_workers", type=int, default=4,  help="Parallel workers for cisTarget (memory-limited)")
    p.add_argument("--n_hvgs",          type=int, default=4000,  help="HVGs to select before GRN (0=all genes)")
    p.add_argument("--top_targets",     type=int, default=100,   help="Top targets per TF to keep before cisTarget (0=all)")
    p.add_argument("--module_chunksize",type=int, default=100,   help="Dask chunk size for cisTarget")
    p.add_argument("--seed",            type=int, default=42)
    return p.parse_args()


def load_mapquery_cell_types(cell_types_csv_path):
    """Load MapQuery cell types from ANNOTATE_SEURAT's mapquery_labels.csv.

    Expected columns: barcode, predicted.celltype.l2, predicted.celltype.l2.score
    Returns a pandas Series indexed by barcode -> predicted.celltype.l2.
    """
    if not os.path.exists(cell_types_csv_path):
        print(f"[SCENIC] WARNING: Cell types CSV not found: {cell_types_csv_path}")
        return None

    try:
        df = pd.read_csv(cell_types_csv_path, index_col="barcode")
        label_col = "predicted.celltype.l2"
        if label_col not in df.columns:
            print(f"[SCENIC] WARNING: '{label_col}' not found in {cell_types_csv_path} "
                  f"(columns: {list(df.columns)})")
            return None
        print(f"[SCENIC]   Loaded {len(df)} cells with {df[label_col].nunique()} MapQuery l2 cell types")
        return df[label_col]
    except Exception as e:
        print(f"[SCENIC] ERROR loading cell types: {e}")
        return None


def load_and_prepare(h5ad_path, n_cells, seed, n_hvgs=0, label="subsample"):
    adata = sc.read_h5ad(h5ad_path)
    adata.obs_names_make_unique()
    if "counts" in adata.layers:
        adata.X = adata.layers["counts"].copy()
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    if n_cells > 0 and adata.n_obs > n_cells:
        sc.pp.subsample(adata, n_obs=n_cells, random_state=seed)
        print(f"[SCENIC] {label}: subsampled to {adata.n_obs} cells")
    if n_hvgs > 0 and adata.n_vars > n_hvgs:
        sc.pp.highly_variable_genes(adata, n_top_genes=n_hvgs, flavor="seurat")
        adata = adata[:, adata.var["highly_variable"]]
        print(f"[SCENIC] {label}: selected {adata.n_vars} HVGs")
    return adata


def to_expr_df(adata):
    X = adata.X
    if hasattr(X, "toarray"):
        X = X.toarray()
    return pd.DataFrame(X, index=adata.obs_names, columns=adata.var_names)


def main():
    args = parse_args()

    # ── Step 0: Load MapQuery cell types from CSV ────────────────────────────
    seurat_cell_types = None
    if args.cell_types_csv:
        print("[SCENIC] Step 0: Loading MapQuery cell type labels (predicted.celltype.l2)")
        seurat_cell_types = load_mapquery_cell_types(args.cell_types_csv)
    if seurat_cell_types is None:
        print("[SCENIC] NOTE: Cell types not available; RSS computation will be skipped")

    # ── Step 1: GRN inference ────────────────────────────────────────────────
    print("[SCENIC] Step 1: GRN inference with GRNBoost2")
    adata_grn = load_and_prepare(args.adata, args.n_cells, args.seed, n_hvgs=args.n_hvgs, label="GRN")
    expr_df = to_expr_df(adata_grn)
    expr_df.to_csv("expr_for_grn.tsv", sep="\t")
    del adata_grn

    subprocess.run(
        [
            "arboreto_with_multiprocessing.py",
            "expr_for_grn.tsv",
            args.tfs,
            "--method", "grnboost2",
            "--output", "adj.tsv",
            "--num_workers", str(args.num_workers),
            "--seed", str(args.seed),
        ],
        check=True,
    )
    os.remove("expr_for_grn.tsv")

    # ── Step 2: cisTarget pruning ─────────────────────────────────────────────
    print("[SCENIC] Step 2: cisTarget regulon pruning")
    adjacencies = pd.read_csv("adj.tsv", sep="\t")
    print(f"[SCENIC] GRN adjacencies before filtering: {len(adjacencies)} edges")
    if args.top_targets > 0:
        adjacencies = (
            adjacencies
            .sort_values("importance", ascending=False)
            .groupby("TF", sort=False)
            .head(args.top_targets)
            .reset_index(drop=True)
        )
        print(f"[SCENIC] GRN adjacencies after top-{args.top_targets} filter: {len(adjacencies)} edges")

    modules = list(modules_from_adjacencies(adjacencies, expr_df))
    del expr_df
    gc.collect()
    print(f"[SCENIC]   Created {len(modules)} TF modules")

    db_fnames = sorted(glob.glob(os.path.join(args.db_dir, "*.feather")))
    if not db_fnames:
        raise FileNotFoundError(f"No .feather databases found in {args.db_dir}")
    print(f"[SCENIC]   Using {len(db_fnames)} ranking database(s): {[os.path.basename(f) for f in db_fnames]}")
    dbs = [FeatherRankingDatabase(fname, name=os.path.basename(fname)) for fname in db_fnames]

    # Load motif annotations once (modules2df expects a DataFrame, not a filename).
    motif_annotations = load_motif_annotations(args.motifs)
    print(f"[SCENIC]   Loaded {len(motif_annotations)} motif annotations", flush=True)

    # Run cisTarget sequentially in-process to avoid fork()-based OOM.
    # prune2df always uses scheduler="processes" (fork), which doubles RSS from
    # the large adata loaded in Step 1.  Direct loop uses the same internals
    # without spawning any child processes.
    transform = _partial(
        modules2df,
        module2features_func=_partial(
            module2features_auc1st_impl,
            rank_threshold=1500,
            auc_threshold=0.05,
            nes_threshold=3.0,
            filter_for_annotation=True,
        ),
        weighted_recovery=False,
    )
    results = []
    chunk_size = args.module_chunksize
    for db in dbs:
        for i in range(0, len(modules), chunk_size):
            chunk = modules[i : i + chunk_size]
            results.append(transform(db, chunk, motif_annotations))
            if (i // chunk_size) % 5 == 0:
                print(f"[SCENIC]   cisTarget: processed {min(i + chunk_size, len(modules))}/{len(modules)} modules", flush=True)
    df = pd.concat(results) if results else pd.DataFrame()
    regulons = list(df2regulons(df))
    print(f"[SCENIC]   Found {len(regulons)} regulons")

    with open("regulons.pkl", "wb") as fh:
        pickle.dump(regulons, fh)

    if not regulons:
        print("[SCENIC] WARNING: No regulons found — check TF list and database compatibility")
        return

    # ── Step 3: AUCell on full dataset ───────────────────────────────────────
    print("[SCENIC] Step 3: AUCell scoring (all cells)")
    # Apply the same HVG filter used for GRN — regulons only contain HVG target
    # genes, so filtering to HVGs preserves all regulon gene rankings while
    # keeping the dense expression matrix to n_hvgs columns (~4000) rather than
    # all 27k genes.  Without this, to_expr_df converts 100k × 27k → ~21 GB.
    adata_full = load_and_prepare(args.adata, n_cells=0, seed=args.seed,
                                  n_hvgs=args.n_hvgs, label="full")

    # Overwrite adata.obs["cell_type"] with MapQuery predicted.celltype.l2 labels.
    # Build a fresh object-dtype column (mapped by barcode) rather than assigning
    # into any pre-existing categorical column — assigning new categories into a
    # Categorical raises "Cannot setitem on a Categorical with a new category".
    if seurat_cell_types is not None:
        mapped = adata_full.obs_names.map(seurat_cell_types)
        n_labeled = int(mapped.notna().sum())
        n_total = adata_full.n_obs
        if n_labeled > 0:
            adata_full.obs["cell_type"] = pd.Series(
                mapped, index=adata_full.obs_names, dtype="object"
            )
            print(f"[SCENIC]   Assigned MapQuery cell types to {n_labeled}/{n_total} cells")
        else:
            print("[SCENIC] WARNING: No cell barcode overlap between MapQuery labels and AnnData")

    expr_full = to_expr_df(adata_full)
    print(f"[SCENIC]   AUCell expression: {expr_full.shape[0]} cells × {expr_full.shape[1]} genes")

    auc_mtx = aucell(expr_full, regulons, auc_threshold=0.05, num_workers=1)
    auc_mtx.to_csv("auc_matrix.csv")
    print(f"[SCENIC]   AUC matrix: {auc_mtx.shape[0]} cells × {auc_mtx.shape[1]} regulons")

    # ── Step 4: Regulon Specificity Scores per cell type ─────────────────────
    if "cell_type" in adata_full.obs.columns and adata_full.obs["cell_type"].notna().any():
        print("[SCENIC] Step 4: Computing RSS per MapQuery cell type")
        # Restrict to labelled cells so unlabelled cells don't form a spurious
        # "nan" cell type; align the AUC matrix to the same barcodes.
        labelled = adata_full.obs["cell_type"].dropna()
        cells = auc_mtx.index.intersection(labelled.index)
        cell_types = labelled.loc[cells].astype(str)
        auc_labelled = auc_mtx.loc[cells]
        n_types = cell_types.nunique()
        print(f"[SCENIC]   Found {n_types} cell types across {len(cells)} labelled cells")

        rss = regulon_specificity_scores(auc_labelled, cell_types)
        rss.to_csv("rss_scores.csv")
        print(f"[SCENIC]   Computed RSS for {len(rss.columns)} cell types × {len(rss)} regulons")

        top_regs = (
            rss.apply(lambda x: x.nlargest(3).index.tolist(), axis=1)
            .explode()
            .dropna()
            .unique()
            .tolist()
        )
        if top_regs:
            n_ct  = len(rss.columns)
            n_reg = len(top_regs)
            fig, ax = plt.subplots(figsize=(max(8, n_ct * 0.9), max(6, n_reg * 0.35)))
            sns.heatmap(rss.T.loc[top_regs], cmap="viridis", ax=ax, linewidths=0)
            ax.set_title("Top regulons per MapQuery cell type (RSS)")
            plt.tight_layout()
            fig.savefig("rss_heatmap.png", dpi=150)
            plt.close()
            print(f"[SCENIC]   Saved rss_heatmap.png ({n_reg} regulons × {n_ct} cell types)")
    else:
        print("[SCENIC] NOTE: Seurat cell types not available — skipping RSS")

    # ── Step 5: Save enriched AnnData ────────────────────────────────────────
    print("[SCENIC] Step 5: Saving adata_scenic.h5ad")
    # Align AUC rows to adata obs order (obs names must match)
    common = adata_full.obs_names.intersection(auc_mtx.index)
    adata_full = adata_full[common]
    adata_full.obsm["X_scenic"] = auc_mtx.loc[common].values
    adata_full.uns["scenic_regulon_names"] = auc_mtx.columns.tolist()
    adata_full.write_h5ad("adata_scenic.h5ad", compression="gzip")
    print(f"[SCENIC] Done. {adata_full.n_obs} cells, {len(regulons)} regulons")


if __name__ == "__main__":
    main()
