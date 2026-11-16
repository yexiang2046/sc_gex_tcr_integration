#!/usr/bin/env python3
"""Integrate 10x V(D)J TCR annotations into a GEX AnnData with scirpy."""

from __future__ import annotations

import argparse
import json
import os
import warnings
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import mudata as mu
import numpy as np
import pandas as pd
import scanpy as sc
import scirpy as ir
import seaborn as sns

warnings.filterwarnings("ignore")
sc.settings.verbosity = 1

CONTIG_CANDIDATES = (
    "filtered_contig_annotations.csv",
    "all_contig_annotations.csv",
)


def parse_args():
    p = argparse.ArgumentParser(
        description="Join Cell Ranger V(D)J TCR calls onto annotated GEX AnnData"
    )
    p.add_argument("--adata", required=True, help="Annotated GEX AnnData (.h5ad)")
    p.add_argument(
        "--mapquery_labels",
        default=None,
        help="Optional mapquery_labels.csv from ANNOTATE_SEURAT",
    )
    p.add_argument(
        "--tcr_map",
        required=True,
        help='JSON list of {"sample": "...", "path": "..."} TCR outs directories',
    )
    p.add_argument("--seed", type=int, default=42)
    return p.parse_args()


def normalize_barcode(bc: str) -> str:
    """Strip GEM suffix variants so GEX and TCR barcodes can match."""
    bc = str(bc)
    # Drop AnnData/MapQuery uniquify suffixes first (e.g. AAAC-1.1 → AAAC-1).
    if "." in bc:
        bc = bc.rsplit(".", 1)[0]
    if "-" in bc:
        return bc.split("-", 1)[0]
    return bc


def ensure_unique_obs_names(adata: sc.AnnData, join: str = ".") -> int:
    """Uniquify GEX obs_names the way MapQuery does (AnnData join='.').

    Merged HT objects reuse the same 10x barcode string across samples, which
    breaks scirpy chain_qc. Preserve the original barcode in obs['barcode_raw']
    for TCR matching, then make obs_names unique in place.
    """
    if "barcode_raw" not in adata.obs.columns:
        adata.obs["barcode_raw"] = adata.obs_names.astype(str)
    n_unique = int(pd.Index(adata.obs_names).nunique())
    n_dup = adata.n_obs - n_unique
    if n_dup <= 0:
        return 0
    adata.obs_names_make_unique(join=join)
    print(
        f"[TCR] Uniquified GEX obs_names: {n_dup} collisions "
        f"({n_unique} unique strings → {adata.n_obs} cells, join={join!r})"
    )
    return n_dup


def find_contig_table(tcr_dir: Path) -> Path | None:
    for base in (tcr_dir, tcr_dir / "outs"):
        if not base.is_dir() and base != tcr_dir:
            continue
        for name in CONTIG_CANDIDATES:
            candidate = base / name if base.is_dir() else tcr_dir / name
            if candidate.is_file() and candidate.stat().st_size > 0:
                return candidate
    # Also allow tcr_dir itself to be a file path edge-case
    if tcr_dir.is_file() and tcr_dir.name in CONTIG_CANDIDATES:
        return tcr_dir
    return None


def load_mapquery_labels(path: str | None, adata: sc.AnnData) -> None:
    if not path or not os.path.exists(path):
        return
    df = pd.read_csv(path)
    if "barcode" not in df.columns or "predicted.celltype.l2" not in df.columns:
        print(
            f"[TCR] WARNING: mapquery labels missing expected columns "
            f"(got {list(df.columns)})"
        )
        return
    series = df.set_index("barcode")["predicted.celltype.l2"]
    overlap = adata.obs_names.intersection(series.index)
    if len(overlap):
        adata.obs["predicted.celltype.l2"] = series.reindex(adata.obs_names)
    else:
        norm = pd.Series(series.values, index=[normalize_barcode(b) for b in series.index])
        norm = norm[~norm.index.duplicated(keep="first")]
        adata.obs["predicted.celltype.l2"] = [
            norm.get(normalize_barcode(b), pd.NA) for b in adata.obs_names
        ]
    print(f"[TCR] Loaded MapQuery labels for {series.shape[0]} cells")


def load_sample_tcr(sample: str, tcr_path: str):
    tcr_dir = Path(tcr_path)
    contig = find_contig_table(tcr_dir)
    if contig is None:
        print(f"[TCR] {sample}: no contig annotations under {tcr_path}; skipping")
        return None, "missing"

    filtered = contig.name.startswith("filtered_")
    print(f"[TCR] {sample}: reading {contig} (filtered={filtered})")
    adata_tcr = ir.io.read_10x_vdj(str(contig), filtered=True)
    if adata_tcr.n_obs == 0:
        print(f"[TCR] {sample}: zero cells after import; skipping")
        return None, "empty"

    raw = adata_tcr.obs_names.astype(str)
    adata_tcr.obs["barcode_raw"] = raw.values
    adata_tcr.obs["barcode_norm"] = [normalize_barcode(b) for b in raw]
    adata_tcr.obs["sample"] = sample
    status = "filtered" if filtered else "all_contig"
    return adata_tcr, status


def match_airr_to_gex(adata: sc.AnnData, airr_by_sample: dict[str, sc.AnnData]):
    """Return AIRR AnnData whose obs_names match GEX cells (matched subset only).

    Matching uses obs['barcode_raw'] (pre-uniquify 10x barcode) + sample.
    AIRR obs_names are set to the current (unique) GEX obs_names so MuData /
    scirpy chain_qc see a 1:1 index.
    """
    adata.obs["has_tcr"] = False
    adata.obs["tcr_status"] = "no_tcr_input"
    if "barcode_raw" not in adata.obs.columns:
        adata.obs["barcode_raw"] = adata.obs_names.astype(str)
    adata.obs["barcode_norm"] = [
        normalize_barcode(b) for b in adata.obs["barcode_raw"].astype(str)
    ]

    if "sample" not in adata.obs.columns:
        raise ValueError("GEX AnnData must contain obs['sample'] for TCR matching")

    matched_airr = []

    for sample, adata_tcr in airr_by_sample.items():
        mask = adata.obs["sample"].astype(str).to_numpy() == sample
        if not mask.any():
            print(f"[TCR] WARNING: sample {sample} has TCR but no GEX cells")
            continue

        gex_names = adata.obs_names.to_numpy()[mask]
        gex_raw = adata.obs["barcode_raw"].to_numpy()[mask]
        gex_norm = adata.obs["barcode_norm"].to_numpy()[mask]

        # Prefer exact raw barcode, then normalized barcode
        by_raw = {
            str(raw): i for i, raw in enumerate(adata_tcr.obs["barcode_raw"].astype(str))
        }
        by_norm = {}
        for i, norm in enumerate(adata_tcr.obs["barcode_norm"].astype(str)):
            by_norm.setdefault(norm, i)

        hit_flags = np.zeros(mask.sum(), dtype=bool)
        airr_indices = []
        gex_hit_names = []
        for j, (gname, graw, gnorm) in enumerate(zip(gex_names, gex_raw, gex_norm)):
            idx = by_raw.get(str(graw))
            if idx is None:
                idx = by_norm.get(str(gnorm))
            if idx is None:
                continue
            hit_flags[j] = True
            airr_indices.append(idx)
            gex_hit_names.append(gname)

        status = np.where(hit_flags, "matched", "gex_only")
        adata.obs.loc[mask, "tcr_status"] = status
        adata.obs.loc[mask, "has_tcr"] = hit_flags

        print(
            f"[TCR] {sample}: matched {int(hit_flags.sum())}/{mask.sum()} GEX cells "
            f"to {adata_tcr.n_obs} TCR barcodes"
        )

        if airr_indices:
            sub = adata_tcr[airr_indices].copy()
            sub.obs_names = pd.Index(gex_hit_names, dtype=object)
            matched_airr.append(sub)

    # Samples without TCR keep no_tcr_input unless already set above
    has_input = set(airr_by_sample)
    other = ~adata.obs["sample"].astype(str).isin(has_input)
    adata.obs.loc[other, "tcr_status"] = "no_tcr_input"
    adata.obs.loc[other, "has_tcr"] = False

    if not matched_airr:
        return None

    airr = sc.concat(matched_airr, join="outer", merge="same")
    if not airr.obs_names.is_unique:
        raise ValueError(
            "AIRR obs_names are not unique after GEX matching; "
            "ensure ensure_unique_obs_names() ran before match_airr_to_gex()"
        )
    return airr


def run_scirpy(mdata: mu.MuData):
    ir.pp.index_chains(mdata)
    ir.tl.chain_qc(mdata)
    ir.pp.ir_dist(mdata)
    ir.tl.define_clonotypes(mdata, receptor_arms="all", dual_ir="primary_only")
    ir.tl.clonal_expansion(mdata)
    try:
        ir.tl.clonotype_network(mdata, min_cells=2)
    except Exception as exc:  # noqa: BLE001 — optional viz helper
        print(f"[TCR] clonotype_network skipped: {exc}")


def pull_useful_obs(adata: sc.AnnData, mdata: mu.MuData):
    """Copy key AIRR columns onto GEX obs for Scanpy plotting."""
    candidates = [
        "airr:receptor_type",
        "airr:receptor_subtype",
        "airr:chain_pairing",
        "airr:clone_id",
        "airr:clone_id_size",
        "airr:clonal_expansion",
        "receptor_type",
        "receptor_subtype",
        "chain_pairing",
        "clone_id",
        "clone_id_size",
        "clonal_expansion",
    ]
    src = mdata.obs
    for col in candidates:
        if col not in src.columns:
            continue
        dest = col.split(":", 1)[-1]
        # Align to full GEX index; unmatched cells stay NA
        values = src[col]
        if not values.index.equals(adata.obs_names):
            values = values.reindex(adata.obs_names)
        adata.obs[dest] = values

    if "chain_pairing" in adata.obs.columns:
        pairing = adata.obs["chain_pairing"]
        bad = pairing.isin(["no IR", "orphan VJ", "orphan VDJ", "ambiguous"]) | pairing.isna()
        adata.obs["has_productive_tcr"] = ~bad
    else:
        adata.obs["has_productive_tcr"] = adata.obs["has_tcr"]


def summarize(adata: sc.AnnData) -> pd.DataFrame:
    rows = []
    for sample, sub in adata.obs.groupby("sample", observed=True):
        n = len(sub)
        n_tcr = int(sub["has_tcr"].sum())
        n_prod = (
            int(sub["has_productive_tcr"].sum())
            if "has_productive_tcr" in sub.columns
            else n_tcr
        )
        n_clones = (
            int(sub.loc[sub["has_tcr"], "clone_id"].nunique())
            if "clone_id" in sub.columns and n_tcr
            else 0
        )
        rows.append(
            {
                "sample": sample,
                "patient_id": sub["patient_id"].iloc[0]
                if "patient_id" in sub.columns
                else "",
                "n_gex_cells": n,
                "n_tcr_matched": n_tcr,
                "pct_tcr_matched": round(100.0 * n_tcr / n, 2) if n else 0.0,
                "n_productive_tcr": n_prod,
                "n_clonotypes": n_clones,
                "tcr_status": ",".join(sorted(sub["tcr_status"].astype(str).unique())),
            }
        )
    return pd.DataFrame(rows)


def clonotype_table(adata: sc.AnnData) -> pd.DataFrame:
    empty_cols = ["sample", "clone_id", "n_cells", "chain_pairing", "cell_type"]
    if "clone_id" not in adata.obs.columns:
        return pd.DataFrame(columns=empty_cols)
    optional = [
        c
        for c in ("chain_pairing", "cell_type", "predicted.celltype.l2", "patient_id")
        if c in adata.obs.columns
    ]
    df = adata.obs.loc[adata.obs["has_tcr"], ["sample", "clone_id", *optional]].copy()
    if df.empty:
        return pd.DataFrame(columns=["sample", "clone_id", "n_cells", *optional])

    records = []
    for (sample, clone), sub in df.groupby(["sample", "clone_id"], observed=True):
        rec = {"sample": sample, "clone_id": clone, "n_cells": len(sub)}
        for c in optional:
            mode = sub[c].mode()
            rec[c] = mode.iloc[0] if len(mode) else sub[c].iloc[0]
        records.append(rec)
    out = pd.DataFrame.from_records(records)
    return out.sort_values(["sample", "n_cells"], ascending=[True, False])


def plot_umap_if_present(adata: sc.AnnData, color: str, filename: str):
    if "X_umap" not in adata.obsm or color not in adata.obs.columns:
        return
    sc.pl.umap(adata, color=color, show=False, frameon=False)
    plt.tight_layout()
    plt.savefig(filename, dpi=150, bbox_inches="tight")
    plt.close("all")


def plot_summary_bars(summary: pd.DataFrame, filename: str):
    if summary.empty:
        return
    fig, ax = plt.subplots(figsize=(max(6, len(summary) * 0.8), 4))
    x = np.arange(len(summary))
    ax.bar(x, summary["n_gex_cells"], color="#9aa5b1", label="GEX cells")
    ax.bar(x, summary["n_tcr_matched"], color="#2f6fed", label="TCR matched")
    ax.set_xticks(x)
    ax.set_xticklabels(summary["sample"], rotation=45, ha="right")
    ax.set_ylabel("Cells")
    ax.set_title("GEX cells with matched TCR")
    ax.legend(frameon=False)
    fig.tight_layout()
    fig.savefig(filename, dpi=150)
    plt.close(fig)


def plot_chain_pairing(adata: sc.AnnData, filename: str):
    if "chain_pairing" not in adata.obs.columns:
        return
    sub = adata.obs.loc[adata.obs["has_tcr"]]
    if sub.empty:
        return
    counts = (
        sub.groupby(["sample", "chain_pairing"], observed=True)
        .size()
        .reset_index(name="n")
    )
    fig, ax = plt.subplots(figsize=(8, 4))
    sns.barplot(data=counts, x="sample", y="n", hue="chain_pairing", ax=ax)
    ax.set_title("TCR chain pairing by sample")
    ax.set_ylabel("Cells")
    ax.legend(bbox_to_anchor=(1.02, 1), loc="upper left", frameon=False)
    fig.tight_layout()
    fig.savefig(filename, dpi=150, bbox_inches="tight")
    plt.close(fig)


def plot_clonotype_sizes(adata: sc.AnnData, filename: str):
    if "clone_id_size" not in adata.obs.columns:
        return
    sub = adata.obs.loc[adata.obs["has_tcr"] & adata.obs["clone_id_size"].notna()]
    if sub.empty:
        return
    fig, ax = plt.subplots(figsize=(6, 4))
    for sample, s in sub.groupby("sample", observed=True):
        ax.hist(s["clone_id_size"].astype(float), bins=30, alpha=0.5, label=str(sample))
    ax.set_xlabel("Clonotype size")
    ax.set_ylabel("Cells")
    ax.set_title("Clonotype size distribution")
    ax.legend(frameon=False)
    fig.tight_layout()
    fig.savefig(filename, dpi=150)
    plt.close(fig)


def write_outputs(adata, mdata, load_status):
    summary = summarize(adata)
    clones = clonotype_table(adata)
    summary.to_csv("tcr_summary.csv", index=False)
    clones.to_csv("clonotypes_by_sample.csv", index=False)
    pd.DataFrame(load_status).to_csv("tcr_input_status.csv", index=False)
    plot_summary_bars(summary, "tcr_match_counts.png")
    plot_chain_pairing(adata, "tcr_chain_pairing.png")
    plot_clonotype_sizes(adata, "tcr_clonotype_sizes.png")
    plot_umap_if_present(adata, "has_tcr", "umap_has_tcr.png")
    plot_umap_if_present(adata, "chain_pairing", "umap_chain_pairing.png")
    plot_umap_if_present(adata, "clonal_expansion", "umap_clonal_expansion.png")
    adata.write_h5ad("adata_gex_tcr.h5ad", compression="gzip")
    mdata.write("mudata_gex_tcr.h5mu")
    print(
        f"[TCR] Done. Matched {int(adata.obs['has_tcr'].sum())}/{adata.n_obs} cells"
    )


def main():
    args = parse_args()
    tcr_map = json.loads(args.tcr_map)
    if isinstance(tcr_map, str):
        tcr_map = json.loads(tcr_map)

    adata = sc.read_h5ad(args.adata)
    print(f"[TCR] GEX: {adata.n_obs} cells × {adata.n_vars} genes")
    # Merged objects reuse 10x barcodes across samples. Uniquify before scirpy
    # (chain_qc) and before MapQuery join; keep barcode_raw for TCR matching.
    ensure_unique_obs_names(adata, join=".")
    load_mapquery_labels(args.mapquery_labels, adata)

    airr_by_sample = {}
    load_status = []
    for entry in tcr_map:
        sample = entry["sample"]
        path = entry["path"]
        adata_tcr, status = load_sample_tcr(sample, path)
        load_status.append({"sample": sample, "path": path, "status": status})
        if adata_tcr is not None:
            airr_by_sample[sample] = adata_tcr

    if not airr_by_sample:
        adata.obs["has_tcr"] = False
        adata.obs["tcr_status"] = "no_tcr_input" if not tcr_map else "missing"
        adata.obs["has_productive_tcr"] = False
        print("[TCR] No usable TCR libraries; writing GEX-only stub outputs")
        write_outputs(adata, mu.MuData({"gex": adata}), load_status)
        return

    airr = match_airr_to_gex(adata, airr_by_sample)
    if airr is None:
        adata.obs["has_productive_tcr"] = False
        print("[TCR] No TCR cells matched GEX; writing unmatched outputs")
        write_outputs(adata, mu.MuData({"gex": adata}), load_status)
        return

    mdata = mu.MuData({"gex": adata, "airr": airr})
    run_scirpy(mdata)
    pull_useful_obs(adata, mdata)
    mdata.mod["gex"] = adata
    write_outputs(adata, mdata, load_status)


if __name__ == "__main__":
    main()
