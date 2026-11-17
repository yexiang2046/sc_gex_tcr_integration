#!/usr/bin/env python3
"""Smoke test for bin/tcr_integrate.py using synthetic 10x V(D)J + GEX fixtures."""

from __future__ import annotations

import json
import subprocess
import sys
import tempfile
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "bin" / "tcr_integrate.py"


def _write_fixtures(tmpdir: Path):
    rng = np.random.default_rng(0)
    barcodes_s1 = [f"AAAAAAAAAAAAAAA{i}-1" for i in range(6)]
    barcodes_s2 = [f"BBBBBBBBBBBBBBB{i}-1" for i in range(8)]
    obs = pd.DataFrame(
        {
            "sample": ["sample1"] * 6 + ["sample2"] * 8,
            "patient_id": ["HD31"] * 6 + ["HD39"] * 8,
            "cell_type": ["CD8 TEM"] * 14,
        },
        index=barcodes_s1 + barcodes_s2,
    )
    adata = ad.AnnData(
        X=rng.poisson(1.0, size=(14, 20)).astype(np.float32),
        obs=obs,
        var=pd.DataFrame(index=[f"G{i}" for i in range(20)]),
    )
    adata.obsm["X_umap"] = rng.normal(size=(14, 2))
    adata_path = tmpdir / "adata_annotated.h5ad"
    adata.write_h5ad(adata_path)

    mq = tmpdir / "mapquery_labels.csv"
    pd.DataFrame(
        {
            "barcode": adata.obs_names,
            "predicted.celltype.l2": "CD8 TEM",
            "predicted.celltype.l2.score": 0.9,
        }
    ).to_csv(mq, index=False)

    tcr_dir = tmpdir / "tcr_sample2"
    tcr_dir.mkdir()
    rows = []
    for i in range(5):
        bc = barcodes_s2[i]
        clone = f"clonotype{i + 1}"
        rows.append(
            {
                "barcode": bc,
                "is_cell": "true",
                "contig_id": f"{bc}_contig_1",
                "high_confidence": "true",
                "length": 500,
                "chain": "TRA",
                "v_gene": "TRAV12-2",
                "d_gene": "",
                "j_gene": "TRAJ9",
                "c_gene": "TRAC",
                "full_length": "true",
                "productive": "true",
                "cdr3": f"CAV{i}NYQLIW",
                "cdr3_nt": f"TGTGCC{'A' * (i + 1)}GTAT",
                "reads": 100 + i,
                "umis": 5 + i,
                "raw_clonotype_id": clone,
                "raw_consensus_id": f"{clone}_consensus_1",
                "exact_subclonotype_id": 1,
            }
        )
        rows.append(
            {
                "barcode": bc,
                "is_cell": "true",
                "contig_id": f"{bc}_contig_2",
                "high_confidence": "true",
                "length": 600,
                "chain": "TRB",
                "v_gene": "TRBV20-1",
                "d_gene": "TRBD1",
                "j_gene": "TRBJ2-7",
                "c_gene": "TRBC2",
                "full_length": "true",
                "productive": "true",
                "cdr3": f"CASS{i}YEQYF",
                "cdr3_nt": f"TGTGCC{'T' * (i + 1)}GAG",
                "reads": 200 + i,
                "umis": 8 + i,
                "raw_clonotype_id": clone,
                "raw_consensus_id": f"{clone}_consensus_2",
                "exact_subclonotype_id": 1,
            }
        )
    pd.DataFrame(rows).to_csv(tcr_dir / "filtered_contig_annotations.csv", index=False)
    return adata_path, mq, tcr_dir


def test_tcr_integrate_matches_by_sample_and_barcode():
    with tempfile.TemporaryDirectory() as td:
        tmpdir = Path(td)
        out = tmpdir / "out"
        out.mkdir()
        adata_path, mq, tcr_dir = _write_fixtures(tmpdir)
        tcr_map = json.dumps(
            [
                {"sample": "sample2", "path": str(tcr_dir)},
                {"sample": "sample3", "path": str(tmpdir / "missing_tcr")},
            ]
        )
        proc = subprocess.run(
            [
                sys.executable,
                str(SCRIPT),
                "--adata",
                str(adata_path),
                "--mapquery_labels",
                str(mq),
                "--tcr_map",
                tcr_map,
                "--seed",
                "42",
            ],
            cwd=out,
            check=True,
            capture_output=True,
            text=True,
        )
        assert "Matched 5/" in proc.stdout

        summary = pd.read_csv(out / "tcr_summary.csv")
        s2 = summary.set_index("sample").loc["sample2"]
        assert int(s2["n_tcr_matched"]) == 5
        assert int(summary.set_index("sample").loc["sample1"]["n_tcr_matched"]) == 0

        status = pd.read_csv(out / "tcr_input_status.csv").set_index("sample")
        assert status.loc["sample2", "status"] == "filtered"
        assert status.loc["sample3", "status"] == "missing"

        adata = sc.read_h5ad(out / "adata_gex_tcr.h5ad")
        assert int(adata.obs["has_tcr"].sum()) == 5
        assert "chain_pairing" in adata.obs.columns
        assert (out / "mudata_gex_tcr.h5mu").exists()


def test_tcr_integrate_empty_map():
    with tempfile.TemporaryDirectory() as td:
        tmpdir = Path(td)
        out = tmpdir / "out"
        out.mkdir()
        adata_path, mq, _ = _write_fixtures(tmpdir)
        proc = subprocess.run(
            [
                sys.executable,
                str(SCRIPT),
                "--adata",
                str(adata_path),
                "--mapquery_labels",
                str(mq),
                "--tcr_map",
                "[]",
                "--seed",
                "42",
            ],
            cwd=out,
            check=True,
            capture_output=True,
            text=True,
        )
        assert "No usable TCR" in proc.stdout or "No TCR paths" in proc.stdout
        adata = sc.read_h5ad(out / "adata_gex_tcr.h5ad")
        assert int(adata.obs["has_tcr"].sum()) == 0


def _contig_rows(barcodes):
    rows = []
    for i, bc in enumerate(barcodes):
        clone = f"clonotype{i + 1}"
        for chain, v, j, c, cdr3, nt, reads, umis in (
            ("TRA", "TRAV12-2", "TRAJ9", "TRAC", f"CAV{i}NYQLIW", f"TGT{'A'*(i+1)}", 100 + i, 5 + i),
            ("TRB", "TRBV20-1", "TRBJ2-7", "TRBC2", f"CASS{i}YEQYF", f"TGT{'T'*(i+1)}", 200 + i, 8 + i),
        ):
            rows.append(
                {
                    "barcode": bc,
                    "is_cell": "true",
                    "contig_id": f"{bc}_contig_{chain}",
                    "high_confidence": "true",
                    "length": 500,
                    "chain": chain,
                    "v_gene": v,
                    "d_gene": "TRBD1" if chain == "TRB" else "",
                    "j_gene": j,
                    "c_gene": c,
                    "full_length": "true",
                    "productive": "true",
                    "cdr3": cdr3,
                    "cdr3_nt": nt,
                    "reads": reads,
                    "umis": umis,
                    "raw_clonotype_id": clone,
                    "raw_consensus_id": f"{clone}_consensus_{chain}",
                    "exact_subclonotype_id": 1,
                }
            )
    return rows


def test_tcr_integrate_uniquifies_duplicate_gex_barcodes():
    """Merged HT objects reuse 10x barcodes across samples; chain_qc needs unique names."""
    with tempfile.TemporaryDirectory() as td:
        tmpdir = Path(td)
        out = tmpdir / "out"
        out.mkdir()
        rng = np.random.default_rng(1)
        # Same barcode strings in both samples (the production failure mode).
        shared = [f"SHARED{i:010d}-1" for i in range(4)]
        obs = pd.DataFrame(
            {
                "sample": ["sample1"] * 4 + ["sample2"] * 4,
                "patient_id": ["HD31"] * 4 + ["HD39"] * 4,
                "cell_type": ["CD8 TEM"] * 8,
            },
            index=shared + shared,
        )
        assert not obs.index.is_unique
        adata = ad.AnnData(
            X=rng.poisson(1.0, size=(8, 10)).astype(np.float32),
            obs=obs,
            var=pd.DataFrame(index=[f"G{i}" for i in range(10)]),
        )
        adata.obsm["X_umap"] = rng.normal(size=(8, 2))
        adata_path = tmpdir / "adata_dup.h5ad"
        adata.write_h5ad(adata_path)

        # MapQuery-style labels already uniquified with join='.'
        mq_names = list(shared) + [f"{b}.1" for b in shared]
        mq = tmpdir / "mapquery_labels.csv"
        pd.DataFrame(
            {
                "barcode": mq_names,
                "predicted.celltype.l2": "CD8 TEM",
                "predicted.celltype.l2.score": 0.9,
            }
        ).to_csv(mq, index=False)

        tcr1 = tmpdir / "tcr1"
        tcr2 = tmpdir / "tcr2"
        tcr1.mkdir()
        tcr2.mkdir()
        pd.DataFrame(_contig_rows(shared[:3])).to_csv(
            tcr1 / "filtered_contig_annotations.csv", index=False
        )
        pd.DataFrame(_contig_rows(shared[:3])).to_csv(
            tcr2 / "filtered_contig_annotations.csv", index=False
        )

        tcr_map = json.dumps(
            [
                {"sample": "sample1", "path": str(tcr1)},
                {"sample": "sample2", "path": str(tcr2)},
            ]
        )
        proc = subprocess.run(
            [
                sys.executable,
                str(SCRIPT),
                "--adata",
                str(adata_path),
                "--mapquery_labels",
                str(mq),
                "--tcr_map",
                tcr_map,
                "--seed",
                "42",
            ],
            cwd=out,
            check=True,
            capture_output=True,
            text=True,
        )
        assert "Uniquified GEX obs_names" in proc.stdout
        assert "chain_qc" not in proc.stderr.lower() or proc.returncode == 0

        result = sc.read_h5ad(out / "adata_gex_tcr.h5ad")
        assert result.obs_names.is_unique
        assert int(result.obs["has_tcr"].sum()) == 6  # 3 per sample
        assert "chain_pairing" in result.obs.columns
        summary = pd.read_csv(out / "tcr_summary.csv").set_index("sample")
        assert int(summary.loc["sample1", "n_tcr_matched"]) == 3
        assert int(summary.loc["sample2", "n_tcr_matched"]) == 3


if __name__ == "__main__":
    test_tcr_integrate_matches_by_sample_and_barcode()
    test_tcr_integrate_empty_map()
    test_tcr_integrate_uniquifies_duplicate_gex_barcodes()
    print("OK")
