#!/usr/bin/env python3
"""Preprocess public HCC scRNA-seq cohorts into one AnnData object."""

from __future__ import annotations

import argparse
import logging
from pathlib import Path
from typing import Iterable

import anndata as ad
import numpy as np
import scanpy as sc

PROJECT_ROOT = Path(__file__).resolve().parents[2]
LOG_DIR = PROJECT_ROOT / "logs"


def setup_logging() -> logging.Logger:
    LOG_DIR.mkdir(exist_ok=True)
    logger = logging.getLogger("preprocess_scrna")
    logger.setLevel(logging.INFO)
    if logger.handlers:
        return logger
    formatter = logging.Formatter("%(asctime)s | %(levelname)s | %(message)s")
    fh = logging.FileHandler(LOG_DIR / "preprocess_scrna.log")
    fh.setFormatter(formatter)
    sh = logging.StreamHandler()
    sh.setFormatter(formatter)
    logger.addHandler(fh)
    logger.addHandler(sh)
    return logger


def load_h5ad_files(input_dir: Path, logger: logging.Logger) -> list[ad.AnnData]:
    files = sorted(input_dir.rglob("*.h5ad"))
    datasets: list[ad.AnnData] = []
    for path in files:
        logger.info("Loading %s", path)
        adata = sc.read_h5ad(path)
        adata.obs["cohort"] = path.parent.name
        datasets.append(adata)
    return datasets


def qc_filter(adata: ad.AnnData) -> ad.AnnData:
    adata.var["mt"] = adata.var_names.str.upper().str.startswith("MT-")
    sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], inplace=True)
    return adata[
        (adata.obs["n_genes_by_counts"] >= 200)
        & (adata.obs["n_genes_by_counts"] <= 6000)
        & (adata.obs["pct_counts_mt"] <= 25)
    ].copy()


def annotate_simple_celltypes(adata: ad.AnnData) -> None:
    marker_map = {
        "T_cell": ["CD3D", "CD3E", "TRBC1"],
        "B_cell": ["MS4A1", "CD79A"],
        "Myeloid": ["LYZ", "S100A8", "CTSS"],
        "Endothelial": ["KDR", "PECAM1", "VWF"],
        "Fibroblast": ["COL1A1", "COL3A1", "DCN"],
        "Malignant_like": ["EPCAM", "KRT8", "KRT18"],
    }
    scores = {}
    for label, genes in marker_map.items():
        present = [gene for gene in genes if gene in adata.var_names]
        if present:
            matrix = adata[:, present].X
            scores[label] = np.asarray(matrix.mean(axis=1)).ravel()
        else:
            scores[label] = np.zeros(adata.n_obs)
    score_frame = np.column_stack([scores[label] for label in marker_map])
    labels = list(marker_map)
    best_idx = np.argmax(score_frame, axis=1)
    best_score = score_frame[np.arange(score_frame.shape[0]), best_idx]
    adata.obs["cell_type"] = [
        labels[idx] if score > 0 else "Unassigned" for idx, score in zip(best_idx, best_score, strict=False)
    ]
    adata.obs["cell_type_score"] = best_score


def preprocess(datasets: Iterable[ad.AnnData], logger: logging.Logger) -> ad.AnnData:
    filtered = [qc_filter(dataset) for dataset in datasets]
    union_gene_count = len(set().union(*(set(dataset.var_names) for dataset in filtered)))
    merged = ad.concat(filtered, join="inner", label="batch")
    merged.var_names_make_unique()
    logger.info("Retained %s shared genes from a union of %s genes", merged.n_vars, union_gene_count)
    sc.pp.normalize_total(merged, target_sum=1e4)
    sc.pp.log1p(merged)
    sc.pp.highly_variable_genes(merged, n_top_genes=3000, subset=True, flavor="seurat")
    sc.pp.scale(merged, max_value=10)
    sc.tl.pca(merged, svd_solver="arpack")
    sc.pp.neighbors(merged)
    sc.tl.umap(merged)
    sc.tl.leiden(merged, key_added="cluster")
    sc.tl.rank_genes_groups(merged, "cluster", method="wilcoxon")
    annotate_simple_celltypes(merged)
    logger.info("Preprocessed %s cells and %s genes", merged.n_obs, merged.n_vars)
    return merged


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-dir", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    logger = setup_logging()
    datasets = load_h5ad_files(args.input_dir, logger)
    if not datasets:
        raise FileNotFoundError(f"No h5ad files found under {args.input_dir}")
    processed = preprocess(datasets, logger)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    processed.write(args.output)
    logger.info("Saved %s", args.output)


if __name__ == "__main__":
    main()
