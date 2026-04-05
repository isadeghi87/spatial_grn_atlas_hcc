#!/usr/bin/env python3
"""Preprocess public HCC spatial transcriptomics data."""

from __future__ import annotations

import argparse
import logging
from pathlib import Path

import scanpy as sc

PROJECT_ROOT = Path(__file__).resolve().parents[2]
LOG_DIR = PROJECT_ROOT / "logs"


def setup_logging() -> logging.Logger:
    LOG_DIR.mkdir(exist_ok=True)
    logger = logging.getLogger("preprocess_spatial")
    logger.setLevel(logging.INFO)
    if logger.handlers:
        return logger
    formatter = logging.Formatter("%(asctime)s | %(levelname)s | %(message)s")
    fh = logging.FileHandler(LOG_DIR / "preprocess_spatial.log")
    fh.setFormatter(formatter)
    sh = logging.StreamHandler()
    sh.setFormatter(formatter)
    logger.addHandler(fh)
    logger.addHandler(sh)
    return logger


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-dir", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    logger = setup_logging()
    files = sorted(args.input_dir.rglob("*.h5ad"))
    if not files:
        logger.warning("No spatial h5ad files found under %s; skipping write.", args.input_dir)
        return
    adata = sc.read_h5ad(files[0])
    sc.pp.filter_genes(adata, min_cells=3)
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata, n_top_genes=2000, subset=True, flavor="seurat")
    sc.tl.pca(adata, svd_solver="arpack")
    sc.pp.neighbors(adata)
    sc.tl.umap(adata)
    sc.tl.leiden(adata, key_added="domain")
    args.output.parent.mkdir(parents=True, exist_ok=True)
    adata.write(args.output)
    logger.info("Saved %s", args.output)


if __name__ == "__main__":
    main()
