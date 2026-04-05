#!/usr/bin/env python3
"""Map regulon activity to spatial transcriptomics data."""

from __future__ import annotations

import argparse
import json
import logging
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc

PROJECT_ROOT = Path(__file__).resolve().parents[2]
LOG_DIR = PROJECT_ROOT / "logs"


def setup_logging() -> logging.Logger:
    LOG_DIR.mkdir(exist_ok=True)
    logger = logging.getLogger("map_regulons")
    logger.setLevel(logging.INFO)
    if logger.handlers:
        return logger
    formatter = logging.Formatter("%(asctime)s | %(levelname)s | %(message)s")
    fh = logging.FileHandler(LOG_DIR / "map_regulons.log")
    fh.setFormatter(formatter)
    sh = logging.StreamHandler()
    sh.setFormatter(formatter)
    logger.addHandler(fh)
    logger.addHandler(sh)
    return logger


def score_regulons(adata: sc.AnnData, regulons: dict[str, list[str]]) -> pd.DataFrame:
    scores: dict[str, np.ndarray] = {}
    for tf, genes in regulons.items():
        present = [gene for gene in genes if gene in adata.var_names]
        if not present:
            continue
        matrix = adata[:, present].X
        mean_axis = matrix.mean(axis=1)
        scores[tf] = mean_axis.A1 if hasattr(mean_axis, "A1") else np.asarray(mean_axis).ravel()
    return pd.DataFrame(scores, index=adata.obs_names)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--spatial", type=Path, required=True)
    parser.add_argument("--regulons", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    logger = setup_logging()
    adata = sc.read_h5ad(args.spatial)
    regulons = json.loads(args.regulons.read_text(encoding="utf-8"))
    scores = score_regulons(adata, regulons)
    args.output_dir.mkdir(parents=True, exist_ok=True)
    scores.to_csv(args.output_dir / "regulon_activity.tsv", sep="\t")
    logger.info("Saved scores for %s regulons across %s spots", scores.shape[1], scores.shape[0])


if __name__ == "__main__":
    main()
