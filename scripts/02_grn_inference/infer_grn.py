#!/usr/bin/env python3
"""Infer compartment-specific GRNs from processed scRNA-seq data."""

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
DEFAULT_TF_LIST = PROJECT_ROOT / "data" / "references" / "human_tf_census.txt"


def setup_logging() -> logging.Logger:
    LOG_DIR.mkdir(exist_ok=True)
    logger = logging.getLogger("infer_grn")
    logger.setLevel(logging.INFO)
    if logger.handlers:
        return logger
    formatter = logging.Formatter("%(asctime)s | %(levelname)s | %(message)s")
    fh = logging.FileHandler(LOG_DIR / "infer_grn.log")
    fh.setFormatter(formatter)
    sh = logging.StreamHandler()
    sh.setFormatter(formatter)
    logger.addHandler(fh)
    logger.addHandler(sh)
    return logger


def choose_tfs(var_names: pd.Index, tf_list_path: Path | None) -> list[str]:
    if tf_list_path and tf_list_path.exists():
        tf_list = {
            line.strip().split("\t")[0]
            for line in tf_list_path.read_text(encoding="utf-8").splitlines()
            if line.strip() and not line.startswith("#")
        }
        return [gene for gene in var_names if gene in tf_list]
    candidates = [
        gene
        for gene in var_names
        if gene.isupper() and len(gene) <= 10 and not gene.startswith(("MT-", "RPL", "RPS"))
    ]
    return candidates[:200]


def expression_frame(adata: sc.AnnData) -> pd.DataFrame:
    matrix = adata.X.toarray() if hasattr(adata.X, "toarray") else np.asarray(adata.X)
    return pd.DataFrame(matrix, index=adata.obs_names, columns=adata.var_names)


def build_edges(
    frame: pd.DataFrame,
    tfs: list[str],
    min_abs_corr: float,
    max_targets_per_tf: int,
) -> pd.DataFrame:
    corr = frame.corr(method="spearman")
    rows = []
    for tf in tfs:
        if tf not in corr.index:
            continue
        ranked = corr.loc[tf].abs().drop(labels=[tf], errors="ignore").sort_values(ascending=False)
        for target, score in ranked[ranked >= min_abs_corr].head(max_targets_per_tf).items():
            rows.append({"tf": tf, "target": target, "score": float(score)})
    return pd.DataFrame(rows).drop_duplicates(["tf", "target"])


def export_regulons(edges: pd.DataFrame, output_dir: Path) -> None:
    regulons = {
        tf: sorted(group["target"].tolist())
        for tf, group in edges.sort_values("score", ascending=False).groupby("tf")
    }
    (output_dir / "regulons.json").write_text(json.dumps(regulons, indent=2), encoding="utf-8")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--rna", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--min-abs-corr", type=float, default=0.3)
    parser.add_argument("--max-targets-per-tf", type=int, default=100)
    parser.add_argument("--tf-list", type=Path, default=DEFAULT_TF_LIST)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    logger = setup_logging()
    adata = sc.read_h5ad(args.rna)
    frame = expression_frame(adata)
    tfs = choose_tfs(frame.columns, args.tf_list)
    logger.info("Using %s candidate TFs", len(tfs))
    edges = build_edges(frame, tfs, args.min_abs_corr, args.max_targets_per_tf)
    args.output_dir.mkdir(parents=True, exist_ok=True)
    edges.to_csv(args.output_dir / "grn_edges.tsv", sep="\t", index=False)
    export_regulons(edges, args.output_dir)
    logger.info("Saved %s edges", len(edges))


if __name__ == "__main__":
    main()
