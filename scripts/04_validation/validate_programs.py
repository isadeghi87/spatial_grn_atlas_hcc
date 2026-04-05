#!/usr/bin/env python3
"""Validate candidate TF programs using summary network metrics and optional TCGA files."""

from __future__ import annotations

import argparse
import logging
from pathlib import Path

import pandas as pd

PROJECT_ROOT = Path(__file__).resolve().parents[2]
LOG_DIR = PROJECT_ROOT / "logs"


def setup_logging() -> logging.Logger:
    LOG_DIR.mkdir(exist_ok=True)
    logger = logging.getLogger("validate_programs")
    logger.setLevel(logging.INFO)
    if logger.handlers:
        return logger
    formatter = logging.Formatter("%(asctime)s | %(levelname)s | %(message)s")
    fh = logging.FileHandler(LOG_DIR / "validate_programs.log")
    fh.setFormatter(formatter)
    sh = logging.StreamHandler()
    sh.setFormatter(formatter)
    logger.addHandler(fh)
    logger.addHandler(sh)
    return logger


def summarize_edges(edges: pd.DataFrame) -> pd.DataFrame:
    summary = (
        edges.groupby("tf")
        .agg(targets=("target", "nunique"), mean_score=("score", "mean"), max_score=("score", "max"))
        .sort_values(["targets", "mean_score", "tf"], ascending=[False, False, True])
        .reset_index()
    )
    summary["priority_rank"] = (
        summary[["targets", "mean_score", "max_score"]]
        .apply(tuple, axis=1)
        .rank(method="min", ascending=False)
        .astype(int)
    )
    return summary


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--grn", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    logger = setup_logging()
    edges = pd.read_csv(args.grn, sep="\t")
    summary = summarize_edges(edges)
    args.output_dir.mkdir(parents=True, exist_ok=True)
    summary.to_csv(args.output_dir / "tf_priority_summary.tsv", sep="\t", index=False)
    logger.info("Saved %s prioritized TF summaries", len(summary))


if __name__ == "__main__":
    main()
