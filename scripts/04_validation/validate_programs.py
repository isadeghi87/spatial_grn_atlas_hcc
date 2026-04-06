#!/usr/bin/env python3
"""Validate candidate TF programs using summary network metrics and optional TCGA files."""

from __future__ import annotations

import argparse
import logging
import sys
from pathlib import Path

import pandas as pd

PROJECT_ROOT = Path(__file__).resolve().parents[2]
LOG_DIR = PROJECT_ROOT / "logs"

if str(PROJECT_ROOT) not in sys.path:
    sys.path.append(str(PROJECT_ROOT))

from scripts.shared_analysis import get_config_value, load_analysis_config


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
    if edges.empty:
        return pd.DataFrame(columns=["tf", "targets", "mean_score", "max_score", "priority_rank"])
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


def summarize_label_specificity(edges: pd.DataFrame, label_column: str) -> pd.DataFrame:
    if edges.empty or label_column not in edges.columns:
        return pd.DataFrame(columns=["tf", f"{label_column}_support", f"dominant_{label_column}", f"{label_column}_specificity"])
    grouped = (
        edges.groupby(["tf", label_column])
        .agg(targets=("target", "nunique"), mean_score=("score", "mean"))
        .reset_index()
    )
    grouped["weighted_support"] = grouped["targets"] * grouped["mean_score"]
    dominant = (
        grouped.sort_values(["tf", "weighted_support", label_column], ascending=[True, False, True])
        .groupby("tf", as_index=False)
        .first()
    )
    total_support = grouped.groupby("tf")["weighted_support"].sum().rename("total_weighted_support")
    dominant = dominant.merge(total_support, on="tf", how="left")
    dominant[f"{label_column}_specificity"] = dominant["weighted_support"] / dominant["total_weighted_support"].where(
        dominant["total_weighted_support"] != 0, 1
    )
    return dominant.rename(
        columns={
            label_column: f"dominant_{label_column}",
        }
    )[["tf", f"dominant_{label_column}", "weighted_support", f"{label_column}_specificity"]].merge(
        grouped.groupby("tf")[label_column].nunique().rename(f"{label_column}_support").reset_index(),
        on="tf",
        how="left",
    )


def summarize_spatial_support(spatial_stats: pd.DataFrame, niche_calls: pd.DataFrame) -> pd.DataFrame:
    if spatial_stats.empty:
        return pd.DataFrame(columns=["tf", "morans_i", "moran_p_value", "dominant_niche", "dominant_niche_fraction"])
    if niche_calls.empty or "dominant_regulon" not in niche_calls.columns:
        return spatial_stats[["tf", "morans_i", "moran_p_value"]].copy()
    niche_counts = (
        niche_calls.groupby(["dominant_regulon", "niche_label"])
        .size()
        .rename("spots")
        .reset_index()
        .sort_values(["dominant_regulon", "spots", "niche_label"], ascending=[True, False, True])
    )
    dominant_niche = niche_counts.groupby("dominant_regulon", as_index=False).first()
    total_spots = niche_counts.groupby("dominant_regulon")["spots"].sum().rename("total_spots")
    dominant_niche = dominant_niche.merge(total_spots, on="dominant_regulon", how="left")
    dominant_niche["dominant_niche_fraction"] = dominant_niche["spots"] / dominant_niche["total_spots"].where(
        dominant_niche["total_spots"] != 0, 1
    )
    dominant_niche = dominant_niche.rename(
        columns={
            "dominant_regulon": "tf",
            "niche_label": "dominant_niche",
        }
    )[["tf", "dominant_niche", "dominant_niche_fraction"]]
    return spatial_stats.merge(dominant_niche, on="tf", how="left")


def build_advanced_summary(
    summary: pd.DataFrame,
    compartment_summary: pd.DataFrame,
    cohort_summary: pd.DataFrame,
    spatial_summary: pd.DataFrame,
) -> pd.DataFrame:
    advanced = summary.copy()
    for other in [compartment_summary, cohort_summary, spatial_summary]:
        if not other.empty:
            advanced = advanced.merge(other, on="tf", how="left")
    return advanced.sort_values(["priority_rank", "max_score", "tf"], ascending=[True, False, True]).reset_index(
        drop=True
    )


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--grn", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--config", type=Path, default=PROJECT_ROOT / "config" / "analysis_config.yaml")
    parser.add_argument("--compartment-grn", type=Path)
    parser.add_argument("--cohort-grn", type=Path)
    parser.add_argument("--spatial-stats", type=Path)
    parser.add_argument("--niche-calls", type=Path)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    logger = setup_logging()
    config = load_analysis_config(args.config)
    compartment_column = str(get_config_value(config, "analysis", "grn", "groupby", default="cell_type"))
    edges = pd.read_csv(args.grn, sep="\t")
    summary = summarize_edges(edges)
    compartment_edges = (
        pd.read_csv(args.compartment_grn, sep="\t")
        if args.compartment_grn and args.compartment_grn.exists()
        else pd.DataFrame()
    )
    cohort_edges = (
        pd.read_csv(args.cohort_grn, sep="\t") if args.cohort_grn and args.cohort_grn.exists() else pd.DataFrame()
    )
    spatial_stats = (
        pd.read_csv(args.spatial_stats, sep="\t") if args.spatial_stats and args.spatial_stats.exists() else pd.DataFrame()
    )
    niche_calls = (
        pd.read_csv(args.niche_calls, sep="\t") if args.niche_calls and args.niche_calls.exists() else pd.DataFrame()
    )
    compartment_summary = summarize_label_specificity(compartment_edges, compartment_column)
    cohort_summary = summarize_label_specificity(cohort_edges, "cohort")
    spatial_summary = summarize_spatial_support(spatial_stats, niche_calls)
    advanced_summary = build_advanced_summary(summary, compartment_summary, cohort_summary, spatial_summary)
    args.output_dir.mkdir(parents=True, exist_ok=True)
    summary.to_csv(args.output_dir / "tf_priority_summary.tsv", sep="\t", index=False)
    if not compartment_summary.empty:
        compartment_summary.to_csv(args.output_dir / "tf_compartment_specificity.tsv", sep="\t", index=False)
    if not cohort_summary.empty:
        cohort_summary.to_csv(args.output_dir / "tf_cohort_support.tsv", sep="\t", index=False)
    if not spatial_summary.empty:
        spatial_summary.to_csv(args.output_dir / "tf_spatial_support.tsv", sep="\t", index=False)
    if not advanced_summary.empty:
        advanced_summary.to_csv(args.output_dir / "advanced_validation_summary.tsv", sep="\t", index=False)
    logger.info("Saved %s prioritized TF summaries", len(summary))


if __name__ == "__main__":
    main()
