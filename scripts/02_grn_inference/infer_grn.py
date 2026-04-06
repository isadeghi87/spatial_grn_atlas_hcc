#!/usr/bin/env python3
"""Infer compartment-specific GRNs from processed scRNA-seq data."""

from __future__ import annotations

import argparse
import json
import logging
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc

PROJECT_ROOT = Path(__file__).resolve().parents[2]
LOG_DIR = PROJECT_ROOT / "logs"
DEFAULT_TF_LIST = PROJECT_ROOT / "data" / "references" / "human_tf_census.txt"
FALLBACK_TF_MAX_NAME_LENGTH = 10
FALLBACK_TF_CANDIDATE_LIMIT = 200

if str(PROJECT_ROOT) not in sys.path:
    sys.path.append(str(PROJECT_ROOT))

from scripts.shared_analysis import expression_frame, get_config_value, load_analysis_config


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
        if gene.isupper()
        and len(gene) <= FALLBACK_TF_MAX_NAME_LENGTH
        and not gene.startswith(("MT-", "RPL", "RPS"))
    ]
    return candidates[:FALLBACK_TF_CANDIDATE_LIMIT]


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
    edge_frame = pd.DataFrame(rows)
    if edge_frame.empty:
        return pd.DataFrame(columns=["tf", "target", "score"])
    return edge_frame.drop_duplicates(["tf", "target"])


def export_regulons(edges: pd.DataFrame, output_dir: Path) -> None:
    regulons = {
        tf: sorted(group["target"].tolist())
        for tf, group in edges.sort_values("score", ascending=False).groupby("tf")
    }
    (output_dir / "regulons.json").write_text(json.dumps(regulons, indent=2), encoding="utf-8")


def export_group_regulons(edges: pd.DataFrame, label_column: str, output_path: Path) -> None:
    grouped_regulons: dict[str, dict[str, list[str]]] = {}
    for label, label_edges in edges.groupby(label_column):
        grouped_regulons[str(label)] = {
            tf: sorted(group["target"].tolist())
            for tf, group in label_edges.sort_values("score", ascending=False).groupby("tf")
        }
    output_path.write_text(json.dumps(grouped_regulons, indent=2), encoding="utf-8")


def summarize_edges(edges: pd.DataFrame) -> pd.DataFrame:
    if edges.empty:
        return pd.DataFrame(columns=["tf", "targets", "mean_score", "max_score", "edge_count"])
    return (
        edges.groupby("tf")
        .agg(
            targets=("target", "nunique"),
            mean_score=("score", "mean"),
            max_score=("score", "max"),
            edge_count=("target", "size"),
        )
        .sort_values(["targets", "mean_score", "tf"], ascending=[False, False, True])
        .reset_index()
    )


def build_labeled_edges(
    adata: sc.AnnData,
    label_column: str,
    min_cells: int,
    min_abs_corr: float,
    max_targets_per_tf: int,
    tf_list_path: Path | None,
    logger: logging.Logger,
) -> pd.DataFrame:
    if label_column not in adata.obs:
        logger.info("Skipping %s-specific GRNs because the column is missing.", label_column)
        return pd.DataFrame(columns=[label_column, "tf", "target", "score"])

    labeled_edges: list[pd.DataFrame] = []
    for label, obs_names in adata.obs.groupby(label_column).groups.items():
        subset = adata[obs_names].copy()
        if subset.n_obs < min_cells:
            logger.info("Skipping %s=%s because it has only %s cells.", label_column, label, subset.n_obs)
            continue
        frame = expression_frame(subset)
        tfs = choose_tfs(frame.columns, tf_list_path)
        edges = build_edges(frame, tfs, min_abs_corr, max_targets_per_tf)
        if edges.empty:
            continue
        edges.insert(0, label_column, label)
        labeled_edges.append(edges)
        logger.info("Built %s %s-specific edges for %s.", len(edges), label_column, label)

    if not labeled_edges:
        return pd.DataFrame(columns=[label_column, "tf", "target", "score"])
    return pd.concat(labeled_edges, ignore_index=True)


def summarize_label_support(edges: pd.DataFrame, label_column: str) -> pd.DataFrame:
    if edges.empty:
        return pd.DataFrame(columns=["tf", f"{label_column}_support", f"dominant_{label_column}"])

    grouped = (
        edges.groupby(["tf", label_column])
        .agg(targets=("target", "nunique"), mean_score=("score", "mean"))
        .reset_index()
    )
    grouped["weighted_support"] = grouped["targets"] * grouped["mean_score"]
    dominant = (
        grouped.sort_values(["tf", "weighted_support", label_column], ascending=[True, False, True])
        .groupby("tf", as_index=False)
        .first()[["tf", label_column]]
        .rename(columns={label_column: f"dominant_{label_column}"})
    )
    support = (
        grouped.groupby("tf")
        .agg(
            **{
                f"{label_column}_support": (label_column, "nunique"),
                f"mean_{label_column}_targets": ("targets", "mean"),
            }
        )
        .reset_index()
    )
    return support.merge(dominant, on="tf", how="left")


def build_advanced_summary(
    global_edges: pd.DataFrame,
    group_edges: pd.DataFrame,
    cohort_edges: pd.DataFrame,
    group_column: str,
) -> pd.DataFrame:
    summary = summarize_edges(global_edges)
    if summary.empty:
        return summary

    if not group_edges.empty:
        summary = summary.merge(summarize_label_support(group_edges, group_column), on="tf", how="left")
    if not cohort_edges.empty:
        summary = summary.merge(summarize_label_support(cohort_edges, "cohort"), on="tf", how="left")

    fill_columns = [column for column in summary.columns if column.endswith("_support") or column.startswith("mean_")]
    for column in fill_columns:
        summary[column] = summary[column].fillna(0)

    # Favor TFs with stronger edges, broader target programs, and support across discovered cell groups.
    summary["priority_score"] = (
        summary["mean_score"] * np.log1p(summary["targets"]) * (1 + summary.get(f"{group_column}_support", 0))
    )
    return summary.sort_values(["priority_score", "max_score", "tf"], ascending=[False, False, True]).reset_index(
        drop=True
    )


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--rna", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--config", type=Path, default=PROJECT_ROOT / "config" / "analysis_config.yaml")
    parser.add_argument("--min-abs-corr", type=float)
    parser.add_argument("--max-targets-per-tf", type=int)
    parser.add_argument("--tf-list", type=Path, default=DEFAULT_TF_LIST)
    parser.add_argument("--groupby", type=str)
    parser.add_argument("--cohort-column", type=str)
    parser.add_argument("--min-group-cells", type=int)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    logger = setup_logging()
    config = load_analysis_config(args.config)
    min_abs_corr = float(get_config_value(config, "analysis", "grn", "min_abs_correlation", default=0.3))
    max_targets_per_tf = int(get_config_value(config, "analysis", "grn", "max_targets_per_tf", default=100))
    group_column = str(get_config_value(config, "analysis", "grn", "groupby", default="cell_type"))
    cohort_column = str(get_config_value(config, "analysis", "grn", "cohort_column", default="cohort"))
    min_group_cells = int(get_config_value(config, "analysis", "grn", "minimum_group_cells", default=50))

    if args.min_abs_corr is not None:
        min_abs_corr = args.min_abs_corr
    if args.max_targets_per_tf is not None:
        max_targets_per_tf = args.max_targets_per_tf
    if args.groupby:
        group_column = args.groupby
    if args.cohort_column:
        cohort_column = args.cohort_column
    if args.min_group_cells is not None:
        min_group_cells = args.min_group_cells

    adata = sc.read_h5ad(args.rna)
    frame = expression_frame(adata)
    logger.info("Computing a full gene-gene correlation matrix; this is memory intensive for very large inputs.")
    tfs = choose_tfs(frame.columns, args.tf_list)
    logger.info("Using %s candidate TFs", len(tfs))
    edges = build_edges(frame, tfs, min_abs_corr, max_targets_per_tf)
    group_edges = build_labeled_edges(
        adata,
        group_column,
        min_group_cells,
        min_abs_corr,
        max_targets_per_tf,
        args.tf_list,
        logger,
    )
    cohort_edges = build_labeled_edges(
        adata,
        cohort_column,
        min_group_cells,
        min_abs_corr,
        max_targets_per_tf,
        args.tf_list,
        logger,
    )
    advanced_summary = build_advanced_summary(edges, group_edges, cohort_edges, group_column)
    args.output_dir.mkdir(parents=True, exist_ok=True)
    edges.to_csv(args.output_dir / "grn_edges.tsv", sep="\t", index=False)
    export_regulons(edges, args.output_dir)
    summarize_edges(edges).to_csv(args.output_dir / "tf_priority_summary.tsv", sep="\t", index=False)
    if not group_edges.empty:
        group_edges.to_csv(args.output_dir / f"{group_column}_grn_edges.tsv", sep="\t", index=False)
        export_group_regulons(group_edges, group_column, args.output_dir / f"{group_column}_regulons.json")
    if not cohort_edges.empty:
        cohort_edges.to_csv(args.output_dir / f"{cohort_column}_grn_edges.tsv", sep="\t", index=False)
        export_group_regulons(cohort_edges, cohort_column, args.output_dir / f"{cohort_column}_regulons.json")
    if not advanced_summary.empty:
        advanced_summary.to_csv(args.output_dir / "advanced_tf_support.tsv", sep="\t", index=False)
    logger.info("Saved %s edges", len(edges))


if __name__ == "__main__":
    main()
