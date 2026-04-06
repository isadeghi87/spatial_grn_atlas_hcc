#!/usr/bin/env python3
"""Map regulon activity to spatial transcriptomics data."""

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
DEFAULT_NICHE_MARKERS = {
    "malignant": ["EPCAM", "KRT8", "KRT18", "AFP"],
    "immune": ["PTPRC", "CD3D", "MS4A1", "LYZ"],
    "stromal": ["COL1A1", "COL3A1", "DCN", "ACTA2"],
    "vascular": ["KDR", "PECAM1", "VWF", "EMCN"],
}

if str(PROJECT_ROOT) not in sys.path:
    sys.path.append(str(PROJECT_ROOT))

from scripts.shared_analysis import get_config_value, load_analysis_config, mean_expression, standardize


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


def prioritize_regulons(
    regulons: dict[str, list[str]],
    summary_path: Path | None,
    top_n: int,
) -> dict[str, list[str]]:
    if not regulons:
        return {}
    ranked_tfs = list(regulons)
    if summary_path and summary_path.exists():
        summary = pd.read_csv(summary_path, sep="\t")
        if "tf" in summary.columns:
            ranked_tfs = [tf for tf in summary["tf"].tolist() if tf in regulons]
    return {tf: regulons[tf] for tf in ranked_tfs[:top_n]}


def resolve_spatial_coordinates(adata: sc.AnnData) -> np.ndarray:
    if "spatial" in adata.obsm:
        coords = np.asarray(adata.obsm["spatial"], dtype=float)
    elif {"array_row", "array_col"}.issubset(adata.obs.columns):
        coords = adata.obs[["array_row", "array_col"]].to_numpy(dtype=float)
    elif "X_umap" in adata.obsm:
        coords = np.asarray(adata.obsm["X_umap"], dtype=float)
    else:
        coords = np.column_stack([np.arange(adata.n_obs, dtype=float), np.zeros(adata.n_obs, dtype=float)])
    if coords.ndim == 1:
        coords = coords[:, None]
    return coords[:, :2]


def build_knn_weights(coords: np.ndarray, neighbors: int) -> np.ndarray:
    if len(coords) <= 1:
        return np.zeros((len(coords), len(coords)), dtype=float)
    k = max(1, min(neighbors, len(coords) - 1))
    deltas = coords[:, None, :] - coords[None, :, :]
    distances = np.sqrt(np.square(deltas).sum(axis=2))
    np.fill_diagonal(distances, np.inf)
    nearest = np.argpartition(distances, kth=k - 1, axis=1)[:, :k]
    weights = np.zeros((len(coords), len(coords)), dtype=float)
    row_index = np.arange(len(coords))[:, None]
    # Symmetrize the k-nearest-neighbor graph so Moran's I uses reciprocal local neighborhoods.
    weights[row_index, nearest] = 1.0
    return np.maximum(weights, weights.T)


def morans_i(values: np.ndarray, weights: np.ndarray, permutations: int) -> tuple[float, float]:
    centered = np.asarray(values, dtype=float) - float(np.mean(values))
    denominator = float(np.square(centered).sum())
    total_weight = float(weights.sum())
    if denominator == 0 or total_weight == 0:
        return 0.0, 1.0
    observed = float((len(values) / total_weight) * ((centered @ weights @ centered) / denominator))
    if permutations <= 0:
        return observed, 1.0
    exceedances = 0
    for _ in range(permutations):
        permuted = np.random.permutation(centered)
        permuted_score = float((len(values) / total_weight) * ((permuted @ weights @ permuted) / denominator))
        if abs(permuted_score) >= abs(observed):
            exceedances += 1
    return observed, (exceedances + 1) / (permutations + 1)


def build_spatial_stats(scores: pd.DataFrame, coords: np.ndarray, neighbors: int, permutations: int) -> pd.DataFrame:
    if scores.empty:
        return pd.DataFrame(columns=["tf", "morans_i", "moran_p_value", "mean_activity", "variance"])
    weights = build_knn_weights(coords, neighbors)
    rows = []
    for tf in scores.columns:
        values = scores[tf].to_numpy(dtype=float)
        moran_score, p_value = morans_i(values, weights, permutations)
        rows.append(
            {
                "tf": tf,
                "morans_i": moran_score,
                "moran_p_value": p_value,
                "mean_activity": float(values.mean()),
                "variance": float(values.var()),
            }
        )
    return pd.DataFrame(rows).sort_values(["morans_i", "variance", "tf"], ascending=[False, False, True])


def assign_niches(
    adata: sc.AnnData,
    scores: pd.DataFrame,
    coords: np.ndarray,
    marker_map: dict[str, list[str]],
) -> pd.DataFrame:
    center = coords.mean(axis=0)
    radial_distance = np.sqrt(np.square(coords - center).sum(axis=1))
    radial_rank = pd.Series(radial_distance, index=adata.obs_names).rank(pct=True, method="average")
    malignant_score = mean_expression(adata, marker_map.get("malignant", DEFAULT_NICHE_MARKERS["malignant"]))
    immune_score = mean_expression(adata, marker_map.get("immune", DEFAULT_NICHE_MARKERS["immune"]))
    stromal_score = mean_expression(adata, marker_map.get("stromal", DEFAULT_NICHE_MARKERS["stromal"]))
    vascular_score = mean_expression(adata, marker_map.get("vascular", DEFAULT_NICHE_MARKERS["vascular"]))
    dominant_regulon = scores.idxmax(axis=1) if not scores.empty else pd.Series("NA", index=adata.obs_names)
    dominant_activity = scores.max(axis=1) if not scores.empty else pd.Series(0.0, index=adata.obs_names)

    niche_labels: list[str] = []
    for obs_name in adata.obs_names:
        idx = adata.obs_names.get_loc(obs_name)
        compartment_scores = {
            "vascular_niche": vascular_score[idx],
            "immune_stroma": immune_score[idx] + stromal_score[idx],
            "tumor_core": malignant_score[idx] - radial_rank.loc[obs_name],
            "invasive_margin": malignant_score[idx] + radial_rank.loc[obs_name],
        }
        niche_labels.append(max(compartment_scores, key=lambda label: compartment_scores[label]))

    domain = (
        adata.obs["domain"].astype(str)
        if "domain" in adata.obs
        else pd.Series(["domain_0"] * adata.n_obs, index=adata.obs_names)
    )
    return pd.DataFrame(
        {
            "spot_id": adata.obs_names,
            "regulon_domain": domain.to_numpy(),
            "niche_label": niche_labels,
            "dominant_regulon": dominant_regulon.to_numpy(),
            "dominant_regulon_activity": dominant_activity.to_numpy(),
            "radial_position": radial_rank.to_numpy(),
            "malignant_score": standardize(malignant_score),
            "immune_score": standardize(immune_score),
            "stromal_score": standardize(stromal_score),
            "vascular_score": standardize(vascular_score),
        }
    )


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--spatial", type=Path, required=True)
    parser.add_argument("--regulons", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--config", type=Path, default=PROJECT_ROOT / "config" / "analysis_config.yaml")
    parser.add_argument("--summary", type=Path)
    parser.add_argument("--top-n", type=int)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    logger = setup_logging()
    config = load_analysis_config(args.config)
    top_n = int(get_config_value(config, "analysis", "spatial", "top_regulons_for_mapping", default=50))
    neighbors = int(get_config_value(config, "analysis", "spatial", "moran_neighbors", default=6))
    permutations = int(get_config_value(config, "analysis", "spatial", "moran_permutations", default=99))
    marker_map = get_config_value(config, "analysis", "spatial", "niche_markers", default={})
    if args.top_n is not None:
        top_n = args.top_n
    adata = sc.read_h5ad(args.spatial)
    regulons = json.loads(args.regulons.read_text(encoding="utf-8"))
    selected_regulons = prioritize_regulons(regulons, args.summary, top_n)
    scores = score_regulons(adata, selected_regulons)
    coords = resolve_spatial_coordinates(adata)
    spatial_stats = build_spatial_stats(scores, coords, neighbors, permutations)
    niche_calls = assign_niches(adata, scores, coords, marker_map)
    args.output_dir.mkdir(parents=True, exist_ok=True)
    scores.to_csv(args.output_dir / "regulon_activity.tsv", sep="\t")
    spatial_stats.to_csv(args.output_dir / "regulon_spatial_stats.tsv", sep="\t", index=False)
    niche_calls.to_csv(args.output_dir / "niche_calls.tsv", sep="\t", index=False)
    logger.info("Saved scores for %s regulons across %s spots", scores.shape[1], scores.shape[0])


if __name__ == "__main__":
    main()
