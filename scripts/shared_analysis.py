#!/usr/bin/env python3
"""Shared helpers for config-aware spatial GRN analyses."""

from __future__ import annotations

from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import yaml

PROJECT_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_CONFIG_PATH = PROJECT_ROOT / "config" / "analysis_config.yaml"


def load_analysis_config(config_path: Path | None = None) -> dict[str, Any]:
    path = config_path or DEFAULT_CONFIG_PATH
    if not path.exists():
        return {}
    data = yaml.safe_load(path.read_text(encoding="utf-8"))
    return data if isinstance(data, dict) else {}


def get_config_value(config: dict[str, Any], *keys: str, default: Any = None) -> Any:
    current: Any = config
    for key in keys:
        if not isinstance(current, dict) or key not in current:
            return default
        current = current[key]
    return current


def expression_frame(adata: Any) -> pd.DataFrame:
    matrix = adata.X.toarray() if hasattr(adata.X, "toarray") else np.asarray(adata.X)
    return pd.DataFrame(matrix, index=adata.obs_names, columns=adata.var_names)


def mean_expression(adata: Any, genes: list[str]) -> np.ndarray:
    present = [gene for gene in genes if gene in adata.var_names]
    if not present:
        return np.zeros(adata.n_obs, dtype=float)
    matrix = adata[:, present].X
    mean_axis = matrix.mean(axis=1)
    return mean_axis.A1 if hasattr(mean_axis, "A1") else np.asarray(mean_axis).ravel()


def standardize(values: np.ndarray) -> np.ndarray:
    array = np.asarray(values, dtype=float)
    mean = float(array.mean())
    scale = float(array.std(ddof=0))
    if scale == 0:
        return np.zeros_like(array)
    return (array - mean) / scale
