#!/usr/bin/env python3
"""Generate a figure manifest for the planned manuscript figures."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[2]

FIGURES = [
    ("Figure1", "A", "Single-cell atlas and cohort integration", "UMAP and cohort composition", "data/processed/scrna_processed.h5ad"),
    ("Figure1", "B", "Single-cell atlas and cohort integration", "Compartment marker overview", "out/grn_networks/cell_type_grn_edges.tsv"),
    ("Figure2", "A", "Consensus GRN architecture and hub TFs", "Global hub TF ranking", "out/grn_networks/advanced_tf_support.tsv"),
    ("Figure2", "B", "Consensus GRN architecture and hub TFs", "Compartment-specific hub TFs", "out/grn_networks/cell_type_grn_edges.tsv"),
    ("Figure3", "A", "Spatial regulon maps and regulatory domains", "Top regulon spatial overlays", "out/spatial_mapping/regulon_activity.tsv"),
    ("Figure3", "B", "Spatial regulon maps and regulatory domains", "Moran's I and niche-restricted regulons", "out/spatial_mapping/regulon_spatial_stats.tsv"),
    ("Figure4", "A", "Niche-specific TF programs and immunotype contrasts", "Spot-level niche calls", "out/spatial_mapping/niche_calls.tsv"),
    ("Figure4", "B", "Niche-specific TF programs and immunotype contrasts", "Immune-stroma versus malignant contrasts", "out/validation/tf_compartment_specificity.tsv"),
    ("Figure5", "A", "Clinical and orthogonal public-data validation", "Cross-cohort and spatial support", "out/validation/advanced_validation_summary.tsv"),
    ("Figure5", "B", "Clinical and orthogonal public-data validation", "Spatially consistent regulons", "out/validation/tf_spatial_support.tsv"),
    ("Figure6", "A", "Therapeutic prioritization and working model", "Prioritized TF candidates", "out/validation/advanced_validation_summary.tsv"),
    ("Figure6", "B", "Therapeutic prioritization and working model", "Working model across HCC niches", "prompts/analysis_prompts.md"),
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["figure_id", "panel_id", "title", "analysis_focus", "primary_input"])
        writer.writerows(FIGURES)


if __name__ == "__main__":
    main()
