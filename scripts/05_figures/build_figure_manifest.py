#!/usr/bin/env python3
"""Generate a figure manifest for the planned manuscript figures."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[2]

FIGURES = [
    ("Figure1", "Single-cell atlas and cohort integration"),
    ("Figure2", "Consensus GRN architecture and hub TFs"),
    ("Figure3", "Spatial regulon maps and regulatory domains"),
    ("Figure4", "Niche-specific TF programs and immunotype contrasts"),
    ("Figure5", "Clinical and orthogonal public-data validation"),
    ("Figure6", "Therapeutic prioritization and working model"),
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
        writer.writerow(["figure_id", "title"])
        writer.writerows(FIGURES)


if __name__ == "__main__":
    main()
