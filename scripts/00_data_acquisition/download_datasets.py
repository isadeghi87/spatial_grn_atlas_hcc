#!/usr/bin/env python3
"""Download public datasets required for the spatial GRN atlas workflow."""

from __future__ import annotations

import argparse
import json
import logging
import subprocess
from pathlib import Path
from typing import Iterable

PROJECT_ROOT = Path(__file__).resolve().parents[2]
DATA_DIR = PROJECT_ROOT / "data"
LOG_DIR = PROJECT_ROOT / "logs"

DATASETS = [
    {
        "name": "Ma2020_HCC_scRNA",
        "group": "scrna",
        "accession": "GSE149614",
        "source": "geo",
        "url": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE149nnn/GSE149614/suppl/",
    },
    {
        "name": "Sun2021_HCC_scRNA",
        "group": "scrna",
        "accession": "GSE151530",
        "source": "geo",
        "url": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE151nnn/GSE151530/suppl/",
    },
    {
        "name": "Ma2019_Early_HCC_scRNA",
        "group": "scrna",
        "accession": "GSE125449",
        "source": "geo",
        "url": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE125nnn/GSE125449/suppl/",
    },
    {
        "name": "Zhang2021_HCC_TME_scRNA",
        "group": "scrna",
        "accession": "GSE156625",
        "source": "geo",
        "url": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE156nnn/GSE156625/suppl/",
    },
    {
        "name": "Chen2022_HCC_Visium",
        "group": "spatial",
        "accession": "GSE203612",
        "source": "geo",
        "url": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE203nnn/GSE203612/suppl/",
    },
    {
        "name": "Zhang2022_HCC_ATAC",
        "group": "atac",
        "accession": "GSE182613",
        "source": "geo",
        "url": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE182nnn/GSE182613/suppl/",
    },
    {
        "name": "TCGA_LIHC",
        "group": "tcga",
        "accession": "TCGA-LIHC",
        "source": "portal",
        "url": "https://portal.gdc.cancer.gov/projects/TCGA-LIHC",
    },
]


def setup_logging() -> logging.Logger:
    LOG_DIR.mkdir(exist_ok=True)
    logger = logging.getLogger("download_datasets")
    logger.setLevel(logging.INFO)
    if logger.handlers:
        return logger
    formatter = logging.Formatter("%(asctime)s | %(levelname)s | %(message)s")
    file_handler = logging.FileHandler(LOG_DIR / "download_datasets.log")
    file_handler.setFormatter(formatter)
    stream_handler = logging.StreamHandler()
    stream_handler.setFormatter(formatter)
    logger.addHandler(file_handler)
    logger.addHandler(stream_handler)
    return logger


def run_command(command: list[str], dry_run: bool, logger: logging.Logger) -> None:
    logger.info("CMD %s", " ".join(command))
    if dry_run:
        return
    subprocess.run(command, check=True)


def iter_datasets(groups: set[str]) -> Iterable[dict]:
    for dataset in DATASETS:
        if not groups or dataset["group"] in groups:
            yield dataset


def download_geo_dataset(dataset: dict, dry_run: bool, overwrite: bool, logger: logging.Logger) -> None:
    output_dir = DATA_DIR / dataset["group"] / f"{dataset['accession']}_{dataset['name']}"
    output_dir.mkdir(parents=True, exist_ok=True)
    if not overwrite and any(output_dir.iterdir()):
        logger.info("Skipping populated dataset directory: %s", output_dir)
        return
    command = [
        "wget",
        "-r",
        "-np",
        "-nd",
        "-N",
        "-P",
        str(output_dir),
        "-A",
        "*.gz,*.tar,*.tar.gz,*.h5,*.h5ad,*.tsv,*.csv,*.mtx",
        dataset["url"],
    ]
    run_command(command, dry_run, logger)


def write_manifest(selected: list[dict], manifest_path: Path, logger: logging.Logger) -> None:
    manifest_path.parent.mkdir(parents=True, exist_ok=True)
    manifest_path.write_text(json.dumps(selected, indent=2), encoding="utf-8")
    logger.info("Wrote manifest to %s", manifest_path)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--groups",
        nargs="*",
        default=[],
        choices=sorted({dataset["group"] for dataset in DATASETS}),
        help="Subset of dataset groups to download.",
    )
    parser.add_argument("--dry-run", action="store_true", help="Print actions only.")
    parser.add_argument("--overwrite", action="store_true", help="Overwrite populated directories.")
    parser.add_argument(
        "--manifest-out",
        type=Path,
        default=PROJECT_ROOT / "docs" / "download_manifest.json",
        help="Where to write the selected dataset manifest.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    logger = setup_logging()
    selected = list(iter_datasets(set(args.groups)))
    for dataset in selected:
        if dataset["source"] == "geo":
            download_geo_dataset(dataset, args.dry_run, args.overwrite, logger)
        else:
            logger.info("Manual/API dataset: %s -> %s", dataset["accession"], dataset["url"])
    write_manifest(selected, args.manifest_out, logger)


if __name__ == "__main__":
    main()
