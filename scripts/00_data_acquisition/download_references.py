#!/usr/bin/env python3
"""Download motif and reference resources used for GRN support."""

from __future__ import annotations

import argparse
import logging
import urllib.request
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[2]
REFERENCE_DIR = PROJECT_ROOT / "data" / "references"
LOG_DIR = PROJECT_ROOT / "logs"
RESOURCES = {
    "jaspar2024": "https://jaspar.elixir.no/download/data/2024/CORE/JASPAR2024_CORE_vertebrates_non-redundant_pfms_jaspar.txt",
    "hocomoco_v11": "https://hocomoco11.autosome.org/final_bundle/hocomoco11/full/HUMAN/mono/HOCOMOCOv11_full_HUMAN_mono_meme_format.meme",
}


def setup_logging() -> logging.Logger:
    LOG_DIR.mkdir(exist_ok=True)
    logger = logging.getLogger("download_references")
    logger.setLevel(logging.INFO)
    if logger.handlers:
        return logger
    formatter = logging.Formatter("%(asctime)s | %(levelname)s | %(message)s")
    fh = logging.FileHandler(LOG_DIR / "download_references.log")
    fh.setFormatter(formatter)
    sh = logging.StreamHandler()
    sh.setFormatter(formatter)
    logger.addHandler(fh)
    logger.addHandler(sh)
    return logger


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--dry-run", action="store_true")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    logger = setup_logging()
    REFERENCE_DIR.mkdir(parents=True, exist_ok=True)
    for name, url in RESOURCES.items():
        target = REFERENCE_DIR / Path(url).name
        logger.info("%s -> %s", name, target)
        if not args.dry_run:
            urllib.request.urlretrieve(url, target)
    manifest = REFERENCE_DIR / "reference_sources.tsv"
    manifest.write_text(
        "resource\turl\n" + "\n".join(f"{name}\t{url}" for name, url in RESOURCES.items()) + "\n",
        encoding="utf-8",
    )
    logger.info("Wrote %s", manifest)


if __name__ == "__main__":
    main()
