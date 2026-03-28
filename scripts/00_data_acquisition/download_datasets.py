#!/usr/bin/env python3
"""
Spatial GRN Atlas of HCC — Data Acquisition
Downloads scRNA-seq, scATAC-seq, and Visium datasets for GRN inference.
"""

import subprocess
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[2]
DATA_DIR = PROJECT_ROOT / "data"

DATASETS = {
    "scrna": [
        {"geo": "GSE149614", "study": "Sharma2020_HCC", "type": "scRNA-seq"},
        {"geo": "GSE151530", "study": "Sun2021_HCC", "type": "scRNA-seq"},
        {"geo": "GSE189903", "study": "Lu2022_multiregion", "type": "scRNA-seq"},
    ],
    "scatac": [
        # HCC multiome / scATAC datasets
        {"geo": "GSE243xxx", "study": "HCC_multiome", "type": "10x Multiome",
         "note": "Search GEO for recent HCC ATAC-seq datasets"},
    ],
    "visium": [
        # Use existing HCC_omics Visium data
        {"source": "local", "path": "~/projects/HCC_omics/data/spatial_transcriptome",
         "study": "HCC_omics_cohort", "type": "Visium"},
    ],
    "references": [
        {"name": "JASPAR2024", "url": "https://jaspar.elixir.no/download/data/2024/CORE/JASPAR2024_CORE_vertebrates_non-redundant_pfms_jaspar.txt"},
        {"name": "ENCODE_TF_ChIP", "url": "https://www.encodeproject.org/"},
    ]
}


def download_geo(geo_id: str, output_dir: Path):
    """Download from GEO."""
    output_dir.mkdir(parents=True, exist_ok=True)
    if any(output_dir.iterdir()):
        print(f"  [SKIP] {geo_id}")
        return

    ftp = f"https://ftp.ncbi.nlm.nih.gov/geo/series/{geo_id[:-3]}nnn/{geo_id}/suppl/"
    subprocess.run([
        "wget", "-r", "-np", "-nd", "-P", str(output_dir),
        "-A", "*.tar.gz,*.h5,*.h5ad,*.csv.gz,*.tsv.gz,*.mtx.gz",
        "--timeout=300", "--tries=3", ftp
    ], capture_output=True, timeout=600)


def main():
    print("=" * 70)
    print("Spatial GRN Atlas — Data Acquisition")
    print("=" * 70)

    for data_type, datasets in DATASETS.items():
        print(f"\n[{data_type.upper()}]")
        for ds in datasets:
            if ds.get("geo"):
                geo = ds["geo"]
                ds_dir = DATA_DIR / data_type / f"{geo}_{ds['study']}"
                if "xxx" not in geo:
                    download_geo(geo, ds_dir)
                else:
                    print(f"  [TODO] Search GEO for: {ds.get('note', ds['study'])}")
            elif ds.get("source") == "local":
                print(f"  [LOCAL] Will symlink from {ds['path']}")
            else:
                print(f"  [MANUAL] {ds['name']}: {ds.get('url', 'N/A')}")

    print("\n[DONE] Data acquisition complete.")


if __name__ == "__main__":
    main()
