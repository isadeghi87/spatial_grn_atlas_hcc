# Spatial Gene Regulatory Network Atlas of HCC

![Python](https://img.shields.io/badge/Python-3.10%2B-blue) ![R](https://img.shields.io/badge/R-4.3%2B-276DC3) ![License: MIT](https://img.shields.io/badge/License-MIT-green) ![Status](https://img.shields.io/badge/Status-Active-brightgreen)

> Inferring spatially-resolved gene regulatory networks in HCC to identify transcription factor programs restricted to specific TME niches.

## Highlights

- Multi-method GRN inference (SCENIC+, CellOracle, Dictys)
- Spatial mapping of regulon activities onto Visium coordinates
- Niche-specific TF program discovery (tumor core, margin, stroma, perivascular)
- Immunotype-specific regulatory switch identification
- Therapeutic TF target prioritization with druggability scoring

**Target Journal:** Nature Genetics / Molecular Cell / Cancer Cell

## Repository Structure

```
spatial_grn_atlas_hcc/
├── config/                          # GRN inference parameters
├── data/
│   ├── references/                  # TF binding databases
│   ├── scatac/                      # scATAC-seq data
│   └── scrna/                       # scRNA-seq data
├── docs/                            # Documentation
├── envs/
│   └── environment.yaml             # Conda specs
├── notebooks/                       # Analysis notebooks
├── out/
│   ├── figures/                     # Publication figures
│   ├── grn_networks/                # Inferred GRN outputs
│   └── niche_programs/              # Niche-specific TF programs
├── prompts/                         # Analysis prompt templates
├── scripts/
│   ├── 00_data_acquisition/
│   │   ├── download_datasets.py
│   │   └── download_direct.sh
│   └── 02_grn_inference/
│       └── infer_grn.py
├── .gitignore
├── LICENSE
├── Makefile
└── README.md
```

## Datasets

### scRNA-seq + scATAC-seq
| Study | GEO | Data Type | Cells |
|-------|-----|-----------|-------|
| HCC Multiome | GSE243xxx | 10x Multiome | ~30k |
| TCGA-ATAC | TCGA | ATAC-seq bulk | 371 |
| Liver Cell Atlas | CZ CELLxGENE | scRNA-seq | ~2M |
| Sharma et al. | GSE149614 | scRNA-seq | ~60k |
| Sun et al. | GSE151530 | scRNA-seq | ~56k |

### TF Binding Databases
- ENCODE (ChIP-seq peaks for 500+ TFs)
- CistromeDB (cancer-specific TF binding)
- JASPAR & HOCOMOCO (TF motifs)

## Getting Started

### Prerequisites

- **Python 3.10+**: scenic+, celloracle, dictys, pyscenic, scanpy, squidpy
- **R 4.3+**: Seurat v5, SPATA2, Signac, AUCell, ComplexHeatmap

### Installation

```bash
# Clone the repository
git clone https://github.com/isadeghi87/spatial_grn_atlas_hcc.git
cd spatial_grn_atlas_hcc

# Create conda environment (if environment.yaml exists)
conda env create -f environment.yaml
# OR use envs/ specs
conda env create -f envs/environment.yaml

# Activate
conda activate spatial_grn_atlas_hcc
```

### Quick Start

```bash
# 1. Download data

# 2. Run analysis (if Makefile present)
make all
```

## Contributing

Contributions are welcome! Please:

1. Fork the repository
2. Create a feature branch (`git checkout -b feature/your-feature`)
3. Commit changes (`git commit -am 'Add feature'`)
4. Push to the branch (`git push origin feature/your-feature`)
5. Open a Pull Request

## Citation

If you use this work, please cite:

```bibtex
@article{sadeghi2026spatial_grn_atlas_hcc,
  title={Spatial Gene Regulatory Network Atlas of HCC},
  author={Sadeghi, Iman},
  year={2026},
  note={Manuscript in preparation}
}
```

## License

This project is licensed under the MIT License — see [LICENSE](LICENSE) for details.

## Contact

- **Iman Sadeghi** — [GitHub](https://github.com/isadeghi87)
