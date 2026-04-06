# Spatial Gene Regulatory Network Atlas of HCC

![Python](https://img.shields.io/badge/Python-3.10%2B-blue) ![R](https://img.shields.io/badge/R-4.3%2B-276DC3) ![License: MIT](https://img.shields.io/badge/License-MIT-green) ![Status](https://img.shields.io/badge/Status-Active-brightgreen)

> Inferring spatially resolved gene regulatory networks in HCC to identify transcription factor programs restricted to tumor and microenvironmental niches, with public-data validation expected for a Nature Cancer-level study.

## Highlights

- Multi-method GRN inference (SCENIC+, CellOracle, Dictys)
- Spatial mapping of regulon activities onto Visium coordinates
- Niche-specific TF program discovery (tumor core, margin, stroma, perivascular)
- Immunotype-specific regulatory switch identification
- Therapeutic TF target prioritization with druggability scoring
- Config-driven compartment-aware GRN summaries, spatial autocorrelation metrics, and niche calls

**Target Journal:** Nature Cancer

## Repository Structure

```
spatial_grn_atlas_hcc/
├── config/                          # GRN inference parameters
├── data/
│   ├── references/                  # TF binding databases
│   ├── scatac/                      # scATAC-seq data
│   └── scrna/                       # scRNA-seq data
├── config/                          # Analysis parameters and public data registry
├── docs/                            # Documentation and public dataset manifest
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
│   ├── 01_preprocessing/
│   ├── 02_grn_inference/
│   ├── 03_spatial_mapping/
│   ├── 04_validation/
│   ├── 05_figures/
│   └── shared_analysis.py            # Shared config and scoring helpers
├── .gitignore
├── LICENSE
├── Makefile
└── README.md
```

## Datasets

### Discovery and validation cohorts
| Study | Accession | Data Type | Primary use |
|-------|-----------|-----------|-------------|
| Ma et al. | GSE149614 | scRNA-seq | discovery atlas |
| Sun et al. | GSE151530 | scRNA-seq | independent single-cell validation |
| Ma et al. early HCC | GSE125449 | scRNA-seq | early lesion validation |
| Zhang et al. TME | GSE156625 | scRNA-seq | immune-state validation |
| Chen et al. | GSE203612 | Visium spatial transcriptomics | spatial niche mapping |
| Zhang et al. | GSE182613 | ATAC-seq | motif and enhancer support |
| TCGA-LIHC | TCGA-LIHC | bulk RNA-seq + clinical | survival and clinical validation |
| Human Liver Cell Atlas | CZ CELLxGENE | scRNA-seq reference | annotation transfer |

### External priors and orthogonal support
- ENCODE TF ChIP-seq metadata
- CistromeDB human cancer TF binding metadata
- JASPAR 2024 CORE vertebrate motifs
- HOCOMOCO v11 human motifs

## Getting Started

### Prerequisites

- **Python 3.10+**: scanpy, squidpy, pyscenic, celloracle, scvi-tools, scrublet, lifelines
- **R 4.3+**: Seurat v5, SPATA2, Signac, AUCell, ComplexHeatmap

### Installation

```bash
# Clone the repository
git clone https://github.com/isadeghi87/spatial_grn_atlas_hcc.git
cd spatial_grn_atlas_hcc

conda env create -f envs/environment.yaml

# Activate
conda activate spatial_grn_atlas_hcc
```

### Quick Start

```bash
# 1. Download public data and reference resources
python3 scripts/00_data_acquisition/download_datasets.py --dry-run
python3 scripts/00_data_acquisition/download_references.py --dry-run

# 2. Review configurable analysis thresholds
cat config/analysis_config.yaml

# 3. Run the workflow
make all
```

## Advanced analysis outputs

The workflow now produces additional advanced-analysis tables beyond the baseline GRN edge list:

- `out/grn_networks/tf_priority_summary.tsv`: global TF ranking across the merged atlas
- `out/grn_networks/advanced_tf_support.tsv`: merged view of TF strength, compartment support, and cohort support
- `out/grn_networks/cell_type_grn_edges.tsv`: compartment-aware GRN edges using the `cell_type` labels
- `out/grn_networks/cohort_grn_edges.tsv`: cohort-aware GRN edges using the `cohort` labels
- `out/spatial_mapping/regulon_spatial_stats.tsv`: Moran's I and spot-level variance for prioritized regulons
- `out/spatial_mapping/niche_calls.tsv`: heuristic tumor-core, invasive-margin, immune-stroma, and vascular niche assignments
- `out/validation/advanced_validation_summary.tsv`: integrated validation table spanning GRN rank, compartment specificity, cohort support, and spatial restriction

These outputs are parameterized by `config/analysis_config.yaml`.

## Nature Cancer-level analysis scope

The repository now includes:

- an expanded analysis prompt in `prompts/analysis_prompts.md`
- a public dataset registry in `docs/public_data_registry.tsv`
- CLI entry points for data acquisition, preprocessing, GRN inference, spatial mapping, validation, and figure manifest generation
- config-aware advanced summaries for compartment-specific GRNs, cohort support, spatial structure, and niche calls

These pieces are intended to support:

1. multi-cohort single-cell integration
2. consensus GRN construction with compartment and cohort support
3. spatial regulon mapping with Moran's I ranking and heuristic niche annotation
4. orthogonal public-data validation summaries
5. publication figure planning with panel-level source tracking

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
