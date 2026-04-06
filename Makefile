.PHONY: all download references preprocess analyze spatial validate figures clean

CONFIG=config/analysis_config.yaml

all: download references preprocess analyze spatial validate figures

download:
	python3 scripts/00_data_acquisition/download_datasets.py

references:
	python3 scripts/00_data_acquisition/download_references.py

preprocess:
	python3 scripts/01_preprocessing/preprocess_scrna.py --input-dir data/scrna --output data/processed/scrna_processed.h5ad
	python3 scripts/01_preprocessing/preprocess_spatial.py --input-dir data/spatial --output data/processed/spatial_processed.h5ad

analyze:
	python3 scripts/02_grn_inference/infer_grn.py --config $(CONFIG) --rna data/processed/scrna_processed.h5ad --output-dir out/grn_networks

spatial:
	@if [ -f data/processed/spatial_processed.h5ad ] && [ -f out/grn_networks/regulons.json ]; then \
		python3 scripts/03_spatial_mapping/map_regulons.py --config $(CONFIG) --spatial data/processed/spatial_processed.h5ad --regulons out/grn_networks/regulons.json --summary out/grn_networks/advanced_tf_support.tsv --output-dir out/spatial_mapping; \
	else \
		echo "Skipping spatial mapping because required inputs are missing"; \
	fi

validate:
	python3 scripts/04_validation/validate_programs.py --config $(CONFIG) --grn out/grn_networks/grn_edges.tsv --compartment-grn out/grn_networks/cell_type_grn_edges.tsv --cohort-grn out/grn_networks/cohort_grn_edges.tsv --spatial-stats out/spatial_mapping/regulon_spatial_stats.tsv --niche-calls out/spatial_mapping/niche_calls.tsv --output-dir out/validation

figures:
	python3 scripts/05_figures/build_figure_manifest.py --output out/figures/figure_manifest.tsv

clean:
	rm -rf data/processed/*.h5ad out/
