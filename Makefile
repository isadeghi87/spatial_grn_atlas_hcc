.PHONY: all download preprocess analyze figures clean

all: download preprocess analyze figures

download:
	bash scripts/00_data_acquisition/download_direct.sh

preprocess:
	python3 scripts/01_preprocessing/preprocess_grn.py

analyze:
	python3 scripts/02_analysis/infer_grn.py

figures:
	@echo "Figures generated during analysis steps"

clean:
	rm -rf data/processed/*.h5ad out/
