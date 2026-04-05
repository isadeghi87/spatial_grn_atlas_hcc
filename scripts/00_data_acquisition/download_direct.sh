#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"
DATA_DIR="${PROJECT_ROOT}/data"
LOG_DIR="${PROJECT_ROOT}/logs"
mkdir -p "${LOG_DIR}"
LOG_FILE="${LOG_DIR}/download_direct.log"

DRY_RUN=0
if [[ "${1:-}" == "--dry-run" ]]; then
  DRY_RUN=1
fi

run_cmd() {
  echo "$(date -Iseconds) | $*" | tee -a "${LOG_FILE}"
  if [[ "${DRY_RUN}" -eq 0 ]]; then
    "$@"
  fi
}

fetch_geo() {
  local accession="$1"
  local group="$2"
  local label="$3"
  local prefix="${accession:0:${#accession}-3}nnn"
  local url="https://ftp.ncbi.nlm.nih.gov/geo/series/${prefix}/${accession}/suppl/"
  local out_dir="${DATA_DIR}/${group}/${accession}_${label}"
  mkdir -p "${out_dir}"
  run_cmd wget -r -np -nd -N -P "${out_dir}" -A "*.gz,*.tar,*.tar.gz,*.h5,*.h5ad,*.mtx,*.tsv,*.csv" "${url}"
}

mkdir -p "${DATA_DIR}"/{scrna,spatial,atac,tcga,references}

fetch_geo GSE149614 scrna Ma2020_HCC_scRNA
fetch_geo GSE151530 scrna Sun2021_HCC_scRNA
fetch_geo GSE125449 scrna Ma2019_Early_HCC_scRNA
fetch_geo GSE156625 scrna Zhang2021_HCC_TME_scRNA
fetch_geo GSE203612 spatial Chen2022_HCC_Visium
fetch_geo GSE182613 atac Zhang2022_HCC_ATAC

run_cmd wget -N -P "${DATA_DIR}/references" "https://jaspar.elixir.no/download/data/2024/CORE/JASPAR2024_CORE_vertebrates_non-redundant_pfms_jaspar.txt"
run_cmd wget -N -P "${DATA_DIR}/references" "https://hocomoco11.autosome.org/final_bundle/hocomoco11/full/HUMAN/mono/HOCOMOCOv11_full_HUMAN_mono_meme_format.meme"

echo "$(date -Iseconds) | TCGA-LIHC is available via GDC portal/API: https://portal.gdc.cancer.gov/projects/TCGA-LIHC" | tee -a "${LOG_FILE}"
