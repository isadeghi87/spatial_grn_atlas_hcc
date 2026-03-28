#!/bin/bash
# Spatial GRN Atlas — Download HCC scRNA + spatial data
set -euo pipefail

DATA_DIR="$HOME/projects/spatial_grn_atlas_hcc/data"
LOG="$DATA_DIR/download_log.txt"
echo "$(date): Starting GRN Atlas downloads" > "$LOG"

geo_fetch() {
    local geo_id="$1" study="$2" subdir="$3"
    local out_dir="$DATA_DIR/$subdir/${geo_id}_${study}"
    mkdir -p "$out_dir"

    if [ "$(find "$out_dir" -type f \( -name '*.gz' -o -name '*.h5' -o -name '*.tar' \) 2>/dev/null | wc -l)" -gt 0 ]; then
        echo "  [SKIP] $geo_id ($study)" | tee -a "$LOG"; return
    fi

    echo "  [FETCH] $geo_id ($study)" | tee -a "$LOG"
    local prefix="${geo_id:0:$((${#geo_id}-3))}nnn"
    local base_url="https://ftp.ncbi.nlm.nih.gov/geo/series/${prefix}/${geo_id}/suppl"
    local tmp_idx="/tmp/geo_idx_${geo_id}.html"
    
    wget -q --no-check-certificate "$base_url/" -O "$tmp_idx" 2>/dev/null || { echo "    [FAIL] Cannot access $geo_id" | tee -a "$LOG"; return; }
    local files=$(grep -oP 'href="\K[^"]+' "$tmp_idx" | grep -v 'Parent\|geo/series\|hhs.gov' | grep -v '/$')
    
    for fname in $files; do
        [ -f "$out_dir/$fname" ] && continue
        echo "    [GET] $fname" | tee -a "$LOG"
        wget -q --no-check-certificate -P "$out_dir" "${base_url}/${fname}" 2>&1 || true
    done
    rm -f "$tmp_idx"
    echo "    [OK] $(find "$out_dir" -type f | wc -l) files" | tee -a "$LOG"
}

echo "================================================================="
echo "Spatial GRN Atlas HCC — Data Download"
echo "================================================================="

# --- HCC scRNA-seq ---
echo -e "\n=== HCC SCRNA ===" | tee -a "$LOG"
geo_fetch "GSE149614" "Ma2020_HCC" "scrna"
geo_fetch "GSE151530" "Sun2021_HCC" "scrna"
geo_fetch "GSE125449" "Ma2019_HCC" "scrna"

# --- HCC Spatial (Visium) ---
echo -e "\n=== HCC SPATIAL ===" | tee -a "$LOG"
geo_fetch "GSE203612" "Chen2022_HCC_Visium" "spatial"

# --- Chromatin accessibility for GRN ---
echo -e "\n=== ATAC-SEQ ===" | tee -a "$LOG"
geo_fetch "GSE182613" "Zhang2022_HCC_ATAC" "atac"

echo "$(date): GRN Atlas downloads complete" | tee -a "$LOG"
