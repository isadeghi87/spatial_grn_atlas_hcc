#!/usr/bin/env python3
"""
Spatial GRN Atlas — GRN Inference using SCENIC+ and CellOracle
"""

import scanpy as sc
import numpy as np
import pandas as pd
from pathlib import Path
import warnings
warnings.filterwarnings("ignore")

PROJECT_ROOT = Path(__file__).resolve().parents[2]
DATA_DIR = PROJECT_ROOT / "data"
OUT_DIR = PROJECT_ROOT / "out" / "grn_networks"
OUT_DIR.mkdir(parents=True, exist_ok=True)


def run_scenic_plus(adata_rna, adata_atac=None):
    """Run SCENIC+ for GRN inference."""
    print("[SCENIC+] Starting GRN inference...")

    try:
        import scenicplus
        from scenicplus.scenicplus_class import SCENICPLUS
        from scenicplus.grn_builder.gsea import TF_link_genes_to_targets

        # SCENIC+ requires joint scRNA + scATAC
        if adata_atac is None:
            print("  [WARN] No scATAC data provided. Using motif-only mode.")

        # Step 1: cisTopic for topic modeling on ATAC peaks
        # Step 2: pycisTarget for motif enrichment
        # Step 3: Link TFs to target genes via correlation + motif evidence

        # This is a stub — full implementation requires ~200 lines of SCENIC+ API calls
        print("  [TODO] Implement full SCENIC+ pipeline")
        print("  See: https://scenicplus.readthedocs.io/")

        return None

    except ImportError:
        print("  [WARN] scenicplus not installed. Install via:")
        print("    pip install scenicplus")
        return None


def run_celloracle(adata_rna, adata_atac=None):
    """Run CellOracle for causal GRN inference."""
    print("[CellOracle] Starting GRN inference...")

    try:
        import celloracle as co

        # Step 1: Base GRN from ATAC peaks + motif scanning
        if adata_atac is not None:
            print("  Using ATAC peaks for base GRN construction")
            # Scan for TF binding motifs in ATAC peaks
            # base_GRN = co.motif_analysis.scan_motifs(...)
        else:
            print("  Using default promoter-based base GRN")

        # Step 2: Contextualize with expression
        oracle = co.Oracle()

        # Preprocess
        sc.pp.normalize_total(adata_rna, target_sum=1e4)
        sc.pp.log1p(adata_rna)
        sc.pp.highly_variable_genes(adata_rna, n_top_genes=3000)

        oracle.import_anndata_as_raw_count(
            adata=adata_rna,
            cluster_column_name="cell_type",
            embedding_name="X_umap"
        )

        # Step 3: Fit GRN per cell type
        # oracle.fit_GRN_for_simulation(alpha=10, use_cluster_specific_TFinfo=True)

        print("  [TODO] Complete CellOracle GRN fitting")
        print("  See: https://morris-lab.github.io/CellOracle.documentation/")

        return oracle

    except ImportError:
        print("  [WARN] celloracle not installed. Install via:")
        print("    pip install celloracle")
        return None


def score_regulons_spatial(adata_visium, regulons: dict, out_dir: Path):
    """Score regulon activities per Visium spot using AUCell."""
    print("[AUCell] Scoring regulon activities on spatial data...")

    try:
        from ctxcore.rnkdb import FeatureRankingDatabase
        from pyscenic.aucell import aucell

        # Create gene sets from regulons
        # AUCell scoring per spot
        # This gives regulon × spot activity matrix

        print("  [TODO] Implement AUCell scoring on Visium data")

    except ImportError:
        # Fallback: simple mean expression scoring
        print("  Using mean expression scoring (AUCell not available)")

        if not regulons:
            print("  [SKIP] No regulons provided")
            return

        for tf_name, target_genes in regulons.items():
            available = [g for g in target_genes if g in adata_visium.var_names]
            if available:
                adata_visium.obs[f"regulon_{tf_name}"] = np.mean(
                    adata_visium[:, available].X.toarray() if hasattr(adata_visium.X, 'toarray')
                    else adata_visium[:, available].X,
                    axis=1
                )

        print(f"  Scored {len(regulons)} regulons on {adata_visium.n_obs} spots")


def main():
    print("=" * 70)
    print("Spatial GRN Atlas — GRN Inference Pipeline")
    print("=" * 70)

    # Load data
    processed_dir = DATA_DIR / "processed"
    rna_files = list((DATA_DIR / "scrna").rglob("*.h5ad"))
    atac_files = list((DATA_DIR / "scatac").rglob("*.h5ad"))
    visium_files = list((DATA_DIR / "visium").rglob("*.h5ad"))

    print(f"Found: {len(rna_files)} scRNA, {len(atac_files)} scATAC, {len(visium_files)} Visium files")

    if not rna_files:
        print("[WARN] No scRNA-seq data found. Run data acquisition first.")
        return

    # Load first RNA dataset
    adata_rna = sc.read_h5ad(rna_files[0])
    adata_atac = sc.read_h5ad(atac_files[0]) if atac_files else None

    # Run GRN inference
    run_scenic_plus(adata_rna, adata_atac)
    oracle = run_celloracle(adata_rna, adata_atac)

    # Score regulons on spatial data
    if visium_files:
        adata_vis = sc.read_h5ad(visium_files[0])
        # Example regulons (would come from SCENIC+ output)
        example_regulons = {
            "MYC": ["CCND1", "CDK4", "E2F1", "MCM2", "PCNA"],
            "STAT1": ["IRF1", "GBP1", "CXCL10", "IDO1", "TAP1"],
            "HIF1A": ["VEGFA", "LDHA", "PGK1", "ENO1", "SLC2A1"],
        }
        score_regulons_spatial(adata_vis, example_regulons, OUT_DIR)

    print("\n[DONE] GRN inference pipeline complete.")


if __name__ == "__main__":
    main()
