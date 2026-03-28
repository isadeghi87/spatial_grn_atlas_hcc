# Spatial GRN Atlas — Analysis Prompts

## Prompt 1: GRN Inference
1. Download HCC Multiome (scRNA + scATAC) datasets
2. Run SCENIC+ for joint GRN inference (TF → enhancer → target gene)
3. Run CellOracle as alternative (uses ATAC peaks for causal edges)
4. Run Dictys for dynamic GRN inference per cell type
5. Intersect networks: keep edges with support from ≥2 methods
6. Identify cell-type-specific regulons (top 500 TF-target gene modules)

## Prompt 2: Spatial Regulon Mapping
1. Score regulon activities per Visium spot using AUCell
2. Create spatial activity maps for top 50 TFs
3. Compute spatial autocorrelation of regulon activities (Moran's I)
4. Cluster spots by regulon profiles → "regulatory domains"
5. Overlay regulatory domains on histological annotation
6. Compare regulatory domains with immunotype zones

## Prompt 3: Niche-Specific TF Programs
1. For each spatial niche (tumor core, border, immune stroma, perivascular), identify enriched regulons
2. Compute niche specificity score (entropy-based)
3. Key comparisons: immune-excluded vs infiltrated niche regulons
4. In silico TF perturbation: CellOracle knockout simulation
5. Predict effect of TF activation/repression on immunotype conversion
6. Cross-validate with TCGA-LIHC: do niche-specific TF programs predict survival?
