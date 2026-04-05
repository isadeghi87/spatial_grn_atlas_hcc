# Spatial GRN Atlas — Nature Cancer-Level Analysis Prompt

## Core objective
Build a hepatocellular carcinoma atlas that resolves transcription factor programs across malignant, immune, stromal, and vascular compartments; maps these programs into spatial niches; and validates them with chromatin, public cohorts, and clinical outcome data.

## Public data to use
- Discovery scRNA-seq: GSE149614, GSE151530, GSE125449, GSE156625
- Spatial transcriptomics: GSE203612
- Chromatin support: GSE182613
- Clinical validation: TCGA-LIHC
- Reference priors: JASPAR2024, HOCOMOCO v11, ENCODE, CistromeDB, Human Liver Cell Atlas

## Prompt 1: Multi-cohort single-cell discovery
1. Harmonize all public scRNA-seq cohorts and remove low-quality cells and doublets.
2. Annotate malignant, T, B, myeloid, endothelial, fibroblast, cholangiocyte, and hepatocyte-like states.
3. Derive malignant meta-programs, immune activation/exhaustion programs, and stromal activation states.
4. Quantify cohort consistency of cell states and marker genes.
5. Build lineage and transition analyses for malignant and immune compartments.

## Prompt 2: GRN inference and prioritization
1. Infer cell-type-specific TF-target relationships from expression correlations with optional motif support.
2. Incorporate chromatin evidence from public ATAC-seq for enhancer and motif support.
3. Build regulons for each TF and rank them by effect size, specificity, and support.
4. Derive consensus hub TFs per compartment and across cohorts.
5. Flag candidate oncogenic, immune-evasive, and stromal-remodeling regulators.

## Prompt 3: Spatial regulon mapping
1. Score regulon activity per Visium spot using AUCell when available, otherwise robust mean-expression scoring.
2. Compute Moran's I and neighborhood enrichment to identify spatially structured regulons.
3. Define regulatory domains and align them with histology-informed tumor core, invasive margin, immune stroma, and vascular niches.
4. Link niche-level regulon activity to inferred cell composition and ligand-receptor context.
5. Generate publication-grade spatial overlays for prioritized TF programs.

## Prompt 4: Advanced validations expected at journal level
1. Cross-cohort validation: confirm hub TFs and regulon target overlap in at least one independent single-cell cohort.
2. Chromatin validation: require motif or peak support for prioritized TF-target edges when ATAC evidence exists.
3. Clinical validation: test TCGA-LIHC associations with overall survival, stage, vascular invasion, and mutation context.
4. Specificity validation: compare malignant-specific regulons against immune and stromal compartments.
5. Spatial validation: confirm niche restriction of the same regulon across multiple sections or samples when available.
6. Literature validation: annotate whether top TFs are known in HCC biology or represent novel candidates.

## Prompt 5: Nature Cancer-style deliverables
1. Figure 1: cohort integration, cell atlas, and lineage overview.
2. Figure 2: consensus GRN architecture and compartment-specific hub TFs.
3. Figure 3: spatial distribution of regulon activity and regulatory domains.
4. Figure 4: niche-specific TF programs with immune-excluded versus inflamed contrasts.
5. Figure 5: TCGA and orthogonal validation of prioritized TF programs.
6. Figure 6: therapeutic prioritization and a model linking TF programs to HCC ecosystem states.
7. Extended data: QC, sensitivity analyses, cohort replication, and alternate thresholds.
