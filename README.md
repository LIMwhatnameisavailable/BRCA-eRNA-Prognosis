# TCGA-BRCA eRNA Prognostic Signature Project

## Project Overview
This project aims to identify and validate a prognostic enhancer RNA (eRNA) signature for Breast Cancer (BRCA) using multi-omics data from TCGA. The analysis pipeline includes differential expression analysis, LASSO-Cox modeling, survival analysis, immune/pathway enrichment (GSEA), and multi-omics regulatory mechanism exploration (CNV, Methylation, Mutation, Hi-C).

## Repository Structure
The project is organized into three main directories. Please maintain this structure to ensure code reproducibility.

```text
Prognosis_R/
├── Prognosis_R.Rproj         <-- Double-click this to open the project
├── README.md                 <-- Project documentation
├── Data_Source/              <-- Contains all raw input files
│   ├── TCGA_RPKM_eRNA...csv
│   ├── TCGA-BRCA.star_fpkm.tsv.gz
│   ├── clinical_info.tsv
│   └── ... (other omics data)
├── R_Scripts/                <-- Numbered analysis scripts
│   ├── 01_Differnetially_expressed_eRNA_Recognition.R
│   ├── 02_Risk_Score_Stratification&KM_Curve&ROC.R
│   ├── 03_Independent_Prognostic_Analysis_Cox.R
│   ├── 04_GSEA_Analysis.R
│   ├── 11_Comprehensive_Regulatory_Network_eRNA&TF&Gene.R
│   └── 20_PCA_mRNA_with_eRNA_Tumor&Normal.R
└── Results/                  <-- Automatically generated Figures (.svg/.tiff) and Tables

## Data Preparation (Crucial)

**Note:** Due to file size limits, raw data is NOT included in this repository. 
To reproduce the analysis, please create a folder named `Data_Source` in the root directory and place the following files inside.

** Filenames must match exactly:**

| File Category | Required Filename (Must be Exact) | Source / Description |
| :--- | :--- | :--- |
| **Expression** | `TCGA_RPKM_eRNA_300k_peaks_in_Super_enhancer_BRCA.csv` | eRNA Expression Matrix |
| **Expression** | `TCGA-BRCA.star_fpkm.tsv.gz` | mRNA Expression (Star Counts) |
| **Clinical** | `clinical_info.tsv` | Patient Clinical Data |
| **Mutation** | `PCAWG_WGS_mutations.tsv.gz` | WGS Mutation Data |
| **CNV** | `TCGA.BRCA.sampleMap_Gistic2_CopyNumber_Gistic2_all_thresholded.by_genes.gz` | Gistic2 Copy Number Variation |
| **Methylation** | `TCGA.BRCA.sampleMap_HumanMethylation450.gz` | 450k Methylation Array |
| **Annotation** | `probeMap_illuminaMethyl450_hg19_GPL16304_TCGAlegacy` | Methylation Probe Map |
| **Hi-C** | `loop_info.csv` | Chromatin Loop Data |
| **Reference** | `hg19ToHg38.over.chain.gz` | LiftOver Chain File |

How to Run 
Prerequisite: Please open the project by double-clicking Prognosis_R.Rproj. This ensures all file paths (e.g., Data_Source/...) work correctly relative to the root directory.

Step 1: Differential Expression & Identification
Script: R_Scripts/01_Differnetially_expressed_eRNA_Recognition.R
Input: eRNA expression matrix, Clinical data.
Output: Volcano Plot (Fig_Volcano.svg), Differentially Expressed eRNAs (DEEs).

Step 2: Model Construction & Validation (Core)
Script: R_Scripts/02_Risk_Score_Stratification&KM_Curve&ROC.R
Description: Performs LASSO regression to build the 10-eRNA signature. Validates using KM curves and Time-dependent ROC in Training, Testing, and Combined cohorts.
Output: KM Curves (Fig_KM_*.svg), ROC Curves (Fig_ROC_*.svg), Risk Score Distribution.

Step 3: Clinical Independence Analysis
Script: R_Scripts/03_Independent_Prognostic_Analysis_Cox.R
Description: Univariate and Multivariate Cox regression to confirm the signature is an independent prognostic factor.
Output: Cox Regression Tables (Table2_*.csv), Forest Plots.

Step 4: Functional Enrichment (GSEA)
Script: R_Scripts/04_GSEA_Analysis.R
Description: Performs GSEA to explore activated/suppressed pathways (KEGG) in High vs. Low risk groups.
Output: GSEA Multi-plots, Bubble plots, and Pathway Heatmaps.

Step 5: Multi-Omics Regulatory Mechanism
Script: R_Scripts/11_Comprehensive_Regulatory_Network_eRNA&TF&Gene.R
Description: Integrates WGS mutations, CNV (Gistic2), Methylation (450k), and Hi-C loops to construct the upstream regulatory network.
Output: Network node/edge files for Cytoscape (Network_*.csv), CNV Boxplots.

(Optional) Quality Control
Script: R_Scripts/20_PCA_mRNA_with_eRNA_Tumor&Normal.R
Description: Generates PCA plots to visualize sample clustering (Tumor vs. Normal).

Dependencies
The analysis was performed using R (v4.x). Key packages required:
Data Manipulation: data.table, dplyr, tibble
Survival Analysis: survival, survminer, timeROC, glmnet
Visualization: ggplot2, ComplexHeatmap, pheatmap, ggsci, svglite
Bioinformatics: limma, clusterProfiler, org.Hs.eg.db, GenomicRanges, rtracklayer

Notes
Data_Source: Ensure all large files (e.g., PCAWG_WGS_mutations.tsv.gz) are fully downloaded and placed in the Data_Source/ folder before running Step 5.

Results: All figures are saved in vector format (SVG/PDF) for publication quality.
