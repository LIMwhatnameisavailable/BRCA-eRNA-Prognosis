# TCGA-BRCA eRNA Prognostic Signature Project

## 📌 Project Overview
This project identifies and validates a prognostic enhancer RNA (eRNA) signature for Breast Cancer (BRCA) using multi-omics data. The analysis pipeline integrates differential expression analysis, LASSO-Cox modeling, survival analysis, GSEA pathway enrichment, and multi-omics regulatory mechanism exploration (CNV, Methylation, Mutation, Hi-C).

## 📂 Repository Structure
The project is organized into three main directories. 
**Note:** `Data_Source` is excluded from the repo due to size limits.

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

## 💾 Data Preparation (Crucial)

**⚠️ Action Required:** Raw data is **NOT** included in this repository. To reproduce the analysis, please create a folder named `Data_Source` in the root directory and download the following files.

**Filenames must match EXACTLY as listed below:**

| File Category | Required Filename (Must be Exact) | Source / Database |
| :--- | :--- | :--- |
| **eRNA Expr** | `TCGA_RPKM_eRNA_300k_peaks_in_Super_enhancer_BRCA.csv` | **TCeA** (The Cancer eRNA Atlas) |
| **mRNA Expr** | `TCGA-BRCA.star_fpkm.tsv.gz` | **UCSC Xena** (TCGA-BRCA) |
| **Clinical** | `TCGA-BRCA.clinical.tsv` | **UCSC Xena** (Phenotype data) |
| **Clinical** | `clinical_info.tsv` | **UCSC Xena** (Survival data) |
| **Mutation** | `PCAWG_WGS_mutations.tsv.gz` | **UCSC Xena** |
| **CNV** | `TCGA.BRCA.sampleMap_Gistic2_CopyNumber_Gistic2_all_thresholded.by_genes.gz` | **UCSC Xena** (Gistic2) |
| **Methylation** | `TCGA.BRCA.sampleMap_HumanMethylation450.gz` | **UCSC Xena** (450k Array) |
| **Annotation** | `probeMap_illuminaMethyl450_hg19_GPL16304_TCGAlegacy` | **UCSC Xena** (Platform Map) |
| **Hi-C** | `loop_info.csv` | **ENCODE / 4D Nucleome** |
| **Reference** | `hg19ToHg38.over.chain.gz` | **UCSC Genome Browser** |
| **Reference** | `gencode.v36.annotation.gtf.gene.probemap` | **UCSC Xena** |

> **Note:** Scripts `01-03` primarily use `TCGA-BRCA.clinical.tsv`, while Script `11` uses `clinical_info.tsv`. Please ensure both are present.

## 🚀 How to Run

**Prerequisite:** Open the project by double-clicking **`Prognosis_R.Rproj`**. This sets the working directory correctly.

### Step 1: Differential Expression & Identification
* **Script:** `R_Scripts/01_Differnetially_expressed_eRNA_Recognition.R`
* **Output:** Volcano Plot (`Fig_Volcano.svg`), Differentially Expressed eRNAs.

### Step 2: Model Construction & Validation (Core)
* **Script:** `R_Scripts/02_Risk_Score_Stratification&KM_Curve&ROC.R`
* **Function:** LASSO regression (10-eRNA signature), KM Survival Curves, Time-dependent ROC.
* **Output:** `Fig_KM_*.svg`, `Fig_ROC_*.svg`, `Fig_stratification_*.svg`.

### Step 3: Clinical Independence Analysis
* **Script:** `R_Scripts/03_Independent_Prognostic_Analysis_Cox.R`
* **Function:** Univariate/Multivariate Cox regression, Forest plots.
* **Output:** `Fig_Forest_Plot_Cox.svg`.

### Step 4: Functional Enrichment (GSEA)
* **Script:** `R_Scripts/04_GSEA_Analysis.R`
* **Function:** KEGG pathway enrichment (Activated vs Suppressed).
* **Output:** `Fig_GSEA_Dotplot.svg`, `Fig_GSEA_Multiplot_*.svg`.

### Step 5: Multi-Omics Regulatory Mechanism
* **Script:** `R_Scripts/11_Comprehensive_Regulatory_Network_eRNA&TF&Gene.R`
* **Function:** Integrates WGS mutations, CNV, Methylation, and Hi-C loops to build a regulatory network.
* **Output:** `Network_Edges_for_Cytoscape.csv`, `Fig_CNV_Boxplots.svg`.

### (Optional) Quality Control
* **Script:** `R_Scripts/20_PCA_mRNA_with_eRNA_Tumor&Normal.R`
* **Output:** `Fig_PCA_mRNA&eRNA.svg`.

## 🛠 Dependencies
* **R Version:** 4.x
* **Key Packages:** `data.table`, `dplyr`, `survival`, `survminer`, `glmnet`, `timeROC`, `ggplot2`, `ComplexHeatmap`, `limma`, `clusterProfiler`, `rtracklayer`.

