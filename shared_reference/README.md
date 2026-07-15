# Prognostic eRNA Reference Files

This directory contains BED files for the 10 prognostic eRNAs identified via LASSO-Cox regression.

## Files

| File | Description | Window | Genome | Status |
|------|-------------|--------|--------|--------|
| `prognostic_eRNA_10.bed` | Summit ±500 bp (1000 bp window) | ±500 bp | **hg38** ✅ | 主文件 |
| `prognostic_eRNA_10_2kb.bed` | Extended window for deepTools visualization | ±2000 bp | **hg38** ✅ | 主文件 |
| `prognostic_eRNA_10_summit.bed` | Raw summit single-base coordinates (motif analysis) | 1 bp | **hg38** ✅ | 主文件 |

## eRNA Summit Coordinates (hg38, liftOver-corrected)

| # | eRNA ID | hg38 Summit | hg19 Original | Name Tag |
|---|---------|-------------|---------------|----------|
| 1 | chr1:155158995 | chr1:155186519 | chr1:155158995 | eRNA_chr1_155158995 |
| 2 | chr3:11236700 | chr3:11195014 | chr3:11236700 | eRNA_chr3_11236700 |
| 3 | chr8:22624675 | chr8:22767162 | chr8:22624675 | eRNA_chr8_22624675 |
| 4 | chr3:138070534 | chr3:138351692 | chr3:138070534 | eRNA_chr3_138070534 |
| 5 | chr9:71398719 | chr9:68783803 | chr9:71398719 | eRNA_chr9_71398719 |
| 6 | chr9:114689796 | chr9:111927516 | chr9:114689796 | eRNA_chr9_114689796 |
| 7 | chr10:5531356 | chr10:5489393 | chr10:5531356 | eRNA_chr10_5531356 |
| 8 | chr9:71398939 | chr9:68784023 | chr9:71398939 | eRNA_chr9_71398939 |
| 9 | chr12:13371038 | chr12:13218104 | chr12:13371038 | eRNA_chr12_13371038 |
| 10 | chr10:5528926 | chr10:5486963 | chr10:5528926 | eRNA_chr10_5528926 |

## Table 2. LASSO Cox Regression Coefficients and Hazard Statistics

Source: LASSO-Cox regression analysis of prognostic eRNAs.

| eRNA | Coefficient | Hazard Ratio (HR) | P-value | Direction |
|------|-------------|--------------------|---------|-----------|
| chr1:155158995 | -0.123381 | 0.825276 | 0.0016 | Protective (HR < 1) |
| chr3:11236700 | 0.218885 | 1.480299 | 0.0019 | Risk (HR > 1) |
| chr8:22624675 | -0.241738 | 0.738873 | 0.0026 | Protective (HR < 1) |
| chr3:138070534 | -0.852960 | 0.549234 | 0.0040 | Protective (HR < 1) |
| chr9:71398719 | -0.205501 | 0.622034 | 0.0045 | Protective (HR < 1) |
| chr9:114689796 | -0.209878 | 0.701644 | 0.0046 | Protective (HR < 1) |
| chr10:5531356 | 0.216798 | 1.240605 | 0.0053 | Risk (HR > 1) |
| chr9:71398939 | -0.037728 | 0.658367 | 0.0084 | Protective (HR < 1) |
| chr12:13371038 | 0.359465 | 1.381464 | 0.0098 | Risk (HR > 1) |
| chr10:5528926 | 0.064139 | 1.264055 | 0.0100 | Risk (HR > 1) |

## Notes

- BED format (0-based, half-open): column 1 = chromosome, column 2 = start, column 3 = end, column 4 = name, column 5 = score (0), column 6 = strand (.).
- Coordinates are in **hg38** (corrected from original hg19 via liftOver on 2026-07-05).
- The two chr9 entries (71398719 and 71398939, hg19 notation) are distinct, independent eRNA loci.
