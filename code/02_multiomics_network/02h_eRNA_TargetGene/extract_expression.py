#!/usr/bin/env python3
"""
step2_extract_expression.py
Extract eRNA and mRNA expression values, align samples.
Output: eRNA_mRNA_expression.csv
"""
import pandas as pd
import numpy as np
import gzip
import os
import warnings
warnings.filterwarnings('ignore')

BASE_DIR = "../../../../Data_Source"
OUT_DIR = "analysis"
os.makedirs(OUT_DIR, exist_ok=True)

# === Load nearest gene mapping ===
mapping = pd.read_csv(os.path.join(OUT_DIR, "nearest_gene_mapping.csv"))
erna_ids = mapping['eRNA_name'].tolist()
gene_symbols = mapping['nearest_gene'].tolist()
gene_ens_ids = mapping['gene_ensembl_id'].tolist()

# eRNA row IDs in the expression matrix (hg19 1-based summit)
erna_row_ids = [
    'chr1:155158995', 'chr10:5528926', 'chr10:5531356',
    'chr12:13371038', 'chr3:11236700', 'chr3:138070534',
    'chr8:22624675', 'chr9:71398719', 'chr9:71398939', 'chr9:114689796'
]

unique_symbols = []
unique_ens_ids = []
for i, sid in enumerate(gene_ens_ids):
    if sid not in unique_ens_ids:
        unique_ens_ids.append(sid)
        unique_symbols.append(gene_symbols[i])

print(f"=== Step 2: Extract Expression ===")
print(f"eRNA row IDs: {erna_row_ids}")
print(f"Unique gene Ensembl IDs: {unique_ens_ids}")
print(f"Unique gene symbols: {unique_symbols}")

# === Load eRNA expression ===
print("\n>>> Loading eRNA expression matrix...")
erna_all = pd.read_csv(
    f"{BASE_DIR}/TCGA/TCGA_RPKM_eRNA_300k_peaks_in_Super_enhancer_BRCA.txt",
    sep='\t', index_col=0, nrows=None
)
# Extract only tumor samples
tumor_cols = [c for c in erna_all.columns if c.endswith('_tumor')]
print(f"   Total samples: {len(erna_all.columns)}, Tumor: {len(tumor_cols)}")

# Find rows matching our 10 eRNAs
erna_found = {}
for row_id in erna_row_ids:
    if row_id in erna_all.index:
        erna_found[row_id] = True
    else:
        erna_found[row_id] = False

print(f"   eRNA rows found: {sum(erna_found.values())}/{len(erna_found)}")
for rid, found in erna_found.items():
    print(f"     {rid}: {'✓' if found else '✗'}")

# Extract eRNA expression
erna_expr = erna_all.loc[[rid for rid in erna_row_ids if erna_found[rid]], tumor_cols]
erna_expr.index = [mapping[mapping['eRNA_name'] == f'eRNA_{rid.replace(":", "_")}']['eRNA_name'].values[0]
                    if f'eRNA_{rid.replace(":", "_")}' in mapping['eRNA_name'].values
                    else rid for rid in erna_expr.index]

# Rename columns to 12-char TCGA barcode
erna_expr.columns = [c.replace('_tumor', '')[:12] for c in erna_expr.columns]
print(f"   eRNA matrix shape: {erna_expr.shape}")

# === Load mRNA FPKM expression ===
print("\n>>> Loading mRNA FPKM matrix...")
all_rows = []
target_gene_ids = set(unique_ens_ids)
with gzip.open(f"{BASE_DIR}/TCGA/TCGA-BRCA.star_fpkm.tsv.gz", 'rt') as f:
    header = next(f).strip().split('\t')
    print(f"   FPKM total samples: {len(header) - 1}")
    for line in f:
        parts = line.strip().split('\t')
        gene_id = parts[0].split('.')[0]
        if gene_id in target_gene_ids:
            all_rows.append(parts)
            target_gene_ids.discard(gene_id)
            if not target_gene_ids:
                break

# Build mRNA expression DataFrame
mrna_cols = header[1:]
mrna_data = {}
for parts in all_rows:
    gene_id = parts[0].split('.')[0]
    symbol_idx = unique_ens_ids.index(gene_id)
    symbol = unique_symbols[symbol_idx]
    mrna_data[symbol] = np.array([float(x) if x else 0.0 for x in parts[1:]])

mrna_expr = pd.DataFrame(mrna_data, index=mrna_cols)
# Rename to 12-char TCGA barcode
mrna_expr.index = [x[:12] for x in mrna_expr.index]
# Deduplicate: average if multiple samples per patient
mrna_expr = mrna_expr.groupby(mrna_expr.index).mean()
print(f"   mRNA matrix shape: {mrna_expr.shape}")
print(f"   mRNA genes found: {mrna_expr.columns.tolist()}")

# === Align samples ===
erna_samples = set(erna_expr.columns)
mrna_samples = set(mrna_expr.index)
common_samples = sorted(erna_samples & mrna_samples)
print(f"\n>>> Sample alignment:")
print(f"   eRNA samples: {len(erna_samples)}")
print(f"   mRNA samples: {len(mrna_samples)}")
print(f"   Common: {len(common_samples)}")

erna_aligned = erna_expr[common_samples].T
mrna_aligned = mrna_expr.loc[common_samples]

# Combine
combined = pd.concat([erna_aligned, mrna_aligned], axis=1)
combined.index.name = 'Patient_ID'
combined.to_csv(os.path.join(OUT_DIR, "eRNA_mRNA_expression.csv"))
print(f"\n   Saved: eRNA_mRNA_expression.csv ({combined.shape})")
print(f"   Columns: {combined.columns.tolist()}")
print(f"   First patient: {combined.index[0]}")
print("✓ Step 2 complete")