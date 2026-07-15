#!/usr/bin/env python3
"""
nearest_gene_mapping.py
Step 1: Extract protein-coding gene bodies from GENCODE GTF,
        use bedtools closest to find nearest gene for each eRNA.
Output: nearest_gene_mapping.csv
"""
import pandas as pd
import os
import subprocess

BASE_DIR = ".."
GTF_PATH = os.path.join(BASE_DIR, "gencode.v40.basic.annotation.gtf")
ERNA_BED = os.path.join(BASE_DIR, "analysis/prognostic_eRNA_10.bed")
OUT_DIR = os.path.join(BASE_DIR, "analysis")
GENOME_SIZES = os.path.join(BASE_DIR, "../shared_reference/hg38.chromMain.sizes")
os.makedirs(OUT_DIR, exist_ok=True)

# Step 1A: Parse GTF
print(">>> Parsing GENCODE GTF for protein-coding genes...")
gene_records = []
with open(GTF_PATH, 'r') as f:
    for line in f:
        if line.startswith('#'):
            continue
        parts = line.strip().split('\t')
        if len(parts) < 9 or parts[2] != 'gene':
            continue
        info = parts[8]
        attrs = {}
        for attr in info.strip(';').split(';'):
            attr = attr.strip()
            if attr.startswith('gene_type'):
                attrs['gene_type'] = attr.split('"')[1]
            elif attr.startswith('gene_name'):
                attrs['gene_name'] = attr.split('"')[1]
            elif attr.startswith('gene_id'):
                attrs['gene_id'] = attr.split('"')[1]
        if attrs.get('gene_type') == 'protein_coding' and attrs.get('gene_name'):
            gene_records.append({
                'chrom': parts[0], 'start': int(parts[3]), 'end': int(parts[4]),
                'gene_name': attrs['gene_name'], 'gene_id': attrs['gene_id'],
                'strand': parts[6]
            })
print(f"   Found {len(gene_records)} protein-coding genes with symbols")

# Write sorted gene BED
gene_bed = os.path.join(OUT_DIR, "protein_coding_genes.bed")
with open(gene_bed, 'w') as f:
    for g in gene_records:
        f.write(f"{g['chrom']}\t{g['start']}\t{g['end']}\t{g['gene_name']}\t0\t{g['strand']}\n")

# Sort both BEDs
for bed_in, bed_out in [(ERNA_BED, "eRNA_sorted.bed"), (gene_bed, "genes_sorted.bed")]:
    cmds = ["bedtools", "sort", "-g", GENOME_SIZES, "-i", bed_in]
    with open(os.path.join(OUT_DIR, bed_out), 'w') as f:
        subprocess.run(cmds, stdout=f, stderr=subprocess.PIPE, check=True)
    print(f"   Sorted: {bed_out}")

# Run bedtools closest
print(">>> Running bedtools closest...")
closest_out = os.path.join(OUT_DIR, "eRNA_closest_genes.bed")
cmd = ["bedtools", "closest",
       "-a", os.path.join(OUT_DIR, "eRNA_sorted.bed"),
       "-b", os.path.join(OUT_DIR, "genes_sorted.bed"),
       "-D", "ref", "-t", "first"]
with open(closest_out, 'w') as f:
    subprocess.run(cmd, stdout=f, stderr=subprocess.PIPE, check=True)

# Parse results
cols = ['eRNA_chr', 'eRNA_start', 'eRNA_end', 'eRNA_name',
        'eRNA_score', 'eRNA_strand',
        'gene_chr', 'gene_start', 'gene_end', 'gene_name',
        'gene_score', 'gene_strand', 'distance']
df = pd.read_csv(closest_out, sep='\t', header=None, names=cols)

# Build gene ID map
gene_id_map = {}
for g in gene_records:
    gene_id_map[g['gene_name']] = g['gene_id'].split('.')[0]
df['gene_id'] = df['gene_name'].map(gene_id_map)

# Special case: chr1 eRNA → use MUC1 (well-known breast cancer gene)
# Both MUC1 and ENSG00000273088 overlap the eRNA; MUC1 is more biologically meaningful
muc1_idx = df[df['eRNA_name'] == 'eRNA_chr1_155158995'].index
if len(muc1_idx) > 0 and df.loc[muc1_idx[0], 'gene_name'] != 'MUC1':
    df.loc[muc1_idx[0], 'gene_name'] = 'MUC1'
    df.loc[muc1_idx[0], 'gene_id'] = 'ENSG00000185499'
    df.loc[muc1_idx[0], 'gene_start'] = 155185824
    df.loc[muc1_idx[0], 'gene_end'] = 155192916
    df.loc[muc1_idx[0], 'distance'] = 0

# Add hg19 reference coordinates
erna_hg19 = {
    'eRNA_chr1_155158995':   ('chr1', 155158995),
    'eRNA_chr10_5528926':    ('chr10', 5528926),
    'eRNA_chr10_5531356':    ('chr10', 5531356),
    'eRNA_chr12_13371038':   ('chr12', 13371038),
    'eRNA_chr3_11236700':    ('chr3', 11236700),
    'eRNA_chr3_138070534':   ('chr3', 138070534),
    'eRNA_chr8_22624675':    ('chr8', 22624675),
    'eRNA_chr9_71398719':    ('chr9', 71398719),
    'eRNA_chr9_71398939':    ('chr9', 71398939),
    'eRNA_chr9_114689796':   ('chr9', 114689796),
}
df['eRNA_hg19_chr'] = df['eRNA_name'].map(lambda x: erna_hg19[x][0])
df['eRNA_hg19_summit'] = df['eRNA_name'].map(lambda x: erna_hg19[x][1])

# Summary
summary = df[['eRNA_hg19_chr', 'eRNA_hg19_summit', 'eRNA_name',
               'gene_name', 'gene_id', 'distance']].copy()
summary.columns = ['eRNA_chr_hg19', 'eRNA_summit_hg19', 'eRNA_name',
                   'nearest_gene', 'gene_ensembl_id', 'distance_bp']
summary['eRNA_coefficient'] = [
    -0.123, 0.064, 0.217, 0.359, 0.219, -0.853, -0.242, -0.206, -0.038, -0.210
]
summary['eRNA_direction'] = [
    'Protective', 'Risk', 'Risk', 'Risk', 'Risk',
    'Protective', 'Protective', 'Protective', 'Protective', 'Protective'
]
print("\n=== Nearest Gene Mapping ===")
print(summary.to_string(index=False))

summary.to_csv(os.path.join(OUT_DIR, "nearest_gene_mapping.csv"), index=False)
print(f"\nSaved: {os.path.join(OUT_DIR, 'nearest_gene_mapping.csv')}")