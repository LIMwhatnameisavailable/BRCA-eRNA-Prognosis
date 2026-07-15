#!/usr/bin/env python3
"""
step3_spearman_correlation.py
Spearman correlation between eRNA and nearest mRNA expression.
Output: spearman_correlation_results.csv, Figure_correlation.{pdf,png}
"""
import pandas as pd
import numpy as np
from scipy.stats import spearmanr
from statsmodels.stats.multitest import multipletests
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
import os, warnings
warnings.filterwarnings('ignore')

OUT_DIR = "analysis"

# Color palette (manuscript style)
RISK_COLOR = '#E87A5D'   # warm orange for Risk eRNAs
PROT_COLOR = '#5B7B9A'   # blue-grey for Protective eRNAs

# === Load data ===
df = pd.read_csv(os.path.join(OUT_DIR, "eRNA_mRNA_expression.csv"), index_col=0)
print(f"=== Step 3: Spearman Correlation ===")
print(f"Loaded {df.shape[0]} samples, {df.shape[1]} features")

# Define pairs: eRNA column name → mRNA column name
erna_gene_map = {
    'eRNA_chr1_155158995': 'MUC1',
    'eRNA_chr10_5528926': 'CALML5',
    'eRNA_chr10_5531356': 'CALML5',
    'eRNA_chr12_13371038': 'EMP1',
    'eRNA_chr3_11236700': 'HRH1',
    'eRNA_chr3_138070534': 'MRAS',
    'eRNA_chr8_22624675': 'PEBP4',
    'eRNA_chr9_71398719': 'PIP5K1B',
    'eRNA_chr9_71398939': 'PIP5K1B',
    'eRNA_chr9_114689796': 'UGCG',
}

# eRNA direction and coefficient
erna_info = {
    'eRNA_chr1_155158995': ('chr1:155158995', -0.123, 'Protective'),
    'eRNA_chr10_5528926': ('chr10:5528926', 0.064, 'Risk'),
    'eRNA_chr10_5531356': ('chr10:5531356', 0.217, 'Risk'),
    'eRNA_chr12_13371038': ('chr12:13371038', 0.359, 'Risk'),
    'eRNA_chr3_11236700': ('chr3:11236700', 0.219, 'Risk'),
    'eRNA_chr3_138070534': ('chr3:138070534', -0.853, 'Protective'),
    'eRNA_chr8_22624675': ('chr8:22624675', -0.242, 'Protective'),
    'eRNA_chr9_71398719': ('chr9:71398719', -0.206, 'Protective'),
    'eRNA_chr9_71398939': ('chr9:71398939', -0.038, 'Protective'),
    'eRNA_chr9_114689796': ('chr9:114689796', -0.210, 'Protective'),
}

# === Compute Spearman correlations ===
results = []
for erna_col, mrna_gene in erna_gene_map.items():
    x = df[erna_col].values
    y = df[mrna_gene].values
    rho, pval = spearmanr(x, y)
    coord, coef, direction = erna_info[erna_col]
    results.append({
        'eRNA_ID': coord,
        'eRNA_column': erna_col,
        'nearest_gene': mrna_gene,
        'spearman_rho': rho,
        'p_value': pval,
        'LASSO_coefficient': coef,
        'eRNA_direction': direction,
        'n_samples': len(x)
    })

res_df = pd.DataFrame(results)

# FDR correction
reject, p_corrected, _, _ = multipletests(res_df['p_value'], method='fdr_bh')
res_df['p_adjusted_BH'] = p_corrected
res_df['significant'] = reject

res_df.to_csv(os.path.join(OUT_DIR, "spearman_correlation_results.csv"), index=False)
print("\n=== Spearman Correlation Results ===")
print(res_df[['eRNA_ID', 'nearest_gene', 'spearman_rho', 'p_value', 'p_adjusted_BH', 'significant']].to_string(index=False))

# === Figure 1: Correlation scatter plot (2×5 grid) ===
fig, axes = plt.subplots(2, 5, figsize=(20, 8))
fig.suptitle('eRNA vs Nearest Gene mRNA Expression (Spearman Correlation)', fontsize=16, y=1.02)

for idx, (erna_col, mrna_gene) in enumerate(erna_gene_map.items()):
    row = idx // 5
    col = idx % 5
    ax = axes[row, col]

    x = df[erna_col].values
    y = df[mrna_gene].values
    rho, pval = spearmanr(x, y)
    coord, coef, direction = erna_info[erna_col]

    # Color by direction
    color = RISK_COLOR if direction == 'Risk' else PROT_COLOR
    ax.scatter(x, y, alpha=0.3, s=5, color=color, edgecolors='none')

    # Regression line
    mask = ~(np.isnan(x) | np.isnan(y))
    if mask.sum() > 2:
        m, b = np.polyfit(x[mask], y[mask], 1)
        x_line = np.linspace(x[mask].min(), x[mask].max(), 100)
        ax.plot(x_line, m * x_line + b, color='#333333', linewidth=1.5, linestyle='-')

    # Format P-value
    if pval < 0.001:
        p_str = 'P < 0.001'
    elif pval < 0.01:
        p_str = f'P = {pval:.3f}'
    else:
        p_str = f'P = {pval:.2f}'

    ax.set_xlabel(f'{coord} eRNA (RPKM)', fontsize=9)
    ax.set_ylabel(f'{mrna_gene} mRNA (FPKM)', fontsize=9)
    ax.set_title(f'{coord} → {mrna_gene}\nρ={rho:.3f}, {p_str}', fontsize=10)
    ax.tick_params(labelsize=7)

plt.tight_layout()
fig.savefig(os.path.join(OUT_DIR, 'Figure_correlation.pdf'), bbox_inches='tight', dpi=300)
fig.savefig(os.path.join(OUT_DIR, 'Figure_correlation.png'), bbox_inches='tight', dpi=300)
print(f"\nSaved: Figure_correlation.pdf/png")

# === Figure 2: Correlation heatmap ===
fig2, ax2 = plt.subplots(1, 1, figsize=(8, 6))

# Build correlation matrix (eRNAs × unique genes)
erna_cols = list(erna_gene_map.keys())
gene_cols = list(set(erna_gene_map.values()))
corr_matrix = pd.DataFrame(index=[erna_info[c][0] for c in erna_cols],
                            columns=gene_cols, dtype=float)
pval_matrix = pd.DataFrame(index=[erna_info[c][0] for c in erna_cols],
                            columns=gene_cols, dtype=float)
for erna_col, mrna_gene in erna_gene_map.items():
    erna_label = erna_info[erna_col][0]
    rho_row = res_df[res_df['eRNA_column'] == erna_col]
    if not rho_row.empty:
        corr_matrix.loc[erna_label, mrna_gene] = rho_row['spearman_rho'].values[0]
        pval_matrix.loc[erna_label, mrna_gene] = rho_row['p_value'].values[0]

g = sns.heatmap(corr_matrix.astype(float), annot=True, fmt='.3f',
                cmap='RdBu_r', center=0, vmin=-1, vmax=1,
                linewidths=0.5, square=True,
                cbar_kws={'label': "Spearman's ρ"},
                ax=ax2)
ax2.set_title('eRNA-mRNA Spearman Correlation Matrix', fontsize=14)
ax2.set_xlabel('mRNA Gene', fontsize=11)
ax2.set_ylabel('eRNA Locus', fontsize=11)
plt.xticks(rotation=45, ha='right')
plt.yticks(rotation=0)
plt.tight_layout()
fig2.savefig(os.path.join(OUT_DIR, 'Figure_correlation_heatmap.pdf'), bbox_inches='tight', dpi=300)
fig2.savefig(os.path.join(OUT_DIR, 'Figure_correlation_heatmap.png'), bbox_inches='tight', dpi=300)
print(f"Saved: Figure_correlation_heatmap.pdf/png")
print("✓ Step 3 complete")