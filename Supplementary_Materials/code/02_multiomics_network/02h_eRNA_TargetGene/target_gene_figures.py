#!/usr/bin/env python3
"""
step5_figures_km_roc.py
Generate mRNA model KM curve and ROC AUC plots.
Output: Figure_model_comparison.{pdf,png}
"""
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from lifelines import KaplanMeierFitter
from lifelines.statistics import multivariate_logrank_test
from sklearn.metrics import roc_auc_score
import os, warnings
warnings.filterwarnings('ignore')

OUT_DIR = "analysis"
COLOR_HIGH = '#E87A5D'  # warm orange
COLOR_LOW = '#5B7B9A'   # blue-grey

# === Load data ===
merged = pd.read_csv(os.path.join(OUT_DIR, "eRNA_mRNA_expression.csv"), index_col=0)
BASE_DIR = "../../../../Data_Source"
cohort = pd.read_csv(os.path.join(BASE_DIR, "TCGA/Combined_Cohort_RiskScore_Group.csv"))
clinical = pd.read_excel(os.path.join(BASE_DIR, "TCGA/TableS5_Clinical_and_Survival_Data.xlsx"))

# mRNA risk scores
risk_df = pd.read_csv(os.path.join(OUT_DIR, "mRNA_risk_scores.csv"))

# Merge
merged_for_merge = merged.reset_index()
merged_data = merged_for_merge.merge(cohort, on='Patient_ID', how='inner')
merged_data = merged_data.merge(clinical[['Patient_ID', 'Data_Split']], on='Patient_ID', how='inner')

# Add risk scores
merged_data = merged_data.merge(risk_df[['Patient_ID', 'mRNA_risk_score', 'mRNA_risk_group']],
                                on='Patient_ID', how='inner')
print(f"Total: {len(merged_data)}")

# Get mRNA risk group assignments
all_time = merged_data['OS_time_month'].values
all_event = merged_data['OS_status'].values
all_risk = merged_data['mRNA_risk_score'].values
all_group = merged_data['mRNA_risk_group'].values

# For training, testing and combined
train_data = merged_data[merged_data['Data_Split'] == 'Training Cohort']
test_data = merged_data[merged_data['Data_Split'] == 'Testing Cohort']

datasets = {
    'Combined Cohort (N=1072)': (all_time, all_event, all_risk, all_group),
}

# === Figure 1: KM Curves (Combined Cohort) ===
print("\n>>> Generating KM curve...")
fig, axes = plt.subplots(2, 2, figsize=(14, 10))
fig.suptitle('mRNA Model: Kaplan-Meier Survival Curves', fontsize=16, y=1.01)

for idx, (label, (t, e, risk, group)) in enumerate(datasets.items()):
    row, col = idx // 2, idx % 2
    ax = axes[row, col]

    high_mask = group == 'High'
    low_mask = group == 'Low'

    kmf_high = KaplanMeierFitter()
    kmf_low = KaplanMeierFitter()

    kmf_high.fit(t[high_mask], e[high_mask], label=f'High Risk (N={high_mask.sum()})')
    kmf_low.fit(t[low_mask], e[low_mask], label=f'Low Risk (N={low_mask.sum()})')

    kmf_high.plot_survival_function(ax=ax, color=COLOR_HIGH, linewidth=2, ci_show=True)
    kmf_low.plot_survival_function(ax=ax, color=COLOR_LOW, linewidth=2, ci_show=True)

    # Log-rank test
    groups_vec = np.where(high_mask, 'High', 'Low')
    lr_result = multivariate_logrank_test(t, groups_vec, e)
    lr_p = lr_result.p_value
    p_text = f'P = {lr_p:.2e}' if lr_p < 0.001 else f'P = {lr_p:.4f}'

    ax.set_xlabel('Time (months)', fontsize=12)
    ax.set_ylabel('Survival Probability', fontsize=12)
    ax.set_title(f'{label}', fontsize=13)
    ax.text(0.6, 0.15, p_text, transform=ax.transAxes, fontsize=12,
            bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
    ax.legend(fontsize=10)
    ax.set_ylim(0, 1)
    ax.grid(alpha=0.3)

# Add Training + Testing KM curves
for idx, (dataset_label, data) in enumerate([
    ('Training Cohort (N=752)', train_data),
    ('Testing Cohort (N=320)', test_data)
]):
    row = 1 if idx == 0 else row
    col = 0 if idx == 0 else 1
    ax = axes[row, col]

    t = data['OS_time_month'].values
    e = data['OS_status'].values
    group = data['mRNA_risk_group'].values

    high_mask = group == 'High'
    low_mask = group == 'Low'

    kmf_high = KaplanMeierFitter()
    kmf_low = KaplanMeierFitter()
    kmf_high.fit(t[high_mask], e[high_mask], label=f'High Risk (N={high_mask.sum()})')
    kmf_low.fit(t[low_mask], e[low_mask], label=f'Low Risk (N={low_mask.sum()})')
    kmf_high.plot_survival_function(ax=ax, color=COLOR_HIGH, linewidth=2, ci_show=True)
    kmf_low.plot_survival_function(ax=ax, color=COLOR_LOW, linewidth=2, ci_show=True)

    groups_vec = np.where(high_mask, 'High', 'Low')
    lr_result = multivariate_logrank_test(t, groups_vec, e)
    lr_p = lr_result.p_value
    p_text = f'P = {lr_p:.2e}' if lr_p < 0.001 else f'P = {lr_p:.4f}'

    ax.set_xlabel('Time (months)', fontsize=12)
    ax.set_ylabel('Survival Probability', fontsize=12)
    ax.set_title(f'{dataset_label}', fontsize=13)
    ax.text(0.6, 0.15, p_text, transform=ax.transAxes, fontsize=12,
            bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
    ax.legend(fontsize=10)
    ax.set_ylim(0, 1)
    ax.grid(alpha=0.3)

plt.tight_layout()
fig.savefig(os.path.join(OUT_DIR, 'Figure_model_comparison.pdf'), bbox_inches='tight', dpi=300)
fig.savefig(os.path.join(OUT_DIR, 'Figure_model_comparison.png'), bbox_inches='tight', dpi=300)
print("Saved: Figure_model_comparison.pdf/png")

# === Figure 2: ROC AUC comparison bar chart ===
print("\n>>> Generating ROC AUC bar chart...")
perf = pd.read_csv(os.path.join(OUT_DIR, "mRNA_model_performance.csv"))

# Focus on mRNA model
mrna_perf = perf[perf['Model'] == 'mRNA']
erna_perf = perf[perf['Model'] == 'eRNA']

fig2, axes2 = plt.subplots(1, 2, figsize=(12, 5))

# Left: C-index comparison (mRNA only)
ax = axes2[0]
cohorts = mrna_perf['Cohort'].tolist()
c_indices = mrna_perf['C_index'].tolist()
bars = ax.bar(cohorts, c_indices, color=[COLOR_HIGH, COLOR_LOW, '#7A9B7A'])
ax.axhline(0.5, color='gray', linestyle='--', alpha=0.5, label='Random (C=0.5)')
for bar, val in zip(bars, c_indices):
    ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 0.01,
            f'{val:.3f}', ha='center', fontsize=10)
ax.set_ylabel('C-index (Concordance)', fontsize=12)
ax.set_title('mRNA Model: Concordance Index', fontsize=13)
ax.set_ylim(0, 0.8)
ax.legend(fontsize=10)
ax.grid(alpha=0.3, axis='y')

# Right: Time-dependent ROC AUC (3 bars per cohort)
ax = axes2[1]
x = np.arange(len(mrna_perf))
width = 0.25
for i, (year, color) in enumerate([('1yr', '#E8B57D'), ('3yr', '#D4956A'), ('5yr', '#C0755A')]):
    vals = mrna_perf[f'AUC_{year}'].values
    offset = (i - 1) * width
    bars = ax.bar(x + offset, vals, width, label=f'{year.replace("yr", "-year")}', color=color)
    for bar, val in zip(bars, vals):
        ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 0.01,
                f'{val:.3f}', ha='center', fontsize=9)

ax.axhline(0.5, color='gray', linestyle='--', alpha=0.5, label='Random (0.5)')
ax.set_xticks(x)
ax.set_xticklabels(mrna_perf['Cohort'].tolist(), fontsize=11)
ax.set_ylabel('AUC (Time-dependent ROC)', fontsize=12)
ax.set_title('mRNA Model: Time-dependent ROC AUC', fontsize=13)
ax.set_ylim(0, 1.0)
ax.legend(fontsize=9)
ax.grid(alpha=0.3, axis='y')

plt.tight_layout()
fig2.savefig(os.path.join(OUT_DIR, 'Figure_ROC_comparison.pdf'), bbox_inches='tight', dpi=300)
fig2.savefig(os.path.join(OUT_DIR, 'Figure_ROC_comparison.png'), bbox_inches='tight', dpi=300)
print("Saved: Figure_ROC_comparison.pdf/png")

print("\n✓ Step 5 complete")
print(f"\nAll output files in {OUT_DIR}/")
import os as _os
for f in sorted(_os.listdir(OUT_DIR)):
    if f.endswith('.csv') or f.endswith('.pdf') or f.endswith('.png'):
        print(f"  - {f} ({_os.path.getsize(_os.path.join(OUT_DIR, f))/1024:.1f} KB)")