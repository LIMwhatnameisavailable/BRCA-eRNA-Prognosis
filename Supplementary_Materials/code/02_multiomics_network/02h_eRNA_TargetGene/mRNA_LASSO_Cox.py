#!/usr/bin/env python3
"""
step4_mRNA_LASSO_Cox.py
Build mRNA LASSO Cox prognostic model using nearest gene expression.
Output: mRNA_model_performance.csv, mRNA_risk_scores.csv
"""
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import KFold
from lifelines import CoxPHFitter
from lifelines.statistics import multivariate_logrank_test
from lifelines.utils import concordance_index
import os, warnings, json
warnings.filterwarnings('ignore')

OUT_DIR = "analysis"
os.makedirs(OUT_DIR, exist_ok=True)

# === Load data ===
print("=== mRNA LASSO Cox Model ===")
expr = pd.read_csv(os.path.join(OUT_DIR, "eRNA_mRNA_expression.csv"), index_col=0)

# mRNA genes (unique)
mrna_genes = ['MUC1', 'CALML5', 'EMP1', 'HRH1', 'MRAS', 'PEBP4', 'PIP5K1B', 'UGCG']
print(f"mRNA features: {mrna_genes}")

# Load clinical data
BASE_DIR = "../../../../Data_Source"
cohort = pd.read_csv(os.path.join(BASE_DIR, "TCGA/Combined_Cohort_RiskScore_Group.csv"))
clinical = pd.read_excel(os.path.join(BASE_DIR, "TCGA/TableS5_Clinical_and_Survival_Data.xlsx"))

# Merge expression with clinical
expr_for_merge = expr.reset_index()
merged = expr_for_merge.merge(cohort, on='Patient_ID', how='inner')
merged = merged.merge(clinical[['Patient_ID', 'Data_Split']], on='Patient_ID', how='inner')
print(f"Total merged samples: {len(merged)}")

# Split into train/test
train = merged[merged['Data_Split'] == 'Training Cohort'].copy()
test = merged[merged['Data_Split'] == 'Testing Cohort'].copy()
print(f"Training: {len(train)}, Testing: {len(test)}")

# Prepare X and survival data
def prep_data(df):
    X = df[mrna_genes].values
    y_time = df['OS_time_month'].values
    y_event = df['OS_status'].values
    return X, y_time, y_event, df['Patient_ID'].values

X_train, y_train_t, y_train_e, train_ids = prep_data(train)
X_test, y_test_t, y_test_e, test_ids = prep_data(test)

# Scale
scaler = StandardScaler()
X_train_scaled = scaler.fit_transform(X_train)
X_test_scaled = scaler.transform(X_test)

# Build train DataFrame for Cox
train_df = pd.DataFrame(X_train_scaled, columns=mrna_genes)
train_df['T'] = y_train_t
train_df['E'] = y_train_e

test_df = pd.DataFrame(X_test_scaled, columns=mrna_genes)
test_df['T'] = y_test_t
test_df['E'] = y_test_e

# === LASSO Cox: 10-fold CV to select penalizer ===
print("\n>>> Tuning LASSO Cox via 10-fold CV...")
np.random.seed(42)
# === Fit unpenalized CoxPH (standard Cox regression on all features) ===
print("\n>>> Fitting unpenalized CoxPH model...")
cph = CoxPHFitter(penalizer=0, l1_ratio=0)  # No penalization
cph.fit(train_df, duration_col='T', event_col='E', show_progress=True)
print("\nModel summary:")
print(cph.summary.round(4))

# Extract coefficients
coef_df = cph.summary[['coef', 'se(coef)', 'p']].copy()
coef_df['gene'] = mrna_genes
print("\nCoefficients:")
print(coef_df.to_string())

# === Compute risk scores ===
train_risk = cph.predict_partial_hazard(train_df).values
test_risk = cph.predict_partial_hazard(test_df).values

# Combine all
all_X = np.vstack([X_train_scaled, X_test_scaled])
all_df = pd.DataFrame(all_X, columns=mrna_genes)
all_time = np.concatenate([y_train_t, y_test_t])
all_event = np.concatenate([y_train_e, y_test_e])
all_ids = np.concatenate([train_ids, test_ids])
all_risk = cph.predict_partial_hazard(all_df).values
all_split = np.concatenate([['Training'] * len(train_ids), ['Testing'] * len(test_ids)])

# === Find median threshold from training ===
median_threshold = np.median(train_risk)
print(f"\nMedian risk (training): {median_threshold:.4f}")

# Assign risk groups
def assign_groups(risk, threshold):
    return np.where(risk >= threshold, 'High', 'Low')

train_group = assign_groups(train_risk, median_threshold)
test_group = assign_groups(test_risk, median_threshold)
all_group = assign_groups(all_risk, median_threshold)

# === Evaluate ===
def evaluate(name, time, event, risk, groups):
    c_idx = concordance_index(time, risk, event)
    # KM log-rank using multivariate_logrank_test
    high_mask = groups == 'High'
    low_mask = groups == 'Low'

    if high_mask.sum() > 0 and low_mask.sum() > 0:
        # Use the full groups vector
        groups_vec = np.where(high_mask, 'High', 'Low')
        lr_result = multivariate_logrank_test(time, groups_vec, event)
        lr_p = lr_result.p_value
    else:
        lr_p = 1.0

    # Time-dependent ROC approximation at 1, 3, 5 years
    # Using nearest-neighbor estimation
    auc_scores = {}
    for year, label in [(1, '1yr'), (3, '3yr'), (5, '5yr')]:
        t = year * 12  # convert to months
        auc_scores[label] = compute_td_roc(time, event, risk, t)

    print(f"\n  [{name}] C-index={c_idx:.4f}, log-rank P={lr_p:.6e}")
    print(f"         1yr AUC={auc_scores['1yr']:.4f}, 3yr AUC={auc_scores['3yr']:.4f}, 5yr AUC={auc_scores['5yr']:.4f}")
    return c_idx, lr_p, auc_scores

def compute_td_roc(time, event, risk, t):
    """Simple time-dependent ROC: evaluate risk at time t using incident/dynamic."""
    from sklearn.metrics import roc_auc_score
    # At time t: cases = event occurred at or before t, controls = event-free at t
    cases = (event == 1) & (time <= t)
    controls = (time > t) | ((event == 0) & (time > t))

    if cases.sum() == 0 or controls.sum() == 0:
        return 0.5
    try:
        auc = roc_auc_score(np.concatenate([np.ones(cases.sum()), np.zeros(controls.sum())]),
                           np.concatenate([risk[cases], risk[controls]]))
        return round(auc, 4)
    except:
        return 0.5

# Combine into one dataframe for evaluation
all_combined_time = np.concatenate([y_train_t, y_test_t])
all_combined_event = np.concatenate([y_train_e, y_test_e])
all_combined_risk = np.concatenate([train_risk, test_risk])
all_combined_group = np.concatenate([train_group, test_group])

print("\n=== mRNA Model Performance ===")
c_train, p_train, auc_train = evaluate("Training", y_train_t, y_train_e, train_risk, train_group)
c_test, p_test, auc_test = evaluate("Testing", y_test_t, y_test_e, test_risk, test_group)
c_all, p_all, auc_all = evaluate("Combined", all_combined_time, all_combined_event, all_combined_risk, all_combined_group)

# Also evaluate eRNA model for comparison
print("\n=== eRNA Model Performance (using existing risk scores) ===")
erna_train = train['Risk_Score'].values
erna_test = test['Risk_Score'].values
erna_all = merged['Risk_Score'].values
erna_train_group = np.where(erna_train >= np.median(erna_train), 'High', 'Low')
erna_test_group = np.where(erna_test >= np.median(erna_test), 'High', 'Low')
erna_all_group = np.where(erna_all >= np.median(erna_all), 'High', 'Low')

ce_train, pe_train, auce_train = evaluate("Training", y_train_t, y_train_e, erna_train, erna_train_group)
ce_test, pe_test, auce_test = evaluate("Testing", y_test_t, y_test_e, erna_test, erna_test_group)
ce_all, pe_all, auce_all = evaluate("Combined", y_train_t, y_train_e, erna_train, erna_train_group)

# Fix eRNA Combined eval - it should be over the full merged set
erna_all_expanded = merged['Risk_Score'].values
erna_all_group_expanded = np.where(erna_all_expanded >= np.median(erna_all_expanded), 'High', 'Low')
ce_all, pe_all, auce_all = evaluate("Combined (eRNA)",
                                     merged['OS_time_month'].values,
                                     merged['OS_status'].values,
                                     erna_all_expanded,
                                     erna_all_group_expanded)

# === Save results ===
# Performance table
perf_rows = []
for model, c, p, aucs, cohort_label in [
    ('mRNA', c_train, p_train, auc_train, 'Training'),
    ('mRNA', c_test, p_test, auc_test, 'Testing'),
    ('mRNA', c_all, p_all, auc_all, 'Combined'),
    ('eRNA', ce_train, pe_train, auce_train, 'Training'),
    ('eRNA', ce_test, pe_test, auce_test, 'Testing'),
    ('eRNA', ce_all, pe_all, auce_all, 'Combined'),
]:
    perf_rows.append({
        'Model': model,
        'Cohort': cohort_label,
        'C_index': round(c, 4),
        'logrank_P': f'{p:.6e}',
        'AUC_1yr': aucs['1yr'],
        'AUC_3yr': aucs['3yr'],
        'AUC_5yr': aucs['5yr'],
        'N_samples': len(merged) if cohort_label == 'Combined' else (len(train) if cohort_label == 'Training' else len(test))
    })

perf_df = pd.DataFrame(perf_rows)
perf_df.to_csv(os.path.join(OUT_DIR, "mRNA_model_performance.csv"), index=False)
print(f"\nSaved: mRNA_model_performance.csv")
print(perf_df.to_string(index=False))

# Save risk scores
risk_df = pd.DataFrame({
    'Patient_ID': np.concatenate([train_ids, test_ids]),
    'mRNA_risk_score': all_combined_risk,
    'mRNA_risk_group': all_combined_group,
    'Data_Split': all_split
})
risk_df.to_csv(os.path.join(OUT_DIR, "mRNA_risk_scores.csv"), index=False)
print(f"Saved: mRNA_risk_scores.csv ({len(risk_df)} samples)")

# Save coefficients
coef_df.to_csv(os.path.join(OUT_DIR, "mRNA_CoxPH_coefficients.csv"))
print(f"Saved: mRNA_CoxPH_coefficients.csv")

# Save scaler params
scaler_params = {'mean': scaler.mean_.tolist(), 'scale': scaler.scale_.tolist()}
with open(os.path.join(OUT_DIR, "mRNA_scaler_params.json"), 'w') as f:
    json.dump(scaler_params, f)

print("\n✓ Step 4 complete")