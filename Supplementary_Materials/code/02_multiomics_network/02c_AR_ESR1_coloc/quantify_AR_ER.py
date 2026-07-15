#!/usr/bin/env python3
"""
quantify_AR_ER.py — Quantitative analysis of AR/ER co-binding at prognostic eRNA loci.

Reads AR_ER_per_eRNA_counts.tab from multiBigwigSummary,
computes enrichment ratios, identifies co-bound loci, and generates a scatter plot.
"""

import os
import csv
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

WORKDIR = ".."
INFILE = os.path.join(WORKDIR, "analysis", "AR_ER_per_eRNA_counts.tab")
OUT_CSV = os.path.join(WORKDIR, "analysis", "AR_ER_enrichment_table.csv")
OUT_FIG = os.path.join(WORKDIR, "analysis", "Fig_2.4_AR_ER_scatter.pdf")

# ── 1. Read data ──────────────────────────────────────────────────────────
samples = []
with open(INFILE) as f:
    header = f.readline().strip()
    col_names = [c.strip("'") for c in header.split('\t')]
    # columns: chr, start, end, DMSO_AR.bw, DMSO_ER.bw, DMSO_Input.bw
    for line in f:
        parts = line.strip().split('\t')
        if len(parts) < 6:
            continue
        chrom = parts[0].strip("'")
        start = int(parts[1])
        end = int(parts[2])
        ar_signal = float(parts[3])
        er_signal = float(parts[4])
        input_signal = float(parts[5])
        # Use the eRNA name from the iGenomics-style naming
        samples.append({
            'chrom': chrom,
            'start': start,
            'end': end,
            'ar': ar_signal,
            'er': er_signal,
            'input': input_signal,
        })

# ── 2. Compute enrichment and determine co-binding ────────────────────────
results = []
for s in samples:
    ar_enrich = s['ar'] / s['input'] if s['input'] > 0 else float('inf')
    er_enrich = s['er'] / s['input'] if s['input'] > 0 else float('inf')
    co_bound = ar_enrich > 1.5 and er_enrich > 1.5
    results.append({
        'eRNA': f"eRNA_{s['chrom']}_{(s['start'] + s['end']) // 2}",
        'AR_signal': round(s['ar'], 4),
        'ER_signal': round(s['er'], 4),
        'Input_signal': round(s['input'], 4),
        'AR_enrichment': round(ar_enrich, 3),
        'ER_enrichment': round(er_enrich, 3),
        'co_bound': co_bound,
    })

n_co = sum(1 for r in results if r['co_bound'])
n_total = len(results)

# ── 3. Write CSV table ────────────────────────────────────────────────────
with open(OUT_CSV, 'w', newline='') as f:
    writer = csv.DictWriter(f, fieldnames=[
        'eRNA', 'AR_signal', 'ER_signal', 'Input_signal',
        'AR_enrichment', 'ER_enrichment', 'co_bound',
    ])
    writer.writeheader()
    writer.writerows(results)

print(f"  ✓ Wrote {OUT_CSV}")
print(f"    Total eRNA loci: {n_total}")
print(f"    AR-enriched (>1.5× Input): {sum(1 for r in results if r['AR_enrichment'] > 1.5)}")
print(f"    ER-enriched (>1.5× Input): {sum(1 for r in results if r['ER_enrichment'] > 1.5)}")
print(f"    Co-bound loci: {n_co}")

# ── 4. Scatter plot: AR enrichment vs ER enrichment ───────────────────────
fig, ax = plt.subplots(figsize=(8, 7))

for r in results:
    color = 'red' if r['co_bound'] else 'steelblue'
    ax.scatter(r['AR_enrichment'], r['ER_enrichment'],
               c=color, s=80, edgecolors='black', linewidths=0.5, zorder=5)
    ax.annotate(r['eRNA'].replace('eRNA_', ''),
                (r['AR_enrichment'], r['ER_enrichment']),
                textcoords="offset points", xytext=(6, 6),
                fontsize=7, alpha=0.8)

# Reference lines at 1.5× enrichment
ax.axhline(y=1.5, color='gray', linestyle='--', linewidth=0.8, alpha=0.6)
ax.axvline(x=1.5, color='gray', linestyle='--', linewidth=0.8, alpha=0.6)

# Labels and title
ax.set_xlabel('AR Enrichment (vs Input)', fontsize=12)
ax.set_ylabel('ESR1 Enrichment (vs Input)', fontsize=12)
ax.set_title('AR vs ESR1 Binding Enrichment at Prognostic eRNA Loci', fontsize=13, fontweight='bold')

# Legend
from matplotlib.patches import Patch
legend_elements = [
    Patch(facecolor='red', edgecolor='black', label=f'Co-bound (n={n_co})'),
    Patch(facecolor='steelblue', edgecolor='black', label=f'Single/None (n={n_total - n_co})'),
]
ax.legend(handles=legend_elements, loc='lower right', fontsize=10)

# Axis limits with padding
all_ar = [r['AR_enrichment'] for r in results]
all_er = [r['ER_enrichment'] for r in results]
max_val = max(max(all_ar), max(all_er)) * 1.2
min_val = 0
ax.set_xlim(min_val, max_val)
ax.set_ylim(min_val, max_val)

# Diagonal line y=x for reference
ax.plot([0, max_val], [0, max_val], color='gray', linestyle=':', linewidth=0.5, alpha=0.4)

ax.grid(alpha=0.3)

plt.tight_layout()
plt.savefig(OUT_FIG, dpi=200)
print(f"  ✓ Saved {OUT_FIG}")

# ── 5. Summary for stdout ─────────────────────────────────────────────────
print()
print("=" * 60)
print("  AR-ESR1 Co-binding Summary")
print("=" * 60)
print(f"  {'eRNA':<30} {'AR/Input':>10} {'ER/Input':>10} {'Co-bound':>10}")
print(f"  {'-'*30} {'-'*10} {'-'*10} {'-'*10}")
for r in results:
    cb = "✓ YES" if r['co_bound'] else "  —"
    print(f"  {r['eRNA']:<30} {r['AR_enrichment']:>10.3f} {r['ER_enrichment']:>10.3f} {cb:>10}")
print(f"  {'-'*60}")
print(f"  Total co-bound: {n_co}/{n_total} loci")
print(f"  Mean AR enrichment: {sum(r['AR_enrichment'] for r in results)/n_total:.3f}")
print(f"  Mean ER enrichment: {sum(r['ER_enrichment'] for r in results)/n_total:.3f}")
print()