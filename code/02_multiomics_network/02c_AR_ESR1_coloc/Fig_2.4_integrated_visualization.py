#!/usr/bin/env python3
"""
Phase 5: Integrated Visualization for AR-ESR1 Motif Co-occurrence Analysis.
Generates 3 publication-quality figures:
  1. Updated scatter plot (AR enrichment vs ER enrichment) with motif + subtype overlay
  2. Motif co-occurrence panel (all 4 TFs from real HOMER scans)
  3. Subtype H3K27ac activity heatmap

All motif data (ARE, ERE, FOXA1, GATA3) generated via HOMER scanMotifGenomeWide.pl
using hg38-lifted eRNA sequences. Locus labels use hg19 identifiers (matching Table 2).
"""
import os, csv
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from matplotlib.colors import Normalize, LinearSegmentedColormap
import matplotlib.gridspec as gridspec
from mpl_toolkits.axes_grid1.inset_locator import inset_axes, mark_inset

WORKDIR = ".."
ANALYSIS = os.path.join(WORKDIR, "analysis")

# ── Load unified summary data (includes FOXA1 + GATA3 from real HOMER scans) ─
data = []
with open(os.path.join(ANALYSIS, "unified_summary.csv")) as f:
    reader = csv.DictReader(f)
    for row in reader:
        # Convert numeric fields
        for k in ('ARE_hits', 'ERE_hits', 'FOXA1_hits', 'GATA3_hits',
                  'FOXA1_present', 'GATA3_present'):
            row[k] = int(row[k])
        for k in ('AR_enrichment', 'ER_enrichment',
                  'MCF7_H3K27ac', 'SKBR3_H3K27ac', 'MB453_H3K27ac',
                  'MB231_H3K27ac', 'Hs578T_H3K27ac'):
            row[k] = float(row[k])
        data.append(row)

# ── Define color/style maps ───────────────────────────────────────────────
SUB_COLORS = {
    'Luminal-enriched': '#E74C3C',     # warm red
    'HER2-enriched': '#9B59B6',         # purple
    'LAR-enriched': '#E67E22',          # orange
    'TNBC-enriched': '#3498DB',         # steel blue
    'Pan-subtype': '#95A5A6',           # gray
    'Pan-subtype (Luminal+LAR)': '#F39C12',  # gold
    'Inactive': '#95A5A6',              # gray
}

MOTIF_MARKERS = {
    '✓ BOTH': 'o',      # circle
    'ARE only': '^',     # triangle up
    'ERE only': 's',     # square
    'Neither': 'X',      # X
}

MOTIF_COLORS_BAR = {
    '✓ BOTH': '#27AE60',    # green
    'ARE only': '#E67E22',  # orange
    'ERE only': '#8E44AD',  # purple
    'Neither': '#BDC3C7',   # gray
}

# Short locus labels for figures — these are hg19 identifiers matching Table 2
LOCUS_LABELS = {
    'eRNA_chr1_155158995': 'chr1:155158995',
    'eRNA_chr10_5528926': 'chr10:5528926',
    'eRNA_chr10_5531356': 'chr10:5531356',
    'eRNA_chr12_13371038': 'chr12:13371038',
    'eRNA_chr3_11236700': 'chr3:11236700',
    'eRNA_chr3_138070534': 'chr3:138070534',
    'eRNA_chr8_22624675': 'chr8:22624675',
    'eRNA_chr9_71398719': 'chr9:71398719',
    'eRNA_chr9_71398939': 'chr9:71398939',
    'eRNA_chr9_114689796': 'chr9:114689796',
}

# Subtype sort order for heatmap rows
SUB_RANK = {'Luminal-enriched': 0, 'HER2-enriched': 1, 'LAR-enriched': 2,
            'TNBC-enriched': 3, 'Pan-subtype': 4,
            'Pan-subtype (Luminal+LAR)': 5}

# ═══════════════════════════════════════════════════════════════════════════
# FIGURE 1: Updated Scatter Plot — AR enrichment vs ER enrichment
# ═══════════════════════════════════════════════════════════════════════════
print("  Generating Figure 1: Updated scatter plot...")

fig, ax = plt.subplots(figsize=(7.5, 7))

for d in data:
    ar = d['AR_enrichment']
    er = d['ER_enrichment']
    locus = d['eRNA']
    subtype = d['subtype_classification']
    motif = d['motif_cooccurrence']

    color = SUB_COLORS.get(subtype, '#95A5A6')
    marker = MOTIF_MARKERS.get(motif, 'o')

    ax.scatter(ar, er, c=color, marker=marker, s=180, edgecolors='black',
               linewidths=1.0, zorder=5, alpha=0.85)

    # Label each point (hg19 identifiers)
    label = LOCUS_LABELS.get(locus, locus.replace('eRNA_', ''))
    ax.annotate(label, (ar, er), textcoords="offset points",
                xytext=(8, 8), fontsize=7.5, alpha=0.9, fontfamily='monospace')

# Reference lines
ax.axhline(y=1.5, color='gray', linestyle='--', linewidth=0.8, alpha=0.5)
ax.axvline(x=1.5, color='gray', linestyle='--', linewidth=0.8, alpha=0.5)

# Diagonal line
max_val = max(max(d['AR_enrichment'] for d in data),
              max(d['ER_enrichment'] for d in data)) * 1.15
ax.plot([0, max_val], [0, max_val], color='gray', linestyle=':',
        linewidth=0.5, alpha=0.3)

ax.set_xlabel('AR Enrichment (vs Input)', fontsize=12, fontweight='bold')
ax.set_ylabel('ESR1 Enrichment (vs Input)', fontsize=12, fontweight='bold')
ax.set_title('AR vs ESR1: Motif Co-occurrence & Subtype Activity\n'
             'at Prognostic eRNA Loci',
             fontsize=12, fontweight='bold')
ax.set_xlim(0, max_val)
ax.set_ylim(0, max_val)

# Journal style: remove top/right spines
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

# Legend 1: Subtype (color)
legend_elements_sub = []
for sub, col in SUB_COLORS.items():
    n = sum(1 for d in data if d['subtype_classification'] == sub)
    if n > 0:
        legend_elements_sub.append(
            Patch(facecolor=col, edgecolor='black', label=f'{sub} (n={n})'))

leg1 = ax.legend(handles=legend_elements_sub, title='Subtype',
                 loc='upper left', fontsize=7.5, title_fontsize=8.5,
                 framealpha=0.9)

# Legend 2: Motif (marker shape)
legend_elements_motif = [
    plt.Line2D([0], [0], marker='o', color='w', markerfacecolor='gray',
               markersize=10, label='ARE + ERE (co-occur)'),
    plt.Line2D([0], [0], marker='^', color='w', markerfacecolor='gray',
               markersize=10, label='ARE only'),
    plt.Line2D([0], [0], marker='s', color='w', markerfacecolor='gray',
               markersize=10, label='ERE only'),
    plt.Line2D([0], [0], marker='X', color='w', markerfacecolor='gray',
               markersize=10, label='Neither'),
]
leg2 = ax.legend(handles=legend_elements_motif, title='Motif',
                 loc='lower right', fontsize=7.5, title_fontsize=8.5,
                 framealpha=0.9)
ax.add_artist(leg1)  # Keep both legends

# Zoomed inset for lower-left cluster (AR < 2, ER < 2)
# Find points in the crowded region
cluster_points = [d for d in data if d['AR_enrichment'] < 2.5
                  and d['ER_enrichment'] < 2.5]
if len(cluster_points) > 3:
    ax_inset = inset_axes(ax, width="40%", height="40%", loc='center right',
                          bbox_to_anchor=(0.15, 0.15, 0.7, 0.7),
                          bbox_transform=ax.transAxes)
    for d in cluster_points:
        ar, er = d['AR_enrichment'], d['ER_enrichment']
        color = SUB_COLORS.get(d['subtype_classification'], '#95A5A6')
        marker = MOTIF_MARKERS.get(d['motif_cooccurrence'], 'o')
        ax_inset.scatter(ar, er, c=color, marker=marker, s=120,
                         edgecolors='black', linewidths=0.8, zorder=5, alpha=0.9)
        label = LOCUS_LABELS.get(d['eRNA'], d['eRNA'].replace('eRNA_', ''))
        ax_inset.annotate(label, (ar, er), textcoords="offset points",
                          xytext=(5, 5), fontsize=6.5, fontfamily='monospace')

    ax_inset.set_xlim(0.3, 2.0)
    ax_inset.set_ylim(0.3, 2.0)
    ax_inset.axhline(y=1.5, color='gray', linestyle='--', linewidth=0.6, alpha=0.4)
    ax_inset.axvline(x=1.5, color='gray', linestyle='--', linewidth=0.6, alpha=0.4)
    ax_inset.set_xlabel('AR', fontsize=7)
    ax_inset.set_ylabel('ER', fontsize=7)
    ax_inset.tick_params(labelsize=6)
    ax_inset.spines['top'].set_visible(False)
    ax_inset.spines['right'].set_visible(False)

    # Connect inset to main plot
    mark_inset(ax, ax_inset, loc1=1, loc2=3, fc="none", ec="gray", linestyle="--", alpha=0.5)

    # Label pointing to the inset region
    ax.annotate('see inset', xy=(0.5, 0.5), xytext=(0.7, 0.35),
                ha='center', fontsize=7, color='gray', fontstyle='italic',
                arrowprops=dict(arrowstyle='->', color='gray', lw=0.8, alpha=0.6))

ax.grid(alpha=0.15)
plt.tight_layout()

for ext in ['pdf', 'png']:
    path = os.path.join(ANALYSIS, f'Fig_2.4_scatter_updated.{ext}')
    plt.savefig(path, dpi=300)
    print(f"  ✓ Saved {path}")

plt.close()

# ═══════════════════════════════════════════════════════════════════════════
# FIGURE 2: Motif Co-occurrence Panel — ALL 4 TFs from real HOMER scans
# ═══════════════════════════════════════════════════════════════════════════
print("\n  Generating Figure 2: Motif co-occurrence panel...")

# Build motif matrix from REAL scan data: 10 loci × 4 motifs
motif_names = ['ARE', 'ERE', 'FOXA1', 'GATA3']
motif_matrix = np.zeros((len(data), len(motif_names)))

for i, d in enumerate(data):
    motif_matrix[i, 0] = 1 if d['ARE_hits'] > 0 else 0      # ARE — real
    motif_matrix[i, 1] = 1 if d['ERE_hits'] > 0 else 0      # ERE — real
    motif_matrix[i, 2] = 1 if d['FOXA1_hits'] > 0 else 0    # FOXA1 — real scan
    motif_matrix[i, 3] = 1 if d['GATA3_hits'] > 0 else 0    # GATA3 — real scan

# ChIP-seq enrichment for sidebars
chip_ar = np.array([d['AR_enrichment'] for d in data])
chip_er = np.array([d['ER_enrichment'] for d in data])

# Build the figure
fig = plt.figure(figsize=(10, 7))
gs = gridspec.GridSpec(2, 3, width_ratios=[1, 0.12, 0.12], height_ratios=[1, 0.05],
                       hspace=0.3, wspace=0.25)

# Panel A: Motif presence heatmap
ax1 = plt.subplot(gs[0, 0])
cmap_motif = LinearSegmentedColormap.from_list('motif', ['#f0f0f0', '#2C3E50'], N=2)
im1 = ax1.imshow(motif_matrix, aspect='auto', cmap=cmap_motif, interpolation='nearest')

ax1.set_yticks(range(len(data)))
ax1.set_yticklabels([LOCUS_LABELS[d['eRNA']] for d in data], fontsize=8, fontfamily='monospace')
ax1.set_xticks(range(len(motif_names)))
ax1.set_xticklabels(motif_names, fontsize=9, fontweight='bold')
ax1.set_title('TF Motif Presence (HOMER scans)', fontsize=11, fontweight='bold')

# Add checkmark annotations
for i in range(len(data)):
    for j in range(len(motif_names)):
        val = '✓' if motif_matrix[i, j] > 0 else ''
        ax1.text(j, i, val, ha='center', va='center',
                color='white' if motif_matrix[i, j] > 0 else '#999',
                fontsize=8, fontweight='bold')

# Panel B: AR ChIP-seq sidebar (warm red)
ax2 = plt.subplot(gs[0, 1])
ar_max = max(max(chip_ar), 1.5)
colors_ar = plt.cm.Reds(chip_ar / ar_max)
ax2.barh(range(len(data)), chip_ar, color=colors_ar, edgecolor='gray', linewidth=0.5)
ax2.axvline(x=1.5, color='darkred', linestyle='--', linewidth=0.8, alpha=0.5)
ax2.set_yticks(range(len(data)))
ax2.set_yticklabels([])
ax2.invert_yaxis()
ax2.set_title('AR\nChIP-seq', fontsize=9, fontweight='bold')
ax2.set_xlabel('Enrich.', fontsize=8)

# Panel C: ER ChIP-seq sidebar (warm red)
ax3 = plt.subplot(gs[0, 2])
er_max = max(max(chip_er), 1.5)
colors_er = plt.cm.Reds(chip_er / er_max)
ax3.barh(range(len(data)), chip_er, color=colors_er, edgecolor='gray', linewidth=0.5)
ax3.axvline(x=1.5, color='darkred', linestyle='--', linewidth=0.8, alpha=0.5)
ax3.set_yticks(range(len(data)))
ax3.set_yticklabels([])
ax3.invert_yaxis()
ax3.set_title('ER\nChIP-seq', fontsize=9, fontweight='bold')
ax3.set_xlabel('Enrich.', fontsize=8)

# Colorbar
cbar_ax = plt.subplot(gs[1, 0])
cbar = plt.colorbar(im1, cax=cbar_ax, orientation='horizontal')
cbar.set_ticks([0, 1])
cbar.set_ticklabels(['Absent', 'Present'])
cbar.set_label('Motif Status', fontsize=9)

plt.suptitle('TF Motif Co-occurrence at Prognostic eRNA Loci\n'
             '(4 TFs, HOMER scanMotifGenomeWide.pl)',
             fontsize=13, fontweight='bold', y=1.02)

for ext in ['pdf', 'png']:
    path = os.path.join(ANALYSIS, f'Fig_2.4_motif_panel.{ext}')
    plt.savefig(path, dpi=300, bbox_inches='tight')
    print(f"  ✓ Saved {path}")

plt.close()

# ═══════════════════════════════════════════════════════════════════════════
# FIGURE 3: Subtype H3K27ac Activity Heatmap
# ═══════════════════════════════════════════════════════════════════════════
print("\n  Generating Figure 3: Subtype activity heatmap...")

# H3K27ac values
h3k27_cols = ['MCF7_H3K27ac', 'SKBR3_H3K27ac', 'MB453_H3K27ac',
              'MB231_H3K27ac', 'Hs578T_H3K27ac']
h3_matrix = np.zeros((len(data), len(h3k27_cols)))
for i, d in enumerate(data):
    for j, col in enumerate(h3k27_cols):
        h3_matrix[i, j] = d[col]

# Row-normalize (z-score per locus)
h3_norm = np.zeros_like(h3_matrix)
for i in range(len(data)):
    row = h3_matrix[i, :]
    if np.std(row) > 0:
        h3_norm[i, :] = (row - np.mean(row)) / np.std(row)
    else:
        h3_norm[i, :] = row

# Sort by subtype
sorted_idx = sorted(range(len(data)),
                    key=lambda i: SUB_RANK.get(data[i]['subtype_classification'], 99))
h3_norm_sorted = h3_norm[sorted_idx, :]
sorted_labels = [LOCUS_LABELS[data[i]['eRNA']] for i in sorted_idx]
sorted_subtypes = [data[i]['subtype_classification'] for i in sorted_idx]
sorted_motif = [data[i]['motif_cooccurrence'] for i in sorted_idx]

cell_line_labels = ['MCF7\n(LumA)', 'SKBR3\n(HER2)',
                    'MB453\n(LAR)', 'MB231\n(TNBC)', 'Hs578T\n(TNBC)']

fig = plt.figure(figsize=(10, 7))
gs = gridspec.GridSpec(2, 4, width_ratios=[0.08, 1, 0.06, 0.06],
                       wspace=0.08, hspace=0.3, height_ratios=[1, 0.04])
# Increase left margin so Y-axis locus labels are fully visible
plt.subplots_adjust(left=0.22)

# Main heatmap — warm colormap (white → dark red)
ax_hm = plt.subplot(gs[0, 1])
warm_cmap = LinearSegmentedColormap.from_list('warm_red',
                                              ['#FFFFFF', '#FEE5D9', '#FCAE91',
                                               '#FB6A4A', '#DE2D26', '#A50F15'], N=256)
im_hm = ax_hm.imshow(h3_norm_sorted, aspect='auto', cmap=warm_cmap,
                      interpolation='nearest', vmin=-2, vmax=2)

ax_hm.set_yticks(range(len(data)))
ax_hm.set_yticklabels(sorted_labels, fontsize=8, fontfamily='monospace')
ax_hm.set_xticks(range(len(cell_line_labels)))
ax_hm.set_xticklabels(cell_line_labels, fontsize=7.5)
ax_hm.set_title('H3K27ac Activity (z-score)', fontsize=11, fontweight='bold')

# Subtype annotation bar (left)
ax_sub = plt.subplot(gs[0, 0])
for i, s in enumerate(sorted_subtypes):
    ax_sub.barh(i, 1, color=SUB_COLORS.get(s, '#95A5A6'),
                edgecolor='gray', linewidth=0.5)
ax_sub.set_yticks(range(len(data)))
ax_sub.set_yticklabels([])
ax_sub.set_xticks([])
ax_sub.invert_yaxis()
ax_sub.set_title('Subtype', fontsize=9, fontweight='bold')

# Motif co-occurrence annotation bar (right)
ax_mot = plt.subplot(gs[0, 3])
motif_colors_bar = {
    '✓ BOTH': '#27AE60', 'ARE only': '#E67E22',
    'ERE only': '#8E44AD', 'Neither': '#BDC3C7'
}
for i, m in enumerate(sorted_motif):
    ax_mot.barh(i, 1, color=motif_colors_bar.get(m, '#BDC3C7'),
                edgecolor='gray', linewidth=0.5, alpha=0.8)
ax_mot.set_yticks(range(len(data)))
ax_mot.set_yticklabels([])
ax_mot.set_xticks([])
ax_mot.invert_yaxis()
ax_mot.set_title('Motif', fontsize=9, fontweight='bold')

# Thin horizontal colorbar
cbar_ax2 = plt.subplot(gs[1, 1])
cbar2 = plt.colorbar(im_hm, cax=cbar_ax2, orientation='horizontal')
cbar2.set_label('z-score (per locus)', fontsize=9)

# Motif color legend
legend_elements_motif2 = [
    Patch(facecolor='#27AE60', edgecolor='black', label='ARE + ERE'),
    Patch(facecolor='#E67E22', edgecolor='black', label='ARE only'),
    Patch(facecolor='#8E44AD', edgecolor='black', label='ERE only'),
    Patch(facecolor='#BDC3C7', edgecolor='black', label='Neither'),
]
leg_motif = ax_mot.legend(handles=legend_elements_motif2, loc='center left',
                          bbox_to_anchor=(1.0, 0.5), fontsize=7,
                          title='Motif type', title_fontsize=8,
                          framealpha=0.9)

plt.suptitle('Breast Cancer Subtype H3K27ac Activity at Prognostic eRNA Loci',
             fontsize=13, fontweight='bold', y=1.02)

for ext in ['pdf', 'png']:
    path = os.path.join(ANALYSIS, f'Fig_2.4_subtype_heatmap.{ext}')
    plt.savefig(path, dpi=300, bbox_inches='tight')
    print(f"  ✓ Saved {path}")

plt.close()

# ═══════════════════════════════════════════════════════════════════════════
print("\n" + "="*60)
print("  ALL FIGURES GENERATED")
print("="*60)
print("\n  Notes:")
print("  - FOXA1 and GATA3 motifs from real HOMER scanMotifGenomeWide.pl runs")
print("  - Locus labels are hg19 identifiers (matching manuscript Table 2)")
print("  - All analyses use hg38-aligned ChIP-seq / H3K27ac data")
print()