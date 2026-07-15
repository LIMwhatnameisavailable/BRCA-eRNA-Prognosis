#!/usr/bin/env python3
"""
Restyled v2 figures matching paper aesthetic:
  1. Fig_2.4_scatter_v2 — AR/ESR1 enrichment scatter with inset zoom
  2. Fig_2.4_heatmap_v2  — Subtype H3K27ac activity heatmap

Generates both as PDF (vector) + PNG 300dpi.
"""
import os, csv
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from matplotlib.colors import LinearSegmentedColormap
import matplotlib.gridspec as gridspec
from mpl_toolkits.axes_grid1.inset_locator import inset_axes, mark_inset

# ── Global style ──────────────────────────────────────────────────────────
plt.rcParams.update({
    'font.family': 'sans-serif',
    'font.size': 9,
    'axes.titlesize': 11,
    'axes.titleweight': 'bold',
    'axes.labelsize': 9,
    'xtick.labelsize': 8,
    'ytick.labelsize': 8,
    'axes.spines.top': False,
    'axes.spines.right': False,
    'axes.grid': False,
    'figure.facecolor': 'white',
    'axes.facecolor': 'white',
})

WORKDIR = ".."
ANALYSIS = os.path.join(WORKDIR, "analysis")

# ── Load unified summary data ─────────────────────────────────────────────
data = []
with open(os.path.join(ANALYSIS, "unified_summary.csv")) as f:
    reader = csv.DictReader(f)
    for row in reader:
        for k in ('ARE_hits', 'ERE_hits', 'FOXA1_hits', 'GATA3_hits',
                  'FOXA1_present', 'GATA3_present'):
            row[k] = int(row[k])
        for k in ('AR_enrichment', 'ER_enrichment',
                  'MCF7_H3K27ac', 'SKBR3_H3K27ac', 'MB453_H3K27ac',
                  'MB231_H3K27ac', 'Hs578T_H3K27ac'):
            row[k] = float(row[k])
        data.append(row)

# ── Subtype colors (paper palette) ───────────────────────────────────────
SUB_COLORS = {
    'Luminal-enriched': '#C0392B',
    'HER2-enriched': '#8E44AD',
    'LAR-enriched': '#C0392B',
    'TNBC-enriched': '#2980B9',
    'Pan-subtype': '#7F8C8D',
    'Pan-subtype (Luminal+LAR)': '#7F8C8D',
    'Inactive': '#7F8C8D',
}

MOTIF_MARKERS = {
    '✓ BOTH': 'o',
    'ARE only': '^',
    'ERE only': 's',
    'Neither': 'x',
}

MOTIF_COLORS_BAR = {
    '✓ BOTH': '#27AE60',
    'ARE only': '#E67E22',
    'ERE only': '#8E44AD',
    'Neither': '#BDC3C7',
}

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

SUB_RANK = {
    'Luminal-enriched': 1, 'HER2-enriched': 2, 'LAR-enriched': 2,
    'TNBC-enriched': 3, 'Pan-subtype': 4,
    'Pan-subtype (Luminal+LAR)': 4,
}

# ═══════════════════════════════════════════════════════════════════════════
# FIGURE 1: Scatter plot v2
# ═══════════════════════════════════════════════════════════════════════════
print("  Generating Fig_2.4_scatter_v2 ...")

fig, ax = plt.subplots(figsize=(5.5, 5.0))

# Plot all points
for d in data:
    ar, er = d['AR_enrichment'], d['ER_enrichment']
    color = SUB_COLORS.get(d['subtype_classification'], '#7F8C8D')
    marker = MOTIF_MARKERS.get(d['motif_cooccurrence'], 'o')
    size = 80 if d['motif_cooccurrence'] == 'Neither' else 120
    ec = 'black' if marker != 'x' else 'none'
    ax.scatter(ar, er, c=color, marker=marker, s=size,
               edgecolors=ec, linewidths=0.6, zorder=5, alpha=0.85)

# Annotate only the two outliers on main axes (chr10 loci)
for d in data:
    if 'chr10_5528926' in d['eRNA'] or 'chr10_5531356' in d['eRNA']:
        ar, er = d['AR_enrichment'], d['ER_enrichment']
        label = LOCUS_LABELS.get(d['eRNA'], d['eRNA'].replace('eRNA_', ''))
        ax.annotate(label, (ar, er), textcoords="offset points",
                    xytext=(8, 8), fontsize=7, alpha=0.9, fontfamily='monospace')

# Dashed reference lines at x=1, y=1
ax.axhline(y=1.0, color='#AAAAAA', linewidth=0.8, linestyle='--', zorder=1)
ax.axvline(x=1.0, color='#AAAAAA', linewidth=0.8, linestyle='--', zorder=1)

# Axis limits
ax.set_xlim(0, 12)
ax.set_ylim(0, 12)
ax.set_xlabel('AR Enrichment (vs Input)')
ax.set_ylabel('ESR1 Enrichment (vs Input)')

# Inset zoom: upper right
ax_inset = inset_axes(ax, width="45%", height="45%",
                       bbox_to_anchor=(0.98, 0.98, 0.55, 0.55),
                       bbox_transform=ax.transAxes, loc='upper right')

# Dashed rectangle on main axes indicating zoom region
from matplotlib.patches import Rectangle
ax.add_patch(Rectangle((0.4, 0.4), 1.6, 1.6, fill=False,
                        edgecolor='#888888', linestyle='--', linewidth=0.8))

# Populate inset
for d in data:
    ar, er = d['AR_enrichment'], d['ER_enrichment']
    color = SUB_COLORS.get(d['subtype_classification'], '#7F8C8C')
    marker = MOTIF_MARKERS.get(d['motif_cooccurrence'], 'o')
    sz = 60 if d['motif_cooccurrence'] == 'Neither' else 80
    ec = 'black' if marker != 'x' else 'none'
    ax_inset.scatter(ar, er, c=color, marker=marker, s=sz,
                     edgecolors=ec, linewidths=0.5, zorder=5, alpha=0.85)

# Label ALL points in the inset
for d in data:
    ar, er = d['AR_enrichment'], d['ER_enrichment']
    label = LOCUS_LABELS.get(d['eRNA'], d['eRNA'].replace('eRNA_', ''))
    ax_inset.annotate(label, (ar, er), textcoords="offset points",
                      xytext=(5, 4), fontsize=6, alpha=0.85,
                      fontfamily='monospace')

ax_inset.set_xlim(0.4, 2.0)
ax_inset.set_ylim(0.4, 2.0)
ax_inset.axhline(y=1.0, color='#AAAAAA', linewidth=0.6, linestyle='--', alpha=0.5)
ax_inset.axvline(x=1.0, color='#AAAAAA', linewidth=0.6, linestyle='--', alpha=0.5)
ax_inset.set_xlabel('AR', fontsize=7)
ax_inset.set_ylabel('ER', fontsize=7)
ax_inset.tick_params(labelsize=6)

# Connect inset with mark_inset
mark_inset(ax, ax_inset, loc1=1, loc2=3, fc="none", ec="#888888",
           linestyle="--", linewidth=0.6, alpha=0.6)

# Legend 1: Subtype (colors) — upper left, no box
seen_subtypes = []
for d in data:
    s = d['subtype_classification']
    if s not in seen_subtypes and s in SUB_COLORS:
        seen_subtypes.append(s)
# Sort by rank
seen_subtypes.sort(key=lambda s: SUB_RANK.get(s, 99))
leg_elements_sub = []
for s in seen_subtypes:
    n = sum(1 for d in data if d['subtype_classification'] == s)
    leg_elements_sub.append(
        Patch(facecolor=SUB_COLORS[s], edgecolor='black',
              linewidth=0.5, label=f'{s} (n={n})'))
leg1 = ax.legend(handles=leg_elements_sub, title='Subtype',
                 loc='upper left', fontsize=6.5, title_fontsize=7.5,
                 frameon=False)
ax.add_artist(leg1)

# Legend 2: Motif (shapes) — also upper left, offset
leg_elements_motif = [
    plt.Line2D([0], [0], marker='o', color='w', markerfacecolor='#7F8C8D',
               markersize=7, label='ARE + ERE'),
    plt.Line2D([0], [0], marker='^', color='w', markerfacecolor='#7F8C8D',
               markersize=7, label='ARE only'),
    plt.Line2D([0], [0], marker='s', color='w', markerfacecolor='#7F8C8D',
               markersize=7, label='ERE only'),
    plt.Line2D([0], [0], marker='x', color='#7F8C8D', markersize=7, label='Neither'),
]
leg2 = ax.legend(handles=leg_elements_motif, title='Motif',
                 loc='upper left', bbox_to_anchor=(0.0, 0.60),
                 fontsize=6.5, title_fontsize=7.5, frameon=False)

plt.subplots_adjust(left=0.15, right=0.95, bottom=0.12, top=0.95)

for ext in ['pdf', 'png']:
    path = os.path.join(ANALYSIS, f'Fig_2.4_scatter_v2.{ext}')
    plt.savefig(path, dpi=300)
    print(f"  ✓ Saved {path} ({os.path.getsize(path)/1024:.0f} KB)")

plt.close()

# ═══════════════════════════════════════════════════════════════════════════
# FIGURE 2: Subtype H3K27ac heatmap v2
# ═══════════════════════════════════════════════════════════════════════════
print("\n  Generating Fig_2.4_heatmap_v2 ...")

h3k27_cols = ['MCF7_H3K27ac', 'SKBR3_H3K27ac', 'MB453_H3K27ac',
              'MB231_H3K27ac', 'Hs578T_H3K27ac']
h3_matrix = np.zeros((len(data), len(h3k27_cols)))
for i, d in enumerate(data):
    for j, col in enumerate(h3k27_cols):
        h3_matrix[i, j] = d[col]

# z-score per locus
h3_norm = np.zeros_like(h3_matrix)
for i in range(len(data)):
    row = h3_matrix[i, :]
    if np.std(row) > 0:
        h3_norm[i, :] = (row - np.mean(row)) / np.std(row)
    else:
        h3_norm[i, :] = row

# Sort by subtype rank, then by MCF7 z-score descending within group
def sort_key(i):
    d = data[i]
    rank = SUB_RANK.get(d['subtype_classification'], 99)
    mcf7_z = h3_norm[i, 0]  # MCF7 z-score
    return (rank, -mcf7_z)

sorted_idx = sorted(range(len(data)), key=sort_key)
h3_norm_sorted = h3_norm[sorted_idx, :]
sorted_labels = [LOCUS_LABELS[data[i]['eRNA']] for i in sorted_idx]
sorted_subtypes = [data[i]['subtype_classification'] for i in sorted_idx]
sorted_motif = [data[i]['motif_cooccurrence'] for i in sorted_idx]

cell_line_labels = ['MCF7\n(LumA)', 'SKBR3\n(HER2)',
                    'MB453\n(LAR)', 'MB231\n(TNBC)', 'Hs578T\n(TNBC)']

# Paper warm colormap
paper_cmap = LinearSegmentedColormap.from_list(
    'paper_warm', ['#FFFDE7', '#F5A623', '#C0392B', '#7B0000'], N=256)

fig = plt.figure(figsize=(5.0, 4.5))
gs = gridspec.GridSpec(2, 4, width_ratios=[0.05, 1, 0.05, 0.05],
                       wspace=0.04, hspace=0.3, height_ratios=[1, 0.04])
plt.subplots_adjust(left=0.28, right=0.78, bottom=0.10, top=0.94)

# Main heatmap
ax_hm = plt.subplot(gs[0, 1])
im_hm = ax_hm.imshow(h3_norm_sorted, aspect='auto', cmap=paper_cmap,
                      interpolation='nearest', vmin=-2, vmax=2)

ax_hm.set_yticks(range(len(data)))
ax_hm.set_yticklabels(sorted_labels, fontsize=8, fontfamily='monospace')
ax_hm.set_xticks(range(len(cell_line_labels)))
ax_hm.set_xticklabels(cell_line_labels, fontsize=8)
ax_hm.set_title('H3K27ac Activity (z-score)')

# Left color strip: Subtype
ax_sub = plt.subplot(gs[0, 0])
for i, s in enumerate(sorted_subtypes):
    ax_sub.barh(i, 1, color=SUB_COLORS.get(s, '#7F8C8D'),
                edgecolor='none', linewidth=0)
ax_sub.set_yticks(range(len(data)))
ax_sub.set_yticklabels([])
ax_sub.set_xticks([])
ax_sub.invert_yaxis()

# Right color strip: Motif
ax_mot = plt.subplot(gs[0, 3])
motif_colors_bar = {
    '✓ BOTH': '#27AE60', 'ARE only': '#E67E22',
    'ERE only': '#8E44AD', 'Neither': '#BDC3C7'
}
for i, m in enumerate(sorted_motif):
    ax_mot.barh(i, 1, color=motif_colors_bar.get(m, '#BDC3C7'),
                edgecolor='none', linewidth=0)
ax_mot.set_yticks(range(len(data)))
ax_mot.set_yticklabels([])
ax_mot.set_xticks([])
ax_mot.invert_yaxis()

# Motif legend (to the right of the strip)
legend_elements_subtype = [
    Patch(facecolor='#C0392B', edgecolor='none', label='Luminal-enriched'),
    Patch(facecolor='#2980B9', edgecolor='none', label='TNBC-enriched'),
    Patch(facecolor='#8E44AD', edgecolor='none', label='HER2-enriched'),
    Patch(facecolor='#7F8C8D', edgecolor='none', label='Pan-subtype'),
]
leg_subtype = ax_mot.legend(handles=legend_elements_subtype, loc='center left',
                            bbox_to_anchor=(1.02, 0.85), fontsize=7,
                            title='Subtype', title_fontsize=6.5,
                            frameon=False, labelspacing=1.2)
ax_mot.add_artist(leg_subtype)

legend_elements_motif2 = [
    Patch(facecolor='#27AE60', edgecolor='none', label='ARE+ERE'),
    Patch(facecolor='#E67E22', edgecolor='none', label='ARE only'),
    Patch(facecolor='#8E44AD', edgecolor='none', label='ERE only'),
    Patch(facecolor='#BDC3C7', edgecolor='none', label='Neither'),
]
leg_motif = ax_mot.legend(handles=legend_elements_motif2, loc='center left',
                          bbox_to_anchor=(1.02, 0.5), fontsize=6,
                          title='Motif', title_fontsize=6.5,
                          frameon=False, labelspacing=1.2)

# Thin horizontal colorbar below heatmap
cbar_ax = plt.subplot(gs[1, 1])
cbar = plt.colorbar(im_hm, cax=cbar_ax, orientation='horizontal')
cbar.set_label('H3K27ac Activity (z-score)', fontsize=8)
cbar.set_ticks([-2, 0, 2])
cbar.set_ticklabels([-2, 0, 2], fontsize=7)

for ext in ['pdf', 'png']:
    path = os.path.join(ANALYSIS, f'Fig_2.4_heatmap_v2.{ext}')
    plt.savefig(path, dpi=300)
    print(f"  ✓ Saved {path} ({os.path.getsize(path)/1024:.0f} KB)")

plt.close()

print("\n  Done.")