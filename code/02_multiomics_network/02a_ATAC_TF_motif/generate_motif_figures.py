#!/usr/bin/env python3
"""
Generate ATAC TF motif dot-plot figure (GSEA style) for manuscript revision.

Integrates MCF7 (GSE298769) and TCGA BRCA (Corces et al., 2018) ATAC-seq
HOMER motif results into a two-panel publication-quality figure.

Output:
  Fig_2.1_ATAC_motif.pdf — two-panel dot plot (MCF7 | TCGA)
  Fig_2.1_ATAC_motif.png — same at 300 dpi
"""

import os, re, math
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.legend_handler import HandlerPatch
import numpy as np

# =====================================================================
# Paths
# =====================================================================
BASE = ".."
MCF7_ALL = os.path.join(BASE, "analysis", "homer_ATAC_all", "knownResults.txt")
TCGA_ALL = os.path.join(BASE, "analysis", "homer_TCGA_all", "knownResults.txt")
OUT_PDF  = os.path.join(BASE, "analysis", "Fig_2.1_ATAC_motif.pdf")
OUT_PNG  = os.path.join(BASE, "analysis", "Fig_2.1_ATAC_motif.png")

# =====================================================================
# Key TFs to highlight with ★ (breast cancer TF network + AP-1)
# =====================================================================
KEY_TF_PREFIXES = [
    "CTCF", "CTCFL", "BORIS",
    "FOXA1", "FOXA2", "FOXA3", "FOXM1", "FOX",
    "GATA3", "GATA",
    "ESR1", "ERa", "ERα", "ERb", "ERβ", "ERE", "Erra", "ERRg",
    "AR", "ARE",
    "Fos", "Jun", "AP-1", "Fra", "Atf", "BATF", "Fosl",
]


def is_key_tf(tf_name):
    """Check if a TF name matches any key TF prefix."""
    upper = tf_name.upper()
    for prefix in KEY_TF_PREFIXES:
        if upper.startswith(prefix.upper()):
            return True
    return False


# =====================================================================
# Motif parsing
# =====================================================================
def parse_motif_name(name):
    """Extract short TF name from HOMER motif name.
    "FOXA1(Forkhead)/MCF7-FOXA1-ChIP-Seq(...)" -> "FOXA1"
    """
    m = re.match(r'^([^(]+)\(', name)
    if m:
        return m.group(1).strip()
    return name.split("/")[0].strip()


def read_known_results(filepath, n_top=15):
    """Read HOMER knownResults.txt, return top N motifs by -log10(p-value)."""
    motifs = []
    with open(filepath) as f:
        header = f.readline().strip().split('\t')

        def find_col(prefix):
            for i, h in enumerate(header):
                if h.startswith(prefix):
                    return i
            raise ValueError(f"Column prefix '{prefix}' not found")

        name_idx = find_col("Motif Name")
        # Use "Log P-value" as the primary significance column
        logp_idx = find_col("Log P-value")
        target_pct_idx = find_col("% of Target Sequences with Motif")
        bg_pct_idx = find_col("% of Background Sequences with Motif")
        target_count_idx = find_col("# of Target Sequences with Motif")

        for line in f:
            parts = line.strip().split('\t')
            if len(parts) <= max(name_idx, logp_idx):
                continue
            try:
                logp = float(parts[logp_idx])
                target_pct = float(parts[target_pct_idx].rstrip('%'))
                bg_pct = float(parts[bg_pct_idx].rstrip('%'))
                n_target = float(parts[target_count_idx])
            except (ValueError, IndexError):
                continue

            if n_target <= 0:
                continue

            # -log10(p-value) from Log P-value (which is natural log)
            neg_log10_pval = -logp / math.log(10) if logp < 0 else 0
            enrichment = target_pct / max(bg_pct, 0.001)
            tf_short = parse_motif_name(parts[name_idx])

            motifs.append({
                'name': parts[name_idx],
                'tf': tf_short,
                'neg_log10_pval': neg_log10_pval,
                'target_pct': target_pct,
                'bg_pct': bg_pct,
                'enrichment': enrichment,
                'n_target': n_target,
                'is_key': is_key_tf(tf_short),
            })

    # Sort by -log10(p-value) descending, take top N
    motifs.sort(key=lambda x: x['neg_log10_pval'], reverse=True)
    motifs = motifs[:n_top]

    # Deduplicate labels: "FOXA1", "FOXA1(2)", "FOXA1(3)", ...
    seen = {}
    for m in motifs:
        base = m['tf']
        if base in seen:
            seen[base] += 1
            m['label'] = f"{base}({seen[base]})"
        else:
            seen[base] = 1
            m['label'] = base

    return motifs


# =====================================================================
# Plot a single dot-plot panel (GSEA style)
# =====================================================================
def plot_panel(ax, motifs, title, n_total_seq, cmap):
    """Single GSEA-style dot-plot panel.

    Color = -log10(p-value)  (red = high, blue = low)
    Size  = % target sequences
    """
    # Motifs are sorted most-significant-first; place at top of y-axis
    y_pos = list(range(len(motifs)))[::-1]  # reverse: rank 0 (best) at top
    x_vals = [m['neg_log10_pval'] for m in motifs]
    pct_vals = [m['target_pct'] for m in motifs]
    is_key = [m['is_key'] for m in motifs]
    labels = [m['label'] for m in motifs]

    # Bubble sizes: scale proportional to % target, clip for visibility
    sizes = [max(v * 30, 50) for v in pct_vals]

    # Colors from colormap
    max_logp = max(x_vals) if x_vals else 1
    norm_vals = [v / max_logp for v in x_vals]
    colors = cmap(norm_vals)

    # Draw scatter
    ax.scatter(x_vals, y_pos, s=sizes, c=colors,
               edgecolors='gray', linewidth=0.4, alpha=0.85,
               zorder=3)

    # ★ overlay for key TFs
    for i, (x, y, k) in enumerate(zip(x_vals, y_pos, is_key)):
        if k:
            ax.annotate('★', (x, y), fontsize=12, color='red',
                        ha='center', va='center', fontweight='bold',
                        zorder=5)

    # Y-axis: motif names, most significant at top
    ax.set_yticks(y_pos)
    ax.set_yticklabels(labels, fontsize=8, style='italic')

    # X-axis
    ax.set_xlabel('−log₁₀(p-value)', fontsize=9)

    # Title
    ax.set_title(title, fontsize=11, fontweight='bold', pad=6, loc='left')

    # Clean style: no grid, left+bottom spines only
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.tick_params(axis='both', labelsize=8)
    ax.set_axisbelow(True)

    # Sample size annotation
    ax.text(0.98, 0.02, f'n = {n_total_seq:,} peaks',
            transform=ax.transAxes, fontsize=7.5, color='#666666',
            ha='right', va='bottom', style='italic',
            zorder=6)

    # ---- Colorbar (outside panel, to the right) ----
    cb_ax = ax.inset_axes([1.02, 0.10, 0.030, 0.65])
    norm_obj = plt.Normalize(vmin=0, vmax=max_logp)
    cb = plt.colorbar(
        plt.cm.ScalarMappable(norm=norm_obj, cmap=cmap),
        cax=cb_ax, orientation='vertical'
    )
    cb.set_label('−log₁₀(p-value)', fontsize=7, labelpad=1)
    cb.ax.tick_params(labelsize=6)
    cb.outline.set_linewidth(0.3)

    # ---- Size legend (inside panel, bottom-left area) ----
    _add_size_legend(ax, pct_vals)


def _add_size_legend(ax, pct_vals):
    """Add a small bubble-size legend inside the panel."""
    # Pick 3 representative sizes from the data range
    max_pct = max(pct_vals) if pct_vals else 20
    if max_pct < 5:
        sizes_pct = [2, 5, 10]
    elif max_pct < 15:
        sizes_pct = [5, 10, 15]
    else:
        sizes_pct = [10, 20, 30]

    # Clamp to data range
    sizes_pct = [s for s in sizes_pct if s <= max_pct * 1.2]
    if not sizes_pct:
        sizes_pct = [round(max_pct * 0.3), round(max_pct * 0.6), round(max_pct * 0.9)]
    if len(sizes_pct) < 2:
        sizes_pct = [sizes_pct[0], sizes_pct[0] * 2, sizes_pct[0] * 3]

    # Position legend in bottom-left corner of the axes
    legend_x = 0.52  # fraction of axes
    legend_y = 0.10

    for i, pct in enumerate(sizes_pct):
        r = max(pct * 30, 50) ** 0.5  # radius in points
        # Convert point radius to axes coordinates (approximate)
        # Place bubbles horizontally
        ax.scatter([], [], s=max(pct * 30, 50),
                   color='gray', edgecolors='gray',
                   linewidth=0.3, alpha=0.4, label=f'{pct}%')

    # Create legend for sizes
    leg = ax.legend(
        title='% target',
        title_fontsize=6,
        fontsize=5.5,
        loc='lower left',
        framealpha=0.8,
        edgecolor='#cccccc',
        handlelength=0.5,
        borderpad=0.3,
        labelspacing=1.0,
        markerscale=0.5,
    )
    leg.get_title().set_fontweight('bold')


# =====================================================================
# Extract target count from header
# =====================================================================
def get_target_count(filepath):
    """Parse the target sequence count from the HOMER header line."""
    with open(filepath) as f:
        header = f.readline()
    m = re.search(r'\(of\s+([\d,]+)\)', header)
    return int(m.group(1).replace(',', '')) if m else 0


# =====================================================================
# Main
# =====================================================================
def main():
    print("=" * 60)
    print("ATAC TF motif dot-plot (GSEA style)")
    print("=" * 60)

    # ---- Read data ----
    print("\nReading HOMER knownResults.txt files...")
    mcf7_motifs = read_known_results(MCF7_ALL, 15)
    tcga_motifs = read_known_results(TCGA_ALL, 15)

    mcf7_n = get_target_count(MCF7_ALL)
    tcga_n = get_target_count(TCGA_ALL)

    print(f"\nMCF7 all-peaks: {len(mcf7_motifs)} motifs, {mcf7_n:,} target peaks")
    print(f"  Top 5 motifs (label | −log10(p) | %target):")
    for m in mcf7_motifs[:5]:
        print(f"    {m['label']:12s}  −log10(p)={m['neg_log10_pval']:8.1f}  %target={m['target_pct']:5.2f}%")

    print(f"\nTCGA all-peaks: {len(tcga_motifs)} motifs, {tcga_n:,} target peaks")
    print(f"  Top 5 motifs (label | −log10(p) | %target):")
    for m in tcga_motifs[:5]:
        print(f"    {m['label']:12s}  −log10(p)={m['neg_log10_pval']:8.1f}  %target={m['target_pct']:5.2f}%")

    # ---- Create figure ----
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13, 6), facecolor='white')

    # Colormap: RdYlBu_r = red=high, blue=low (like GSEA)
    cmap = plt.colormaps['RdYlBu_r']

    plot_panel(ax1, mcf7_motifs, 'A  MCF7 ATAC-seq (GSE298769)', mcf7_n, cmap)
    plot_panel(ax2, tcga_motifs, 'B  TCGA BRCA ATAC-seq (Corces et al., 2018)', tcga_n, cmap)

    # Global footer
    fig.text(0.5, 0.01,
             '★ = Key TF network motif (CTCF, FOXA1, GATA3, AP-1, ESR1/ERα, AR)',
             ha='center', fontsize=7.5, color='gray', style='italic')

    plt.subplots_adjust(wspace=0.55, left=0.07, right=0.90, top=0.92, bottom=0.10)

    # ---- Save ----
    fig.savefig(OUT_PDF, dpi=300, bbox_inches='tight')
    fig.savefig(OUT_PNG, dpi=300, bbox_inches='tight')
    plt.close()

    for fpath in [OUT_PDF, OUT_PNG]:
        sz = os.path.getsize(fpath)
        print(f"\n  Saved: {fpath} ({sz/1024:.1f} KB)")

    print("\nDone.")


if __name__ == '__main__':
    main()
