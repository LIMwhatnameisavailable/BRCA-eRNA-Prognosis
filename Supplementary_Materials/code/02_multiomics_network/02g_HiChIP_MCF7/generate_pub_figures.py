#!/usr/bin/env python3
"""
Generate publication-quality HiChIP browser plots for Figure 3E.
All coordinates in hg38. Replaces polyA RNA-seq with PRO-seq for eRNA detection.
Fixes the genome version mismatch (loops were in hg19, other tracks in hg38).
"""

import os, sys, gzip, math
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.lines import Line2D

try:
    import pyBigWig
    HAS_BIGWIG = True
except ImportError:
    print("ERROR: pyBigWig required but not installed.")
    sys.exit(1)

# ── Paths ──────────────────────────────────────────────────────────────
BASE = ".."
DIR = "../02g_HiChIP_MCF7/analysis"

# bigWig tracks (all hg38)
BW_H3K27AC = "../02b_ChIPseq_TF_tracks/GSE298767_Veh_H3K27ac_merge.bw"
BW_PRO_PLUS = "../02d_RNAseq_eRNA_signal/GSE298770_MCF7_Veh_plus.bw"
BW_PRO_MINUS = "../02d_RNAseq_eRNA_signal/GSE298770_MCF7_Veh_minus.bw"
BW_ER = "../02b_ChIPseq_TF_tracks/GSE298767_ER_Veh.bw"
BW_FOXA1 = "../02b_ChIPseq_TF_tracks/GSE298767_FOXA1_veh_R1_R2.bw"
BW_GATA3 = "../02b_ChIPseq_TF_tracks/GSE298767_GATA3_veh_R1_R2.bw"
BW_ATAC = "../02a_ATAC_TF_motif/data/GSE298769_Veh_merged.bw"

OUT_PREFIX = f"{DIR}/Fig_3.7_HiChIP_pub"

DPI = 300

# ── Gene annotation data (from UCSC refGene hg38) ─────────────────────
# chr10 region: 5,400,000-5,700,000
CHR10_GENES = [
    ("NET1",  "chr10", 5412556, 5459056, "+"),
    ("CALML5","chr10", 5498696, 5499570, "-"),
    ("CALML3","chr10", 5524960, 5526771, "+"),
    ("LASTR", "chr10", 5594984, 5596118, "-"),
    ("ASB13", "chr10", 5638866, 5666595, "-"),
    ("TASOR2","chr10", 5684837, 5763740, "+"),
]

# chr8 region: 22,600,000-22,800,000
CHR8_GENES = [
    ("CCAR2", "chr8", 22604756, 22621514, "+"),
    ("BIN3",  "chr8", 22620430, 22669121, "-"),
    ("EGR3",  "chr8", 22687658, 22693480, "-"),
    ("PEBP4", "chr8", 22713250, 22927914, "-"),
]

def read_track(bw_path, chrom, start, end, max_points=2000):
    """Read bigWig values over a 1D interval, with max-pool downsampling.

    Downsampling to ~max_points bins prevents matplotlib from creating
    PDF paths with hundreds of thousands of vertices, which causes some
    PDF renderers (Adobe Illustrator, Edge) to clip or discard content.
    """
    if not os.path.exists(bw_path):
        return None, None
    bw = pyBigWig.open(bw_path)
    try:
        if chrom not in bw.chroms():
            print(f"  WARNING: {os.path.basename(bw_path)} missing chrom {chrom}")
            return None, None
        vals = bw.values(chrom, start, end)
        vals = np.nan_to_num(vals, nan=0.0)
        positions = np.arange(start, end)

        n = len(vals)
        if n > max_points:
            # Max-pool: divide into max_points bins, take the max value per bin
            # This preserves peak shapes far better than mean-pooling.
            bin_edges = np.linspace(0, n, max_points + 1).astype(int)
            bin_centers = np.empty(max_points)
            downsampled = np.empty(max_points)
            for i in range(max_points):
                lo, hi = bin_edges[i], bin_edges[i + 1]
                downsampled[i] = np.max(vals[lo:hi])
                bin_centers[i] = (positions[lo] + positions[hi - 1]) / 2.0
            return bin_centers, downsampled

        return positions, vals
    finally:
        bw.close()

def draw_gene_track(ax, genes, chrom, start, end, max_y=2.0):
    """Draw gene annotations as arrows."""
    ax.set_xlim(start, end)
    ax.set_ylim(-0.3, max_y + 1.0)
    ax.set_yticks([])
    ax.spines["left"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["top"].set_visible(False)

    y_levels = {}
    next_y = 0.5

    for gname, gchrom, gs, ge, gstrand in genes:
        if gchrom != chrom:
            continue
        if ge < start or gs > end + 20000:  # Allow slight overflow
            continue

        if gname not in y_levels:
            y_levels[gname] = next_y
            next_y += 1.0
            if next_y > max_y:
                next_y = 0.5

        y = y_levels[gname]
        vis_s = max(gs, start)
        vis_e = min(ge, end)
        width = vis_e - vis_s

        # Draw gene body as a line
        ax.plot([vis_s, vis_e], [y, y], color="#333333", linewidth=2.0, zorder=3)

        # Draw arrowhead for direction
        arrow_len = min(width * 0.3, 3000)
        if gstrand == "+":
            ax.annotate("", xy=(vis_e, y), xytext=(vis_e - arrow_len, y),
                        arrowprops=dict(arrowstyle="->", color="#333333", lw=1.5))
        else:
            ax.annotate("", xy=(vis_s, y), xytext=(vis_s + arrow_len, y),
                        arrowprops=dict(arrowstyle="->", color="#333333", lw=1.5))

        # Label
        ax.text((vis_s + vis_e) / 2, y + 0.35, gname, ha="center", va="bottom",
                fontsize=7, fontweight="bold", color="#222222",
                bbox=dict(boxstyle="round,pad=0.1", fc="white", ec="none", alpha=0.7))


def draw_eRNA_track(ax, eRNA_bed, chrom, start, end):
    """Draw eRNA loci."""
    ax.set_xlim(start, end)
    ax.set_ylim(-0.3, 1.8)
    ax.set_yticks([])
    ax.spines["left"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["top"].set_visible(False)

    for ename, echrom, es, ee in eRNA_bed:
        if echrom != chrom:
            continue
        if ee < start or es > end:
            continue
        vis_s = max(es, start)
        vis_e = min(ee, end)
        ax.barh(0.7, vis_e - vis_s, left=vis_s, height=0.8,
                color="#E67E22", edgecolor="#D35400", linewidth=0.8, zorder=3)
        ax.text((vis_s + vis_e) / 2, 0.7, ename, ha="center", va="center",
                fontsize=6.5, fontweight="bold", color="white")


def draw_hic_arcs(ax, loops, chrom, start, end, max_span=200000):
    """Draw HiChIP loops as arcs."""
    ax.set_xlim(start, end)
    ax.set_ylim(0, 1.8)
    ax.set_yticks([])
    ax.spines["left"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["top"].set_visible(False)

    n_drawn = 0
    for loop in loops:
        lchrom, l1s, l1e, rchrom, r1s, r1e = loop[:6]
        cc = int(loop[6]) if len(loop) > 6 else 10

        if lchrom != chrom:
            continue
        if max(l1s, r1s) < start or min(l1e, r1e) > end + 10000:
            continue

        left = min(l1s, r1s)
        right = max(l1e, r1e)
        if right <= left:
            continue

        span = right - left
        rel_h = min(span / max_span * 0.7, 0.85)
        yc = 0.5 + rel_h * 0.3
        width_span = span

        lw = min(max(cc / 4, 1.5), 5.0)
        arc = mpatches.Arc(xy=((left + right) / 2, yc),
                           width=width_span,
                           height=rel_h * 0.7,
                           angle=180, theta1=0, theta2=180,
                           color="#C0392B", linewidth=lw, alpha=0.85)
        ax.add_patch(arc)
        # Add contact count label
        ax.text((left + right) / 2, yc + rel_h * 0.35,
                f"CC={cc}", ha="center", va="bottom",
                fontsize=6, color="#C0392B", fontweight="bold")
        n_drawn += 1

    if n_drawn == 0:
        ax.text((start + end) / 2, 0.9, "No loops in view",
                ha="center", va="center", fontsize=9, color="gray", style="italic")


def make_browser_plot(chrom, center, half, loops, eRNA_bed, genes,
                      track_configs, out_pdf, title="", condition=""):
    """Generate multi-track browser plot."""
    start = center - half
    end = center + half

    n_tracks = len(track_configs) + 3  # signal tracks + HiChIP + eRNA + genes

    height_ratios = []
    for tc in track_configs:
        height_ratios.append(tc.get("height", 2.5))
    height_ratios += [1.8, 0.8, 2.0]  # HiChIP, eRNA, genes

    fig_height = 1.5 + sum(height_ratios) * 0.45

    fig, axes = plt.subplots(len(height_ratios), 1,
                             figsize=(11, fig_height),
                             gridspec_kw={"height_ratios": height_ratios,
                                          "hspace": 0.08},
                             sharex=True)
    ax_list = axes if hasattr(axes, "__iter__") else [axes]
    ax_list = list(ax_list)

    tidx = 0

    # ── 1. Signal tracks (bigWigs) ──
    for tc in track_configs:
        ax = ax_list[tidx]
        pos, vals = read_track(tc["path"], chrom, start, end)
        if pos is not None and vals is not None:
            if tc.get("neg_strand", False):
                # Flip minus strand values to positive for display
                vals = -vals
            ax.fill_between(pos, vals, alpha=0.7, color=tc["color"],
                            linewidth=0.2)
            ax.set_ylabel(tc["label"], fontsize=7.5)
            ax.set_ylim(bottom=0)
            if "ymax" in tc:
                ax.set_ylim(top=tc["ymax"])
        else:
            ax.text((start + end) / 2, 0.5, f"No data\n({tc['label']})",
                    ha="center", va="center", fontsize=7, color="gray")
            ax.set_ylabel(tc["label"], fontsize=7.5)
        ax.tick_params(labelsize=6)
        ax.yaxis.set_major_locator(plt.MaxNLocator(4))
        tidx += 1

    # ── 2. HiChIP arcs ──
    ax_arc = ax_list[tidx]
    draw_hic_arcs(ax_arc, loops, chrom, start, end)
    ax_arc.set_ylabel(f"HiChIP\n{condition}", fontsize=7.5)
    tidx += 1

    # ── 3. eRNA track ──
    ax_erna = ax_list[tidx]
    draw_eRNA_track(ax_erna, eRNA_bed, chrom, start, end)
    ax_erna.set_ylabel("eRNA", fontsize=7.5)
    tidx += 1

    # ── 4. Gene annotations ──
    ax_genes = ax_list[tidx]
    draw_gene_track(ax_genes, genes, chrom, start, end)
    ax_genes.set_ylabel("Genes", fontsize=7.5)
    ax_genes.set_xlabel(f"{chrom} position (Mb)", fontsize=8)
    ax_genes.xaxis.set_major_formatter(
        plt.FuncFormatter(lambda x, p: f"{x / 1e6:.2f}"))

    # ── 5. Title ──
    if title:
        fig.suptitle(title, fontsize=10, fontweight="bold", y=0.99)

    plt.tight_layout()
    plt.savefig(out_pdf, bbox_inches="tight", pad_inches=0.3)
    plt.close()
    print(f"  ✅ Saved: {out_pdf}")


# ════════════════════════════════════════════════════════════════════════
# Generate Chr10 DMSO plot (Figure 3E main panel)
# ════════════════════════════════════════════════════════════════════════
print("=" * 65)
print("  Figure 3E — Chr10 DMSO HiChIP Browser Plot (hg38 corrected)")
print("=" * 65)

CHR10_CENTER = 5560000   # hg38: center of 5.46-5.66 Mb window (shows ASB13 at right edge)
CHR10_HALF = 100000       # ±100kb = 5,460,000-5,660,000

# DMSO loop (hg38 coordinates from liftover)
dmso_loops_hg38 = [
    ("chr10", 5488037, 5493037, "chr10", 5618037, 5623037, 12),
]

# eRNA loci in hg38
erna_bed_chr10 = [
    ("eRNA #1", "chr10", 5486463, 5487463),
    ("eRNA #2", "chr10", 5488893, 5489893),
]

# chr10 genes
chr10_genes = [
    ("NET1",  "chr10", 5412556, 5459056, "+"),
    ("CALML5","chr10", 5498696, 5499570, "-"),
    ("CALML3","chr10", 5524960, 5526771, "+"),
    ("LASTR", "chr10", 5594984, 5596118, "-"),
    ("ASB13", "chr10", 5638866, 5666595, "-"),
]

# Track configs for chr10: H3K27ac, PRO-seq (+), PRO-seq (-), FOXA1
chr10_tracks = [
    {"path": BW_H3K27AC, "color": "#1B7837", "label": "H3K27ac",
     "height": 2.5, "ymax": 120},
    {"path": BW_PRO_PLUS, "color": "#2166AC", "label": "PRO-seq (+)",
     "height": 2.5, "ymax": 15},
    {"path": BW_PRO_MINUS, "color": "#B2182B", "label": "PRO-seq (−)",
     "height": 2.5, "ymax": 15, "neg_strand": True},
    {"path": BW_FOXA1, "color": "#762A83", "label": "FOXA1",
     "height": 2.0, "ymax": 30},
]

make_browser_plot("chr10", CHR10_CENTER, CHR10_HALF,
                  dmso_loops_hg38, erna_bed_chr10, chr10_genes,
                  chr10_tracks,
                  f"{OUT_PREFIX}_chr10_DMSO_v2.pdf",
                  title="Figure 3E — DMSO: eRNA chr10:5,528,926-5,531,356 → ASB13 (128 kb loop)",
                  condition="DMSO")

# ════════════════════════════════════════════════════════════════════════
# Generate Chr10 LY plot (Supplementary)
# ════════════════════════════════════════════════════════════════════════
print("\n" + "=" * 65)
print("  Supplementary — Chr10 LY HiChIP Browser Plot (hg38 corrected)")
print("=" * 65)

# LY loops (hg38): anchor1~5,483,537-5,488,537, anchor2~5,543,537-5,548,537
ly_loops_chr10 = [
    ("chr10", 5483037, 5488037, "chr10", 5543037, 5548037, 18),
    ("chr10", 5488037, 5493037, "chr10", 5543037, 5548037, 28),
]

ly_tracks = [
    {"path": BW_H3K27AC, "color": "#1B7837", "label": "H3K27ac",
     "height": 2.5, "ymax": 120},
    {"path": BW_PRO_PLUS, "color": "#2166AC", "label": "PRO-seq (+)",
     "height": 2.5, "ymax": 15},
    {"path": BW_PRO_MINUS, "color": "#B2182B", "label": "PRO-seq (−)",
     "height": 2.5, "ymax": 15, "neg_strand": True},
]

make_browser_plot("chr10", CHR10_CENTER, 100000,
                  ly_loops_chr10, erna_bed_chr10, chr10_genes,
                  ly_tracks,
                  f"{OUT_PREFIX}_chr10_LY_v2.pdf",
                  title="Supplementary — LY: Chr10 eRNA HiChIP (tighter loops, higher contact counts)",
                  condition="LY")

# ════════════════════════════════════════════════════════════════════════
# Generate Chr8 LY plot (Supplementary)
# ════════════════════════════════════════════════════════════════════════
print("\n" + "=" * 65)
print("  Supplementary — Chr8 LY HiChIP (eRNA → EGR3)")
print("=" * 65)

CHR8_CENTER = 22720000  # Between eRNA and EGR3
CHR8_HALF = 100000

# LY loop chr8: left anchor=22,692,487-22,697,487, right=22,767,487-22,772,487
ly_loops_chr8 = [
    ("chr8", 22692487, 22697487, "chr8", 22767487, 22772487, 14),
]

erna_bed_chr8 = [
    ("eRNA", "chr8", 22766662, 22767662),
]

chr8_genes = [
    ("CCAR2", "chr8", 22604756, 22621514, "+"),
    ("BIN3",  "chr8", 22620430, 22669121, "-"),
    ("EGR3",  "chr8", 22687658, 22693480, "-"),
    ("PEBP4", "chr8", 22713250, 22927914, "-"),
]

chr8_tracks = [
    {"path": BW_H3K27AC, "color": "#1B7837", "label": "H3K27ac",
     "height": 2.5, "ymax": 60},
    {"path": BW_PRO_PLUS, "color": "#2166AC", "label": "PRO-seq (+)",
     "height": 2.5, "ymax": 25},
    {"path": BW_ER, "color": "#D6604D", "label": "ERα",
     "height": 2.0, "ymax": 5},
    {"path": BW_FOXA1, "color": "#762A83", "label": "FOXA1",
     "height": 2.0, "ymax": 30},
]

make_browser_plot("chr8", CHR8_CENTER, CHR8_HALF,
                  ly_loops_chr8, erna_bed_chr8, chr8_genes,
                  chr8_tracks,
                  f"{OUT_PREFIX}_chr8_EGR3_LY_v2.pdf",
                  title="Supplementary — LY: eRNA chr8:22,624,675 → EGR3 (75 kb loop)",
                  condition="LY")

print("\n" + "=" * 65)
print("  All figures generated successfully!")
print("=" * 65)
