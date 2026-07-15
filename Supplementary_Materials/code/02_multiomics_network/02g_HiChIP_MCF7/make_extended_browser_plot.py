#!/usr/bin/env python3
"""Generate browser plots for newly supported eRNAs (extended window search).

Adapted from step2_v2.py but generalized for any eRNA and both DMSO/LY loop files.
"""

import gzip
import os
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Arc
import matplotlib.patches as mpatches
import matplotlib.ticker as ticker

try:
    import pyBigWig
    HAS_BIGWIG = True
except ImportError:
    HAS_BIGWIG = False
    print("WARNING: pyBigWig not available")

# ── Paths ──────────────────────────────────────────────────────────────────
LOOP_FILES = {
    "DMSO": "../GSM4763887_FitHiChIP_MCF7_DMSO.interactions_FitHiC_Q0.05.bed.gz",
    "LY":   "../GSM4763888_FitHiChIP_MCF7_LY.interactions_FitHiC_Q0.05.bed.gz",
}
ERNA_BED = "../../shared_reference/prognostic_eRNA_10.bed"
CHIP_BW = "../../2.2_ChIPseq_TF_tracks/GSE298767_Veh_H3K27ac_merge.bw"
RNA_BW = "../../3.2_RNAseq_eRNA_signal/GSE298771_Veh_Forward.bw"
OUT_PREFIX = "Fig_3.7_HiChIP"

DPI = 200

# ── 1. Load eRNA regions ──────────────────────────────────────────────────
erna_regions = {}
with open(ERNA_BED) as f:
    for line in f:
        line = line.strip()
        if not line:
            continue
        p = line.split("\t")
        chrom, start, end, name = p[0], int(p[1]), int(p[2]), p[3]
        erna_regions[name] = {"chrom": chrom, "start": start, "end": end,
                              "midpoint": (start + end) // 2}

def load_loops_two_conditions(chrom, min_coord, max_coord):
    """Load loops from both DMSO and LY within a genomic region."""
    loops_all = []
    for condition, fpath in LOOP_FILES.items():
        with gzip.open(fpath, "rt") as f:
            next(f)
            for line in f:
                p = line.strip().split("\t")
                if len(p) < 25:
                    continue
                c1, s1, e1 = p[0], int(p[1]), int(p[2])
                c2, s2, e2 = p[3], int(p[4]), int(p[5])
                cc = int(p[6])
                q = float(p[24])
                if c1 == chrom and c2 == chrom:
                    if (min_coord <= s1 <= max_coord or min_coord <= s2 <= max_coord or
                        min_coord <= e1 <= max_coord or min_coord <= e2 <= max_coord):
                        loops_all.append({"a1": s1, "e1": e1, "a2": s2, "e2": e2,
                                         "cc": cc, "q": q, "condition": condition})
    return loops_all

def read_bw(bw_path, chrom, start, end):
    if not HAS_BIGWIG or not os.path.exists(bw_path):
        return None
    bw = pyBigWig.open(bw_path)
    try:
        if chrom not in bw.chroms():
            print(f"  Warning: {chrom} not in {bw_path}")
            return None
        values = bw.values(chrom, start, end)
        positions = np.arange(start, end)
        return positions, np.nan_to_num(values, nan=0)
    finally:
        bw.close()

def make_browser_plot(chrom, center, half, loops, eRNAs, out_pdf, ename=""):
    start = center - half
    end = center + half

    print(f"\n{'='*60}")
    print(f"Plotting: {chrom}:{start:,}-{end:,}  ->  {out_pdf}")
    print(f"{'='*60}")

    chip = read_bw(CHIP_BW, chrom, start, end)
    rna = read_bw(RNA_BW, chrom, start, end)

    n_signal = sum(1 for s in [chip, rna] if s is not None)
    if n_signal == 0:
        n_signal = 2

    height_ratios = [3] * n_signal + [1.5, 0.8]
    fig_height = 2.5 + sum(height_ratios) * 0.45

    fig, axes = plt.subplots(
        len(height_ratios), 1,
        figsize=(10, fig_height),
        gridspec_kw={"height_ratios": height_ratios, "hspace": 0.05},
        sharex=True,
    )
    ax_list = axes if hasattr(axes, "__iter__") else [axes]
    ax_list = list(ax_list)

    tidx = 0

    # ── Track: H3K27ac ──
    if chip is not None:
        ax = ax_list[tidx]
        x, y = chip
        ax.fill_between(x, y, alpha=0.6, color="#1b7837", linewidth=0.3)
        ax.set_ylabel("H3K27ac", fontsize=8)
        ax.set_ylim(bottom=0)
        ax.tick_params(labelsize=7)
        tidx += 1

    # ── Track: RNA-seq ──
    if rna is not None:
        ax = ax_list[tidx]
        x, y = rna
        ax.fill_between(x, y, alpha=0.6, color="#2c7bb6", linewidth=0.3)
        ax.set_ylabel("RNA-seq\nForward", fontsize=8)
        ax.set_ylim(bottom=0)
        ax.tick_params(labelsize=7)
        tidx += 1

    # ── Track: HiChIP arcs (both conditions, colored separately) ──
    ax_arcs = ax_list[tidx]
    ax_arcs.set_xlim(start, end)
    ax_arcs.set_ylim(0, 1.0)
    ax_arcs.set_ylabel("HiChIP\nloops", fontsize=8)
    ax_arcs.tick_params(labelsize=7)
    ax_arcs.set_yticks([])

    color_map = {"DMSO": "#d7191c", "LY": "#2c7bb6"}
    n_drawn = 0

    for loop in loops:
        a1, e1, a2, e2 = loop["a1"], loop["e1"], loop["a2"], loop["e2"]
        left = min(a1, a2)
        right = max(e1, e2)
        if right <= left:
            continue
        span = right - left
        max_span = min(end - start, 500_000)
        rel_h = min(span / max_span * 0.7, 0.85)
        yc = 0.5 + rel_h * 0.3

        color = color_map.get(loop["condition"], "gray")
        arc = mpatches.Arc(
            xy=((left + right) / 2, yc),
            width=span,
            height=rel_h * 0.7,
            angle=180, theta1=0, theta2=180,
            color=color,
            linewidth=min(max(loop["cc"] / 4, 1.5), 5.0),
            alpha=0.85,
        )
        ax_arcs.add_patch(arc)
        n_drawn += 1

    if n_drawn == 0:
        ax_arcs.text((start + end) / 2, 0.5, "No HiChIP loops in view",
                     ha="center", va="center", fontsize=8, color="gray", style="italic")
    else:
        # Add legend for DMSO (red) and LY (blue)
        from matplotlib.lines import Line2D
        legend_elements = [
            Line2D([0], [0], color="#d7191c", lw=2, label="DMSO loops"),
            Line2D([0], [0], color="#2c7bb6", lw=2, label="LY loops"),
        ]
        ax_arcs.legend(handles=legend_elements, loc="upper right", fontsize=6, framealpha=0.7)

    ax_arcs.spines["left"].set_visible(False)
    ax_arcs.spines["right"].set_visible(False)
    ax_arcs.spines["top"].set_visible(False)
    tidx += 1

    # ── Track: eRNA annotation ──
    ax_bed = ax_list[tidx]
    ax_bed.set_xlim(start, end)
    ax_bed.set_ylim(-0.5, 1.5)
    ax_bed.set_yticks([])
    ax_bed.set_ylabel("eRNA", fontsize=8)
    ax_bed.spines["left"].set_visible(False)
    ax_bed.spines["right"].set_visible(False)
    ax_bed.spines["top"].set_visible(False)
    ax_bed.tick_params(labelsize=7)

    for en, er in eRNAs.items():
        if er["chrom"] != chrom:
            continue
        es, ee = er["start"], er["end"]
        if es < end and ee > start:
            vis_s = max(es, start)
            vis_e = min(ee, end)
            mid = (vis_s + vis_e) / 2
            ax_bed.barh(0.5, vis_e - vis_s, left=vis_s, height=0.7,
                        color="#f4a582", edgecolor="#d6604d", linewidth=0.5)
            lbl = en.replace(f"eRNA_{chrom}_", "eRNA_")
            ax_bed.text(mid, 0.5, lbl, ha="center", va="center",
                        fontsize=5, fontweight="bold", color="#333333")

    ax_bed.set_xlabel(f"{chrom} position (bp)", fontsize=8)
    ax_bed.xaxis.set_major_formatter(ticker.FuncFormatter(lambda x, p: f"{x/1e6:.2f}"))
    tidx += 1

    plt.tight_layout()
    plt.savefig(out_pdf, dpi=DPI, bbox_inches="tight")
    plt.close()
    print(f"  ✅ Saved: {out_pdf}")


if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python3 make_extended_browser_plot.py <eRNA_name>")
        print("Example: python3 make_extended_browser_plot.py eRNA_chr1_155158995")
        sys.exit(1)

    ename = sys.argv[1]
    if ename not in erna_regions:
        print(f"Error: {ename} not found in eRNA list")
        sys.exit(1)

    e = erna_regions[ename]
    chrom = e["chrom"]
    center = e["midpoint"]
    print(f"Generating browser plot for {ename}")
    print(f"  Chromosome: {chrom}")
    print(f"  Center: {center:,}")

    # Load loops (±150kb)
    half = 150_000
    loops = load_loops_two_conditions(chrom, center - half, center + half)
    print(f"  Loops found: {len(loops)}")
    for l in loops:
        print(f"    {l['condition']}: {chrom}:{l['a1']:,}-{l['e1']:,} <-> {chrom}:{l['a2']:,}-{l['e2']:,}  cc={l['cc']}")

    label = f"{chrom}_{center//1000}kb"
    pdf_path = f"{OUT_PREFIX}_{label}_extended.pdf"

    # Only this eRNA for the annotation track
    target_erna = {ename: e}

    make_browser_plot(chrom, center, half, loops, target_erna, pdf_path, ename)

    # Also make a wider view
    pdf_path_wide = f"{OUT_PREFIX}_{label}_extended_wide.pdf"
    make_browser_plot(chrom, center, 200_000, loops, target_erna, pdf_path_wide, ename)

    print(f"\n✅ Done. Output: {pdf_path}, {pdf_path_wide}")