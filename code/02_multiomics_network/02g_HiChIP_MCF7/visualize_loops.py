#!/usr/bin/env python3
"""Step 3.7 v2: Generate browser plots for relaxed HiChIP analysis.

Uses LY treatment loops (more hits) and shows both chr10 eRNAs in one view.
Falls back to matplotlib if pyGenomeTracks unavailable.
"""

import gzip
import csv
import os
import subprocess
from collections import defaultdict

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.ticker as ticker

try:
    import pyBigWig
    HAS_BIGWIG = True
except ImportError:
    HAS_BIGWIG = False

# ── Paths ──────────────────────────────────────────────────────────────────
LOOP_FILE_LY = "../GSM4763888_FitHiChIP_MCF7_LY.interactions_FitHiC_Q0.05.bed.gz"
LOOP_FILE_DMSO = "../GSM4763887_FitHiChIP_MCF7_DMSO.interactions_FitHiC_Q0.05.bed.gz"
ERNA_BED = "../../shared_reference/prognostic_eRNA_10.bed"
CHIP_BW = "../../2.2_ChIPseq_TF_tracks/GSE298767_Veh_H3K27ac_merge.bw"
RNA_BW = "../../3.2_RNAseq_eRNA_signal/GSE298771_Veh_Forward.bw"
OUT_PREFIX = "Fig_3.7_HiChIP"
DPI = 200

# ── 1. Load eRNA regions ──────────────────────────────────────────────────
ernas = {}
with open(ERNA_BED) as f:
    for line in f:
        line = line.strip()
        if not line:
            continue
        p = line.split("\t")
        chrom, start, end, name = p[0], int(p[1]), int(p[2]), p[3]
        ernas[name] = {"chrom": chrom, "start": start, "end": end,
                       "midpoint": (start + end) // 2}

# ── 2. Load loops ─────────────────────────────────────────────────────────
def load_loops(fpath, chrom, min_coord, max_coord, q_thresh=0.05, halfwidth=5000):
    """Load loops from a given file. Apply Q<q_thresh filter and expand eRNA window by halfwidth."""
    loops = []
    # First load eRNA windows
    eRNA_windows = []
    for ename, e in ernas.items():
        eRNA_windows.append({
            "chrom": e["chrom"], "start": e["midpoint"] - halfwidth,
            "end": e["midpoint"] + halfwidth, "name": ename,
        })

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
            if q > q_thresh:
                continue
            if c1 == chrom and c2 == chrom:
                in_region = False
                for r in eRNA_windows:
                    if r["chrom"] == chrom:
                        a_start, a_end = r["start"], r["end"]
                        if (a_start <= s1 <= a_end or a_start <= e1 <= a_end or
                            a_start <= s2 <= a_end or a_start <= e2 <= a_end):
                            in_region = True
                            break
                if in_region:
                    loops.append({"a1": s1, "e1": e1, "a2": s2, "e2": e2, "cc": cc, "q": q})
    return loops

# ── 3. BigWig reader ──────────────────────────────────────────────────────
def read_bw(bw_path, chrom, start, end):
    if not HAS_BIGWIG or not os.path.exists(bw_path):
        return None
    bw = pyBigWig.open(bw_path)
    try:
        if chrom not in bw.chroms():
            return None
        values = bw.values(chrom, start, end)
        positions = np.arange(start, end)
        return positions, np.nan_to_num(values, nan=0)
    finally:
        bw.close()

# ── 4. Plot function ──────────────────────────────────────────────────────
def make_browser_plot(chrom, center, half, loops, eRNAs_subset, out_pdf, title=""):
    start = center - half
    end = center + half

    print(f"\n  Plotting {chrom}:{start:,}-{end:,}  →  {out_pdf}")

    chip = read_bw(CHIP_BW, chrom, start, end)
    rna = read_bw(RNA_BW, chrom, start, end)

    n_sig = sum(1 for s in [chip, rna] if s is not None)
    if n_sig == 0:
        n_sig = 2

    height_ratios = [3] * n_sig + [1.0, 0.8]
    fig_height = 2.5 + sum(height_ratios) * 0.45

    fig, axes = plt.subplots(len(height_ratios), 1,
                             figsize=(10, fig_height),
                             gridspec_kw={"height_ratios": height_ratios, "hspace": 0.05},
                             sharex=True)
    ax_list = axes if hasattr(axes, "__iter__") else [axes]
    ax_list = list(ax_list)
    tidx = 0

    # H3K27ac
    if chip is not None:
        ax = ax_list[tidx]
        x, y = chip
        ax.fill_between(x, y, alpha=0.6, color="#1b7837", linewidth=0.3)
        ax.set_ylabel("H3K27ac", fontsize=8)
        ax.set_ylim(bottom=0)
        ax.tick_params(labelsize=7)
        tidx += 1

    # RNA-seq
    if rna is not None:
        ax = ax_list[tidx]
        x, y = rna
        ax.fill_between(x, y, alpha=0.6, color="#2c7bb6", linewidth=0.3)
        ax.set_ylabel("RNA-seq\nForward", fontsize=8)
        ax.set_ylim(bottom=0)
        ax.tick_params(labelsize=7)
        tidx += 1

    # HiChIP arcs
    ax_arcs = ax_list[tidx]
    ax_arcs.set_xlim(start, end)
    ax_arcs.set_ylim(0, 1.0)
    ax_arcs.set_ylabel("HiChIP\nloops(LY)", fontsize=8)
    ax_arcs.tick_params(labelsize=7)
    ax_arcs.set_yticks([])

    n_drawn = 0
    for loop in loops:
        left = min(loop["a1"], loop["a2"])
        right = max(loop["e1"], loop["e2"])
        if right <= left:
            continue
        span = right - left
        max_span = min(end - start, 500_000)
        rel_h = min(span / max_span * 0.7, 0.85)
        yc = 0.5 + rel_h * 0.3
        arc = mpatches.Arc(xy=((left+right)/2, yc), width=span, height=rel_h*0.7,
                           angle=180, theta1=0, theta2=180,
                           color="#d7191c", linewidth=min(max(loop["cc"]/4, 1.5), 5.0),
                           alpha=0.85)
        ax_arcs.add_patch(arc)
        n_drawn += 1
    if n_drawn == 0:
        ax_arcs.text((start+end)/2, 0.5, "No loops in view",
                     ha="center", va="center", fontsize=8, color="gray", style="italic")
    ax_arcs.spines["left"].set_visible(False)
    ax_arcs.spines["right"].set_visible(False)
    ax_arcs.spines["top"].set_visible(False)
    tidx += 1

    # eRNA annotation
    ax_bed = ax_list[tidx]
    ax_bed.set_xlim(start, end)
    ax_bed.set_ylim(-0.5, 1.5)
    ax_bed.set_yticks([])
    ax_bed.set_ylabel("eRNA", fontsize=8)
    ax_bed.spines["left"].set_visible(False)
    ax_bed.spines["right"].set_visible(False)
    ax_bed.spines["top"].set_visible(False)
    ax_bed.tick_params(labelsize=7)

    for ename, er in eRNAs_subset.items():
        if er["chrom"] != chrom:
            continue
        es, ee = er["start"], er["end"]
        if es < end and ee > start:
            vis_s = max(es, start)
            vis_e = min(ee, end)
            ax_bed.barh(0.5, vis_e - vis_s, left=vis_s, height=0.7,
                        color="#f4a582", edgecolor="#d6604d", linewidth=0.5)
            lbl = ename.replace("eRNA_chr10_", "eRNA_")
            ax_bed.text((vis_s+vis_e)/2, 0.5, lbl, ha="center", va="center",
                        fontsize=6, fontweight="bold", color="#333333")

    ax_bed.set_xlabel(f"{chrom} position (bp)", fontsize=8)
    ax_bed.xaxis.set_major_formatter(ticker.FuncFormatter(lambda x, p: f"{x/1e6:.2f}"))
    tidx += 1

    plt.tight_layout()
    plt.savefig(out_pdf, dpi=DPI, bbox_inches="tight")
    plt.close()
    print(f"  ✅ Saved: {out_pdf}")

# ── 5. Generate viz ───────────────────────────────────────────────────────
# chr10 eRNAs
chr10_ernas = {n: r for n, r in ernas.items() if r["chrom"] == "chr10"}
print("=== Generating chr10 combined browser plots (LY loops) ===")

# Load LY loops in chr10 5.3-5.8 Mb region using relaxed criteria
ly_loops = load_loops(LOOP_FILE_LY, "chr10", 5_300_000, 5_800_000,
                       q_thresh=0.05, halfwidth=5000)
print(f"LY loops in chr10 region: {len(ly_loops)}")

# Also load DMSO loops for comparison
dmso_loops = load_loops(LOOP_FILE_DMSO, "chr10", 5_300_000, 5_800_000,
                         q_thresh=0.05, halfwidth=5000)
print(f"DMSO loops in chr10 region: {len(dmso_loops)}")

# Plot with LY loops, ±100kb around the eRNA cluster midpoint
center = 5531356
half = 100_000
chr10_ernas_subset = {n: r for n, r in ernas.items() if r["chrom"] == "chr10"}

# Plot 1: LY loops
make_browser_plot("chr10", center, half, ly_loops, chr10_ernas_subset,
                  f"{OUT_PREFIX}_chr10_5.53Mb_LY_relaxed.pdf",
                  title="LY Q<0.05 eRNA±5kb")

# Plot 2: DMSO loops (for comparison)
make_browser_plot("chr10", center, half, dmso_loops, chr10_ernas_subset,
                  f"{OUT_PREFIX}_chr10_5.53Mb_DMSO_relaxed.pdf",
                  title="DMSO Q<0.05 eRNA±5kb")

# Plot 3: Wider view (±150kb) with LY loops
make_browser_plot("chr10", center, 150_000, ly_loops, chr10_ernas_subset,
                  f"{OUT_PREFIX}_chr10_5.53Mb_LY_wide.pdf")

# ── 6. eRNA_chr8_22624675 has 1 loop in LY — also check if worth plotting ──
e8_loops = load_loops(LOOP_FILE_LY, "chr8", 22_600_000, 22_700_000,
                       q_thresh=0.05, halfwidth=5000)
print(f"\nchr8 eRNA region LY loops: {len(e8_loops)}")
if e8_loops:
    e8_erna = {n: r for n, r in ernas.items() if n == "eRNA_chr8_22624675"}
    make_browser_plot("chr8", 22624675, 100_000, e8_loops, e8_erna,
                      f"{OUT_PREFIX}_chr8_22.62Mb_LY.pdf")

print("\n=== Step 3.7 Complete ===")
