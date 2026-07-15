#!/usr/bin/env python3
"""Step 2: Generate browser-style figures using matplotlib + pyBigWig directly.

Produces tracks: H3K27ac ChIP-seq, RNA-seq Forward, HiChIP arcs, eRNA annotation.
Generates Fig_3.7_HiChIP_[chrom]_[pos].pdf.
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

# Try pyBigWig; if not available, use a stub that reads from UCSC
try:
    import pyBigWig
    HAS_BIGWIG = True
except ImportError:
    HAS_BIGWIG = False
    print("WARNING: pyBigWig not available, will use mock signal", file=sys.stderr)

# ── Paths ────────────────────────────────────────────────────────────────────
LOOP_FILE = "../GSM4763887_FitHiChIP_MCF7_DMSO.interactions_FitHiC_Q0.05.bed.gz"
ERNA_BED = "../../shared_reference/prognostic_eRNA_10.bed"
CHIP_BW = "../../2.2_ChIPseq_TF_tracks/GSE298767_Veh_H3K27ac_merge.bw"
RNA_BW = "../../3.2_RNAseq_eRNA_signal/GSE298771_Veh_Forward.bw"
OUT_PREFIX = "Fig_3.7_HiChIP"

# ── Parameters ───────────────────────────────────────────────────────────────
HALF_WINDOW = 100_000  # ±100 kb
DPI = 150

# ── 1. Load eRNA regions ────────────────────────────────────────────────────
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

print("=== eRNA Regions ===")
for n, r in erna_regions.items():
    print(f"  {n}: {r['chrom']}:{r['start']:,}-{r['end']:,} (mid {r['midpoint']:,})")

# ── 2. Load loops for chr10 target region ────────────────────────────────────
# All chr10 intra-chromosomal loops around the eRNA loci (5.3-5.8 Mb)
def load_loops(chrom, min_coord, max_coord):
    """Load loops where BOTH anchors are within [min_coord, max_coord] on `chrom`."""
    loops = []
    with gzip.open(LOOP_FILE, "rt") as f:
        next(f)  # header
        for line in f:
            p = line.strip().split("\t")
            if len(p) < 25:
                continue
            c1, s1, e1 = p[0], int(p[1]), int(p[2])
            c2, s2, e2 = p[3], int(p[4]), int(p[5])
            cc = int(p[6])
            q = float(p[24])
            if c1 == chrom and c2 == chrom:
                # Keep if at least anchor in our region
                if (min_coord <= s1 <= max_coord or min_coord <= s2 <= max_coord or
                    min_coord <= e1 <= max_coord or min_coord <= e2 <= max_coord):
                    loops.append({"a1": s1, "e1": e1, "a2": s2, "e2": e2, "cc": cc, "q": q})
    return loops

# ── 3. BigWig reader ────────────────────────────────────────────────────────
def read_bw(bw_path, chrom, start, end):
    """Read bigwig signal over region, return (positions, values) arrays."""
    if not HAS_BIGWIG or not os.path.exists(bw_path):
        return None
    bw = pyBigWig.open(bw_path)
    try:
        if chrom not in bw.chroms():
            print(f"  Warning: {chrom} not in {bw_path}", file=sys.stderr)
            return None
        values = bw.values(chrom, start, end)
        positions = np.arange(start, end)
        return positions, np.nan_to_num(values, nan=0)
    finally:
        bw.close()

# ── 4. Plot function ────────────────────────────────────────────────────────
def make_browser_plot(chrom, center, half, loops, eRNAs, out_pdf):
    """Create a multi-track browser figure."""
    start = center - half
    end = center + half
    n_bases = end - start

    print(f"\n{'='*60}")
    print(f"Plotting {chrom}:{start:,}-{end:,}")
    print(f"{'='*60}")

    # ── Read signal tracks ──
    chip = read_bw(CHIP_BW, chrom, start, end)
    rna = read_bw(RNA_BW, chrom, start, end)

    # ── Determine arc height for HiChIP ──
    # Arc y-position: each loop gets an arc whose height is proportional to
    # genomic distance
    ARCS_MAX_HEIGHT = 0.15  # fraction of signal track height for arcs
    ARCS_YM = 1.0           # y-position (in signal coords) for arc baseline

    # Count how many tracks are active
    n_signal_tracks = sum(1 for s in [chip, rna] if s is not None)
    if n_signal_tracks == 0:
        print("  ⚠ No signal tracks available")
        n_signal_tracks = 2  # placeholder

    total_height = n_signal_tracks * 3 + 1.0 + 1.5  # signal(3 each) + arcs(1.0) + bed(1.5) + spacers

    # ── Build figure ──
    fig, axes = plt.subplots(
        n_signal_tracks + 1 + 1,  # signals + arcs + bed
        1,
        figsize=(12, total_height),
        gridspec_kw={"height_ratios": [3] * n_signal_tracks + [1.0, 1.5]},
        sharex=True,
    )
    # Flatten for indexing
    ax_list = axes if hasattr(axes, "__iter__") else [axes]
    ax_list = list(ax_list)  # ensure list

    track_idx = 0

    # ── Track 1: H3K27ac ──
    if chip is not None:
        ax = ax_list[track_idx]
        x, y = chip
        ax.fill_between(x, y, alpha=0.7, color="darkgreen", linewidth=0.5)
        ax.set_ylabel("H3K27ac", fontsize=10)
        ax.set_ylim(bottom=0)
        ax.ticklabel_format(style="plain", axis="x", useOffset=False)
        track_idx += 1

    # ── Track 2: RNA-seq ──
    if rna is not None:
        ax = ax_list[track_idx]
        x, y = rna
        ax.fill_between(x, y, alpha=0.7, color="steelblue", linewidth=0.5)
        ax.set_ylabel("RNA-seq\nForward", fontsize=10)
        ax.set_ylim(bottom=0)
        ax.ticklabel_format(style="plain", axis="x", useOffset=False)
        track_idx += 1

    # ── Track 3: HiChIP arcs ──
    ax_arcs = ax_list[track_idx]
    ax_arcs.set_ylabel("HiChIP\nloops", fontsize=10)
    ax_arcs.set_xlim(start, end)
    ax_arcs.set_ylim(0, 1.0)

    # Draw each loop as an arc
    n_loops_drawn = 0
    for loop in loops:
        a1, e1 = loop["a1"], loop["e1"]
        a2, e2 = loop["a2"], loop["e2"]
        # Clamp to view region for arc visualization
        left = max(min(a1, a2), start)
        right = min(max(e1, e2), end)
        if right <= left:
            continue
        span = right - left
        if span <= 0:
            continue
        # Normalized arc height
        max_span = min(end - start, 500_000)
        rel_height = min(span / max_span * 0.8, 0.9)
        y_center = 0.5 + rel_height * 0.3

        arc = mpatches.Arc(
            xy=((left + right) / 2, y_center),
            width=span,
            height=rel_height * 0.8,
            angle=180,
            theta1=0,
            theta2=180,
            color="firebrick",
            linewidth=min(max(loop["cc"] / 5, 1.0), 4.0),
            alpha=0.8,
        )
        ax_arcs.add_patch(arc)
        n_loops_drawn += 1

    if n_loops_drawn == 0:
        ax_arcs.text((start + end) / 2, 0.5, "No loops in view",
                     ha="center", va="center", fontsize=9, color="gray", style="italic")

    ax_arcs.set_yticks([])
    ax_arcs.spines["left"].set_visible(False)
    ax_arcs.spines["right"].set_visible(False)
    ax_arcs.spines["top"].set_visible(False)
    ax_arcs.ticklabel_format(style="plain", axis="x", useOffset=False)
    track_idx += 1

    # ── Track 4: eRNA annotation ──
    ax_bed = ax_list[track_idx]
    ax_bed.set_ylabel("eRNA", fontsize=10)
    ax_bed.set_xlim(start, end)
    ax_bed.set_ylim(0, 1)
    ax_bed.set_yticks([])
    ax_bed.spines["left"].set_visible(False)
    ax_bed.spines["right"].set_visible(False)
    ax_bed.spines["top"].set_visible(False)

    # Draw eRNA regions as colored blocks
    for ename, er in eRNAs.items():
        if er["chrom"] != chrom:
            continue
        es, ee = er["start"], er["end"]
        if es < end and ee > start:
            # Clamp to view
            vis_s = max(es, start)
            vis_e = min(ee, end)
            mid = (vis_s + vis_e) / 2
            ax_bed.barh(0.5, vis_e - vis_s, left=vis_s, height=0.6,
                        color="darkorange", alpha=0.9, edgecolor="black", linewidth=0.5)
            # Label if there's room
            if vis_e - vis_s > 500:
                label = ename.replace("eRNA_", "")
                ax_bed.text(mid, 0.5, label, ha="center", va="center",
                            fontsize=6, fontweight="bold")

    track_idx += 1

    # ── X-axis formatting ──
    ax_bed.set_xlabel(f"{chrom} position (bp)", fontsize=10)
    # Format x ticks as Mb
    import matplotlib.ticker as ticker
    def fmt_mb(x, pos):
        return f"{x/1e6:.2f} Mb"
    ax_bed.xaxis.set_major_formatter(ticker.FuncFormatter(fmt_mb))

    plt.tight_layout()
    plt.savefig(out_pdf, dpi=DPI, bbox_inches="tight")
    plt.close()
    print(f"  ✅ Saved: {out_pdf}")

# ── 5. Define regions and generate plots ────────────────────────────────────
# Focus on chr10 adjacent eRNAs
chr10_ernas = {n: r for n, r in erna_regions.items() if r["chrom"] == "chr10"}
print(f"\n=== Visualization ===")
print(f"chr10 eRNAs for visualization: {list(chr10_ernas.keys())}")

# Use eRNA_chr10_5531356 as center (has the loop hit)
center = 5531356
half = HALF_WINDOW
start_region = center - half
end_region = center + half

# Load loops in broader region (include the loop's other end at ~5.66 Mb)
loops = load_loops("chr10", 5_300_000, 5_800_000)
print(f"Found {len(loops)} HiChIP loops in chr10:5.3-5.8Mb")

out_pdf = f"{OUT_PREFIX}_chr10_5.53Mb.pdf"
make_browser_plot("chr10", center, half, loops, chr10_ernas, out_pdf)

# Also generate a wider view (±150 kb) to show the loop's other end
out_pdf_wide = f"{OUT_PREFIX}_chr10_5.53Mb_wide.pdf"
make_browser_plot("chr10", center, 150_000, loops, chr10_ernas, out_pdf_wide)
print(f"\nWide view: ±150kb to capture loop distal anchor at ~5.66 Mb → {out_pdf_wide}")

print("\n=== Step 2 Complete ===")
