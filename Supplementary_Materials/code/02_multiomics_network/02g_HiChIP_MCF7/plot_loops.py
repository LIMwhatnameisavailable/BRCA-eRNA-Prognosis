#!/usr/bin/env python3
"""Step 2: Generate pyGenomeTracks browser-style figures for eRNA loci with HiChIP loops."""

import gzip
import os
import subprocess
import csv
from collections import defaultdict

# Paths
LOOP_FILE = "../GSM4763887_FitHiChIP_MCF7_DMSO.interactions_FitHiC_Q0.05.bed.gz"
ERNA_BED = "../../shared_reference/prognostic_eRNA_10.bed"
CHIP_BW = "../../2.2_ChIPseq_TF_tracks/GSE298767_Veh_H3K27ac_merge.bw"
RNA_BW = "../../3.2_RNAseq_eRNA_signal/GSE298771_Veh_Forward.bw"

# ── 1. Read eRNA loop data ──────────────────────────────────────────────────
eRNA_loop_counts = defaultdict(int)
eRNA_info = {}  # name -> {chrom, midpoint}

with open("eRNA_loops_annotated.csv") as f:
    reader = csv.DictReader(f)
    for row in reader:
        name = row["eRNA"]
        eRNA_loop_counts[name] += 1
        if name not in eRNA_info:
            eRNA_info[name] = {
                "chrom": row["eRNA_chrom"],
                "midpoint": int(row["eRNA_midpoint"]),
            }

print("=== eRNA Loop Counts (from annotated CSV) ===")
for name, cnt in sorted(eRNA_loop_counts.items()):
    info = eRNA_info[name]
    print(f"  {name}: {cnt} loops at {info['chrom']}:{info['midpoint']:,}")

# ── 2. Select eRNAs for visualization ───────────────────────────────────────
has_2plus = [n for n, c in eRNA_loop_counts.items() if c >= 2]
print(f"\n=== Selection ===")
print(f"eRNAs with ≥2 loops: {has_2plus if has_2plus else 'NONE — relaxing to ≥1 for visualization'}")

has_1plus = [n for n, c in eRNA_loop_counts.items() if c >= 1]
print(f"eRNAs with ≥1 loop: {has_1plus}")

# ── 3. Prepare data files for pyGenomeTracks ──────────────────────────────
# Use BED file for both chr10 adjacent eRNAs
erna_track_bed = "eRNA_chr10_track.bed"
with open(ERNA_BED, "r") as fin, open(erna_track_bed, "w") as fout:
    for line in fin:
        line = line.strip()
        if not line:
            continue
        parts = line.split("\t")
        chrom = parts[0]
        if chrom == "chr10":
            fout.write(line + "\n")
print(f"\nWrote {erna_track_bed}: chr10 eRNA regions")

# Create BEDPE for chr10 intra-chromosomal loops in 5.3-5.8 Mb region
wide_bedpe = "chr10_viz_loops.bedpe"
with gzip.open(LOOP_FILE, "rt") as fin, open(wide_bedpe, "w") as fout:
    next(fin)  # skip header
    n = 0
    for line in fin:
        parts = line.strip().split("\t")
        if len(parts) < 7:
            continue
        chr1, s1, e1 = parts[0], int(parts[1]), int(parts[2])
        chr2, s2, e2 = parts[3], int(parts[4]), int(parts[5])
        cc = parts[6]
        # chr10 intra-chromosomal, spanning roughly 5.3-5.8 Mb
        if chr1 == "chr10" and chr2 == "chr10" and s1 >= 5_300_000 and e2 <= 5_800_000:
            # BEDPE: chr1 s1 e1 chr2 s2 e2 [name] [score]
            fout.write(f"{chr1}\t{s1}\t{e1}\t{chr2}\t{s2}\t{e2}\tloop_{n}\t{cc}\n")
            n += 1
print(f"{wide_bedpe}: {n} loop records")

# ── 4. Define visualization regions ────────────────────────────────────────
# The chr10 adjacent eRNAs (eRNA_chr10_5531356 and eRNA_chr10_5528926)
# are only ~2.4 kb apart. Visualize as one region ±100 kb from summit.
viz_regions = [
    {
        "label": "chr10_5.53Mb",
        "chrom": "chr10",
        "center": 5531356,
        "half_window": 100_000,
        "desc": "chr10 eRNA region (±100 kb)",
    }
]

# ── 5. Generate plots ──────────────────────────────────────────────────────
for region in viz_regions:
    chrom = region["chrom"]
    center = region["center"]
    half = region["half_window"]
    start = center - half
    end = center + half
    label = region["label"]
    region_str = f"{chrom}:{start:,}-{end:,}"

    ini_file = f"tracks_{label}.ini"
    pdf_file = f"Fig_3.7_HiChIP_{label}.pdf"

    ini_content = f"""
[spacer]

[hic_loops]
file = {wide_bedpe}
file_type = arcs
title = HiChIP loops (Q<0.05)
color = firebrick
line_width = 2.5
type = arcs
overlay_previous = no
show_data_range = no

[spacer]

[h3k27ac]
file = {CHIP_BW}
title = H3K27ac ChIP-seq
height = 3
color = darkgreen
min_value = 0
show_data_range = yes
file_type = bigwig

[spacer]

[rna_fwd]
file = {RNA_BW}
title = RNA-seq Forward
height = 3
color = steelblue
min_value = 0
show_data_range = yes
file_type = bigwig

[spacer]

[erna_regions]
file = {erna_track_bed}
title = eRNA (±500bp)
color = darkorange
height = 1.5
file_type = bed
show_labels = yes
show_data_range = no
line_width = 3
"""

    with open(ini_file, "w") as f:
        f.write(ini_content.strip())

    print(f"\n{'='*60}")
    print(f"Region: {region_str}")
    print(f"Config: {ini_file}")
    print(f"Output: {pdf_file}")
    print(f"{'='*60}")

    # Check pyGenomeTracks available
    r = subprocess.run(["which", "pyGenomeTracks"], capture_output=True, text=True)
    if r.returncode == 0:
        cmd = [
            "pyGenomeTracks",
            "--tracks", ini_file,
            "--region", region_str,
            "--output", pdf_file,
            "--dpi", 150,
        ]
        print(f"Running: {' '.join(cmd)}")
        result = subprocess.run(cmd, capture_output=True, text=True)
        if result.returncode == 0:
            print(f"  ✅ Plot generated: {pdf_file}")
        else:
            print(f"  ❌ Error (stderr):")
            print(f"     {result.stderr[:1000]}")
            print(f"  stdout: {result.stdout[:500]}")
    else:
        print("  ❌ pyGenomeTracks not found in PATH")

print("\n=== Step 2 Complete ===")
