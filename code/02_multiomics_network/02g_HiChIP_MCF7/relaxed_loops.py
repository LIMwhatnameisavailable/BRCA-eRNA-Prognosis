#!/usr/bin/env python3
"""Step 3.7: Relaxed HiChIP loop finder — multi-condition intersection.

Conditions tested per loop file (DMSO / LY):
  - Q<0.05, eRNA ±500bp (original)
  - Q<0.05, eRNA ±5kb
  - Q<0.1,  eRNA ±500bp
  - Q<0.1,  eRNA ±5kb

Output:
  - summary table printed to stdout
  - eRNA_loops_annotated_relaxed.csv (best condition)
  - pyGenomeTracks plots for eRNAs with ≥2 loops
"""

import gzip
import csv
import os
import subprocess
from collections import defaultdict

# ── Paths ──────────────────────────────────────────────────────────────────────
LOOP_FILES = {
    "DMSO": "../GSM4763887_FitHiChIP_MCF7_DMSO.interactions_FitHiC_Q0.05.bed.gz",
    "LY":   "../GSM4763888_FitHiChIP_MCF7_LY.interactions_FitHiC_Q0.05.bed.gz",
}
ERNA_BED = "../../shared_reference/prognostic_eRNA_10.bed"
CHIP_BW = "../../2.2_ChIPseq_TF_tracks/GSE298767_Veh_H3K27ac_merge.bw"
RNA_BW = "../../3.2_RNAseq_eRNA_signal/GSE298771_Veh_Forward.bw"

# ── 1. Load eRNA regions (base ±500bp from original BED) ───────────────────────
ernas_base = []
with open(ERNA_BED) as f:
    for line in f:
        line = line.strip()
        if not line:
            continue
        parts = line.split("\t")
        chrom, start, end, name = parts[0], int(parts[1]), int(parts[2]), parts[3]
        erna_mid = (start + end) // 2
        ernas_base.append({
            "chrom": chrom, "start": start, "end": end,
            "name": name, "midpoint": erna_mid,
        })

# Build expanded windows
def make_windows(halfwidth):
    windows = []
    for e in ernas_base:
        windows.append({
            "chrom": e["chrom"],
            "start": e["midpoint"] - halfwidth,
            "end":   e["midpoint"] + halfwidth,
            "name":  e["name"],
            "midpoint": e["midpoint"],
        })
    return windows

ernas_500bp = make_windows(500)
ernas_5kb   = make_windows(5000)

def overlaps(a_start, a_end, b_start, b_end):
    return a_start < b_end and a_end > b_start

# ── 2. Process each condition ───────────────────────────────────────────────────
# Results: dict keyed by (sample, q_threshold, window_label) -> eRNA_name -> list of loop records
results = defaultdict(lambda: defaultdict(list))

for sample, fpath in LOOP_FILES.items():
    print(f"\n{'='*70}")
    print(f"Processing {sample}: {os.path.basename(fpath)}")
    print(f"{'='*70}")

    with gzip.open(fpath, "rt") as f:
        reader = csv.reader(f, delimiter="\t")
        header = next(reader)
        lines = list(reader)

    total_lines = len(lines)
    print(f"  Total loops in file: {total_lines}")

    for q_thresh in [0.05, 0.1]:
        for win_label, ernas_win in [("500bp", ernas_500bp), ("5kb", ernas_5kb)]:
            key = (sample, q_thresh, win_label)
            n_matched = 0
            for row in lines:
                if len(row) < 25:
                    continue
                qval = float(row[24])
                if qval > q_thresh:
                    continue
                chr1, s1, e1 = row[0], int(row[1]), int(row[2])
                chr2, s2, e2 = row[3], int(row[4]), int(row[5])
                cc = int(row[6])

                for erna in ernas_win:
                    hit = False
                    if chr1 == erna["chrom"] and overlaps(s1, e1, erna["start"], erna["end"]):
                        hit = True
                        hit_anchor = f"{chr1}:{s1:,}-{e1:,}"
                        other_anchor = f"{chr2}:{s2:,}-{e2:,}"
                        other_chrom, other_start, other_end = chr2, s2, e2
                    elif chr2 == erna["chrom"] and overlaps(s2, e2, erna["start"], erna["end"]):
                        hit = True
                        hit_anchor = f"{chr2}:{s2:,}-{e2:,}"
                        other_anchor = f"{chr1}:{s1:,}-{e1:,}"
                        other_chrom, other_start, other_end = chr1, s1, e1

                    if hit:
                        results[key][erna["name"]].append({
                            "sample": sample,
                            "eRNA": erna["name"],
                            "eRNA_chrom": erna["chrom"],
                            "eRNA_midpoint": erna["midpoint"],
                            "loop_anchor1": f"{chr1}:{s1:,}-{e1:,}",
                            "loop_anchor2": f"{chr2}:{s2:,}-{e2:,}",
                            "hit_anchor": hit_anchor,
                            "other_anchor": other_anchor,
                            "other_chrom": other_chrom,
                            "other_start": other_start,
                            "other_end": other_end,
                            "contact_count": cc,
                            "q_value": qval,
                        })
                        n_matched += 1
            print(f"  Q<{q_thresh}  eRNA±{win_label}: {n_matched} overlapping loops "
                  f"({len(results[key])} eRNAs hit)")

# ── 3. Print summary table ────────────────────────────────────────────────────
print(f"\n{'='*110}")
print("SUMMARY: Loop counts per eRNA per condition")
print(f"{'='*110}")
header_row = f"{'eRNA':<28}"
for sample in ["DMSO", "LY"]:
    for q in [0.05, 0.1]:
        for win in ["500bp", "5kb"]:
            header_row += f"  {sample}_{q}_{win:<10}"
print(header_row)
print("-" * 110)

eRNA_names = [e["name"] for e in ernas_base]
total_row = f"{'TOTAL':<28}"
grand_totals = {}
for sample in ["DMSO", "LY"]:
    for q in [0.05, 0.1]:
        for win in ["500bp", "5kb"]:
            key = (sample, q, win)
            grand_totals[key] = sum(len(v) for v in results[key].values())
            total_row += f"  {grand_totals[key]:<20}"

for ename in eRNA_names:
    row = f"{ename:<28}"
    for sample in ["DMSO", "LY"]:
        for q in [0.05, 0.1]:
            for win in ["500bp", "5kb"]:
                key = (sample, q, win)
                cnt = len(results[key].get(ename, []))
                row += f"  {cnt:<20}"
    print(row)

print("-" * 110)
print(total_row)

# ── 4. Write relaxed annotated CSV (DMSO Q<0.05 ±5kb or whichever) ────────────
# Pick the most inclusive condition
best_key = ("DMSO", 0.1, "5kb")
# Or use DMSO Q<0.05 ±5kb as default
default_key = ("DMSO", 0.05, "5kb")
# Let's use the condition with most loops
best_sample, best_q, best_win = max(grand_totals, key=grand_totals.get)
best_key = (best_sample, best_q, best_win)
print(f"\nBest condition for annotation CSV: {best_sample} Q<{best_q} ±{best_win} "
      f"({grand_totals[best_key]} loops)")

relaxed_csv = "eRNA_loops_annotated_relaxed.csv"
fieldnames = [
    "sample", "eRNA", "eRNA_chrom", "eRNA_midpoint",
    "loop_anchor1", "loop_anchor2", "hit_anchor", "other_anchor",
    "other_chrom", "other_start", "other_end",
    "contact_count", "q_value",
]
all_records = []
for key, edict in sorted(results.items()):
    for ename, recs in edict.items():
        all_records.extend(recs)

with open(relaxed_csv, "w", newline="") as f:
    writer = csv.DictWriter(f, fieldnames=fieldnames)
    writer.writeheader()
    for rec in all_records:
        writer.writerow(rec)
print(f"  → Wrote {relaxed_csv} ({len(all_records)} total records across all conditions)")

# Also write per-condition CSVs
for sample in ["DMSO", "LY"]:
    for q in [0.05, 0.1]:
        for win in ["500bp", "5kb"]:
            key = (sample, q, win)
            fname = f"eRNA_loops_{sample}_Q{q}_W{win}.csv"
            recs = []
            for ename, rects in results[key].items():
                recs.extend(rects)
            if recs:
                with open(fname, "w", newline="") as f:
                    w = csv.DictWriter(f, fieldnames=fieldnames)
                    w.writeheader()
                    for r in recs:
                        w.writerow(r)
                print(f"  → Wrote {fname} ({len(recs)} records)")

# ── 5. Find eRNAs with ≥2 loops in any condition ──────────────────────────────
candidates = set()
for key, edict in results.items():
    for ename, recs in edict.items():
        if len(recs) >= 2:
            candidates.add(ename)
            print(f"\nCandidate for viz: {ename} — {len(recs)} loops in {key}")

# ── 6. pyGenomeTracks for candidates ────────────────────────────────────────
def create_bedpe(sample, q_thresh, win_label, chrom, min_c, max_c, out_bedpe):
    """Create BEDPE from a specific condition for visualization."""
    key = (sample, q_thresh, win_label)
    loops_bedpe = []
    seen_pairs = set()
    with gzip.open(LOOP_FILES[sample], "rt") as f:
        next(f)
        for line in f:
            p = line.strip().split("\t")
            if len(p) < 25:
                continue
            c1, s1, e1 = p[0], int(p[1]), int(p[2])
            c2, s2, e2 = p[3], int(p[4]), int(p[5])
            qv = float(p[24])
            if qv > q_thresh:
                continue
            if c1 == chrom and c2 == chrom:
                if (min_c <= s1 <= max_c or min_c <= s2 <= max_c):
                    pair = (min(s1, s2), max(e1, e2))
                    if pair not in seen_pairs:
                        seen_pairs.add(pair)
                        loops_bedpe.append((c1, s1, e1, c2, s2, e2, int(p[6])))
    with open(out_bedpe, "w") as f:
        for i, (c1, s1, e1, c2, s2, e2, cc) in enumerate(loops_bedpe):
            f.write(f"{c1}\t{s1}\t{e1}\t{c2}\t{s2}\t{e2}\tloop_{i}\t{cc}\n")
    return loops_bedpe

for ename in sorted(candidates):
    # Find which eRNA this is
    e_info = [e for e in ernas_base if e["name"] == ename][0]
    chrom = e_info["chrom"]
    mid = e_info["midpoint"]
    label = f"{chrom}_{mid//1000}kb"

    print(f"\n{'='*60}")
    print(f"Generating browser plot for {ename} ({label})")
    print(f"{'='*60}")

    # Create BEDPE for loops in this region (±100kb)
    half = 100_000
    region_start = mid - half
    region_end = mid + half

    # Use DMSO Q<0.05 ±5kb loops for consistency
    viz_bedpe = f"viz_{sample}_{label}.bedpe"
    viz_loops = create_bedpe("DMSO", 0.05, "5kb", chrom, region_start, region_end, viz_bedpe)
    print(f"  Loops in region: {len(viz_loops)}")

    # Create eRNA BED for this eRNA
    erna_bed = f"viz_erna_{label}.bed"
    with open(erna_bed, "w") as f:
        f.write(f"{chrom}\t{e_info['start']}\t{e_info['end']}\t{ename}\t0\t.\n")

    # Create ini file
    ini_file = f"tracks_{label}.ini"
    pdf_file = f"Fig_3.7_HiChIP_{label}.pdf"

    ini_content = f"""
[spacer]

[hic_loops]
file = {viz_bedpe}
file_type = arcs
title = HiChIP DMSO loops
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

[erna]
file = {erna_bed}
title = eRNA region
color = darkorange
height = 1.5
file_type = bed
show_labels = yes
show_data_range = no
line_width = 3
"""
    with open(ini_file, "w") as f:
        f.write(ini_content.strip())

    # Try pyGenomeTracks; fall back to matplotlib-based script
    region_str = f"{chrom}:{region_start:,}-{region_end:,}"
    r = subprocess.run(["which", "pyGenomeTracks"], capture_output=True, text=True)
    if r.returncode == 0:
        cmd = ["pyGenomeTracks", "--tracks", ini_file, "--region", region_str,
               "--output", pdf_file, "--dpi", 150]
        result = subprocess.run(cmd, capture_output=True, text=True)
        if result.returncode == 0:
            print(f"  ✅ pyGenomeTracks plot: {pdf_file}")
        else:
            print(f"  ❌ pyGenomeTracks error: {result.stderr[:300]}")
    else:
        print(f"  ⚠ pyGenomeTracks not available, skipping automated plot")
        print(f"  To generate plot manually:\n    pyGenomeTracks --tracks {ini_file} --region {region_str} --output {pdf_file} --dpi 150")

print(f"\n{'='*70}")
print("All done. Check outputs in the current directory.")
print(f"{'='*70}")
