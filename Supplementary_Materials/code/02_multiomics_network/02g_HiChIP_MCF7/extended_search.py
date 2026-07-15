#!/usr/bin/env python3
"""Step 4: Extended window search for eRNA loop coverage.

Tests multiple windows against all loop files to find which eRNAs have
loop support under any condition.

Conditions:
  - Samples: DMSO, LY
  - Q thresholds: 0.05, 0.1
  - Windows: ±500 bp, ±2 kb, ±5 kb, ±10 kb

Output:
  - eRNA_loop_coverage_extended.csv  — per-eRNA, per-window loop counts
  - eRNA_no_loop_support.txt         — eRNAs with zero support in ANY condition
  - Console summary
"""

import gzip
import csv
import os
import sys
from collections import defaultdict

# ── Paths ──────────────────────────────────────────────────────────────────────
LOOP_FILES = {
    "DMSO": "../GSM4763887_FitHiChIP_MCF7_DMSO.interactions_FitHiC_Q0.05.bed.gz",
    "LY":   "../GSM4763888_FitHiChIP_MCF7_LY.interactions_FitHiC_Q0.05.bed.gz",
}
ERNA_BED = "../../shared_reference/prognostic_eRNA_10.bed"
OUT_CSV = "eRNA_loop_coverage_extended.csv"
OUT_NOTXT = "eRNA_no_loop_support.txt"

# ── 1. Load eRNA regions ──────────────────────────────────────────────────────
ernas = []
with open(ERNA_BED) as f:
    for line in f:
        line = line.strip()
        if not line:
            continue
        parts = line.split("\t")
        if len(parts) < 5:
            continue
        chrom, start, end, name = parts[0], int(parts[1]), int(parts[2]), parts[3]
        midpoint = (start + end) // 2
        ernas.append({
            "chrom": chrom, "start": start, "end": end,
            "name": name, "midpoint": midpoint,
        })

print(f"Loaded {len(ernas)} eRNA regions:")
for e in ernas:
    print(f"  {e['name']}: {e['chrom']}:{e['start']:,}-{e['end']:,}  (mid={e['midpoint']:,})")

# ── 2. Define windows to test ──────────────────────────────────────────────────
WINDOWS = [500, 2000, 5000, 10000]  # half-widths in bp
Q_THRESHOLDS = [0.05, 0.1]
SAMPLES = ["DMSO", "LY"]

def make_windows(ernas_list, halfwidth):
    """Create expanded window regions."""
    windows = []
    for e in ernas_list:
        windows.append({
            "chrom": e["chrom"],
            "start": max(0, e["midpoint"] - halfwidth),
            "end":   e["midpoint"] + halfwidth,
            "name":  e["name"],
            "midpoint": e["midpoint"],
        })
    return windows

def overlaps(a_start, a_end, b_start, b_end):
    return a_start < b_end and a_end > b_start

# Pre-compute window regions for each halfwidth
window_regions = {hw: make_windows(ernas, hw) for hw in WINDOWS}

# ── 3. Scan loops ─────────────────────────────────────────────────────────────
# Result key: (sample, q_thresh, halfwidth, ename) -> list of loops
hits = defaultdict(list)

for sample in SAMPLES:
    fpath = LOOP_FILES[sample]
    print(f"\n{'='*70}")
    print(f"Loading {sample}: {os.path.basename(fpath)}")
    print(f"{'='*70}")

    with gzip.open(fpath, "rt") as f:
        reader = csv.reader(f, delimiter="\t")
        header = next(reader)
        lines = list(reader)

    print(f"  Total loops in file: {len(lines)}")

    for q_thresh in Q_THRESHOLDS:
        for hw in WINDOWS:
            ernas_win = window_regions[hw]
            key_prefix = (sample, q_thresh, hw)
            n_matched = 0
            eRNA_matched = set()

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
                        o_chrom, o_start, o_end = chr2, s2, e2
                    elif chr2 == erna["chrom"] and overlaps(s2, e2, erna["start"], erna["end"]):
                        hit = True
                        hit_anchor = f"{chr2}:{s2:,}-{e2:,}"
                        other_anchor = f"{chr1}:{s1:,}-{e1:,}"
                        o_chrom, o_start, o_end = chr1, s1, e1

                    if hit:
                        key = key_prefix + (erna["name"],)
                        hits[key].append({
                            "sample": sample,
                            "eRNA": erna["name"],
                            "eRNA_chrom": erna["chrom"],
                            "eRNA_midpoint": erna["midpoint"],
                            "window_bp": hw,
                            "q_threshold": q_thresh,
                            "hit_anchor": hit_anchor,
                            "other_anchor": other_anchor,
                            "other_chrom": o_chrom,
                            "other_start": o_start,
                            "other_end": o_end,
                            "contact_count": cc,
                            "q_value": qval,
                        })
                        n_matched += 1
                        eRNA_matched.add(erna["name"])

            win_label = f"{hw//1000}kb" if hw >= 1000 else f"{hw}bp"
            print(f"  Q<{q_thresh}  ±{win_label}: {n_matched:>4} loops, "
                  f"{len(eRNA_matched):>2} eRNAs hit")

# ── 4. Build per-eRNA summary matrix ──────────────────────────────────────────
# For each eRNA across all conditions: count loops
# Also track "best condition" for each window

print(f"\n{'='*110}")
print("PER-ERNA LOOP COVERAGE MATRIX")
print(f"{'='*110}")

# Header
col_labels = []
for s in SAMPLES:
    for q in Q_THRESHOLDS:
        for hw in WINDOWS:
            wl = f"{hw//1000}kb" if hw >= 1000 else f"{hw}bp"
            col_labels.append(f"{s}_Q{q}_{wl}")

header_row = f"{'eRNA':<28} " + "  ".join(f"{c:<16}" for c in col_labels) + "  best_condition"
print(header_row)
print("-" * 110)

summary_rows = []

for e in ernas:
    ename = e["name"]
    row_counts = []
    n_conditions_with_loops = 0
    best_condition = "none"
    max_loops = 0

    for s in SAMPLES:
        for q in Q_THRESHOLDS:
            for hw in WINDOWS:
                key = (s, q, hw, ename)
                cnt = len(hits.get(key, []))
                row_counts.append(cnt)
                if cnt > 0:
                    n_conditions_with_loops += 1
                    if cnt > max_loops:
                        max_loops = cnt
                        wl = f"{hw//1000}kb" if hw >= 1000 else f"{hw}bp"
                        best_condition = f"{s}_Q{q}_{wl}"

    # Format best_condition
    if best_condition == "none":
        best_str = "NO SUPPORT"
    elif max_loops >= 5:
        best_str = f"{best_condition}({max_loops})"
    else:
        best_str = f"{best_condition}({max_loops})"

    row_str = f"{ename:<28} " + "  ".join(f"{c:<16}" for c in row_counts) + f"  {best_str}"
    print(row_str)
    summary_rows.append({
        "eRNA": ename,
        "chrom": e["chrom"],
        "midpoint": e["midpoint"],
        **{f"{s}_Q{q}_{wl}": row_counts[i]
           for i, (s, q, hw) in enumerate([(s,q,hw) for s in SAMPLES for q in Q_THRESHOLDS for hw in WINDOWS])
           for wl in [f"{hw//1000}kb" if hw >= 1000 else f"{hw}bp"]},
        "any_loop_support": "YES" if max_loops > 0 else "NO",
        "best_condition": best_condition if max_loops > 0 else "NONE",
        "max_loops": max_loops,
        "n_conditions_with_loops": n_conditions_with_loops,
    })

# Compute totals per column
total_row_counts = [0] * len(col_labels)
for e in ernas:
    i = 0
    for s in SAMPLES:
        for q in Q_THRESHOLDS:
            for hw in WINDOWS:
                key = (s, q, hw, e["name"])
                total_row_counts[i] += len(hits.get(key, []))
                i += 1

total_str = f"{'TOTAL':<28} " + "  ".join(f"{c:<16}" for c in total_row_counts)
print("-" * 110)
print(total_str)

# ── 5. Previously known vs newly supported ────────────────────────────────────
# Original strict: eRNA_chr10_5531356, eRNA_chr10_5528926, eRNA_chr8_22624675
known_supported = {"eRNA_chr10_5531356", "eRNA_chr10_5528926", "eRNA_chr8_22624675"}

newly_supported = []
still_no_support = []

for r in summary_rows:
    if r["any_loop_support"] == "YES":
        if r["eRNA"] not in known_supported:
            newly_supported.append(r["eRNA"])
            print(f"\n  ✨ NEW SUPPORT: {r['eRNA']} — {r['best_condition']}")
    else:
        still_no_support.append(r["eRNA"])
        print(f"\n  ❌ NO SUPPORT: {r['eRNA']}")

# ── 6. Write summary CSV ──────────────────────────────────────────────────────
col_fieldnames = [f"{s}_Q{q}_{wl}" for s in SAMPLES for q in Q_THRESHOLDS
                  for hw in WINDOWS for wl in [f"{hw//1000}kb" if hw>=1000 else f"{hw}bp"]]
fieldnames = ["eRNA", "chrom", "midpoint"] + col_fieldnames + [
    "any_loop_support", "best_condition", "max_loops", "n_conditions_with_loops"
]

with open(OUT_CSV, "w", newline="") as f:
    writer = csv.DictWriter(f, fieldnames=fieldnames)
    writer.writeheader()
    for r in summary_rows:
        writer.writerow(r)

print(f"\n✅ Wrote: {OUT_CSV} ({len(summary_rows)} rows)")

# ── 7. Write "no support" list ────────────────────────────────────────────────
with open(OUT_NOTXT, "w") as f:
    f.write(f"# eRNAs with NO loop support in ANY condition\n")
    f.write(f"# Samples: DMSO, LY | Q thresholds: 0.05, 0.1 | Windows: ±500bp, ±2kb, ±5kb, ±10kb\n")
    f.write(f"# Total eRNAs: {len(ernas)}\n")
    f.write(f"# Supported (≥1 loop): {len(ernas) - len(still_no_support)}\n")
    f.write(f"# Unsupported: {len(still_no_support)}\n")
    f.write("#\n")
    for ename in still_no_support:
        e = [e for e in ernas if e["name"] == ename][0]
        f.write(f"{e['chrom']}\t{e['start']}\t{e['end']}\t{ename}\n")

print(f"✅ Wrote: {OUT_NOTXT} ({len(still_no_support)} unsupported eRNAs)")

# ── 8. Final summary ──────────────────────────────────────────────────────────
print(f"\n{'='*70}")
print("COVERAGE SUMMARY")
print(f"{'='*70}")
print(f"  Total eRNAs:               {len(ernas)}")
print(f"  Previously supported:      {len(known_supported)} (chr8_22624675, chr10_5528926, chr10_5531356)")
print(f"  Newly supported (extended): {len(newly_supported)}")
if newly_supported:
    for ename in newly_supported:
        print(f"    → {ename}: {[r for r in summary_rows if r['eRNA']==ename][0]['best_condition']}")
print(f"  Still no support (all conditions): {len(still_no_support)}")
for ename in still_no_support:
    print(f"    → {ename}")
print(f"\nCoverage: {len(ernas) - len(still_no_support)}/{len(ernas)} "
      f"({100 * (len(ernas) - len(still_no_support)) / len(ernas):.0f}%)")
print("\nDone.")