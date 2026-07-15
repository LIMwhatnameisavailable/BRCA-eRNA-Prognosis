#!/usr/bin/env python3
"""Step 1: Find HiChIP loops overlapping 10 prognostic eRNA loci (±500 bp)."""

import gzip
import csv
import os
from collections import defaultdict

# Paths
LOOP_FILE = "../GSM4763887_FitHiChIP_MCF7_DMSO.interactions_FitHiC_Q0.05.bed.gz"
ERNA_BED = "../../shared_reference/prognostic_eRNA_10.bed"
OUT_CSV = "eRNA_loops_annotated.csv"

# ── 1. Load eRNA regions (±500 bp) ──────────────────────────────────────────
ernas = []  # list of dicts
with open(ERNA_BED) as f:
    for line in f:
        line = line.strip()
        if not line:
            continue
        parts = line.split("\t")
        if len(parts) < 5:
            continue
        chrom = parts[0]
        start = int(parts[1])
        end = int(parts[2])
        name = parts[3]
        ernas.append({
            "chrom": chrom,
            "start": start,
            "end": end,
            "name": name,
            "midpoint": (start + end) // 2,
        })

print(f"Loaded {len(ernas)} eRNA regions:")
for e in ernas:
    print(f"  {e['name']}: {e['chrom']}:{e['start']:,}-{e['end']:,}")

# ── 2. Parse loop file and find overlapping loops ────────────────────────────
def overlaps(a_start, a_end, b_start, b_end):
    """Return True if intervals [a_start, a_end) and [b_start, b_end) overlap."""
    return a_start < b_end and a_end > b_start

def region_str(chrom, s, e):
    return f"{chrom}:{s:,}-{e:,}"

# Results: eRNA_name -> list of matching loop records
eRNA_hits = defaultdict(list)

with gzip.open(LOOP_FILE, "rt") as f:
    reader = csv.reader(f, delimiter="\t")
    header = next(reader)  # skip header
    # Column indices (0-based): chr1=0, s1=1, e1=2, chr2=3, s2=4, e2=5, cc=6, Q-Value_Bias=24
    for row in reader:
        if len(row) < 25:
            continue
        chr1 = row[0]
        s1 = int(row[1])
        e1 = int(row[2])
        chr2 = row[3]
        s2 = int(row[4])
        e2 = int(row[5])
        cc = int(row[6])
        qval = float(row[24])  # Q-Value_Bias

        if qval > 0.05:
            continue  # keep only significant at q<0.05 (already filtered but be safe)

        # Build the two anchor intervals
        anchor1 = (chr1, s1, e1)
        anchor2 = (chr2, s2, e2)

        for eRNA in ernas:
            hit = False
            hit_anchor = None
            other_anchor = None

            # Check if anchor1 overlaps eRNA
            if chr1 == eRNA["chrom"] and overlaps(s1, e1, eRNA["start"], eRNA["end"]):
                hit = True
                hit_anchor = anchor1
                other_anchor = anchor2
            # Check if anchor2 overlaps eRNA
            elif chr2 == eRNA["chrom"] and overlaps(s2, e2, eRNA["start"], eRNA["end"]):
                hit = True
                hit_anchor = anchor2
                other_anchor = anchor1

            if hit:
                eRNA_hits[eRNA["name"]].append({
                    "eRNA": eRNA["name"],
                    "eRNA_chrom": eRNA["chrom"],
                    "eRNA_midpoint": eRNA["midpoint"],
                    "loop_anchor1": region_str(chr1, s1, e1),
                    "loop_anchor2": region_str(chr2, s2, e2),
                    "hit_anchor": region_str(hit_anchor[0], hit_anchor[1], hit_anchor[2]),
                    "other_anchor": region_str(other_anchor[0], other_anchor[1], other_anchor[2]),
                    "other_chrom": other_anchor[0],
                    "other_start": other_anchor[1],
                    "other_end": other_anchor[2],
                    "contact_count": cc,
                    "q_value": qval,
                })

# ── 3. Write annotated CSV ──────────────────────────────────────────────────
fieldnames = [
    "eRNA", "eRNA_chrom", "eRNA_midpoint",
    "loop_anchor1", "loop_anchor2",
    "hit_anchor", "other_anchor",
    "other_chrom", "other_start", "other_end",
    "contact_count", "q_value",
]

with open(OUT_CSV, "w", newline="") as f:
    writer = csv.DictWriter(f, fieldnames=fieldnames)
    writer.writeheader()
    total_loops = 0
    for eRNA_name in sorted(eRNA_hits.keys()):
        for rec in eRNA_hits[eRNA_name]:
            writer.writerow(rec)
            total_loops += 1

print(f"\n=== RESULTS ===")
print(f"Total overlapping loops found: {total_loops}")
print(f"Output written to: {OUT_CSV}")
print()

for eRNA_name in sorted(eRNA_hits.keys()):
    recs = eRNA_hits[eRNA_name]
    # Other-end coordinate range
    other_starts = [r["other_start"] for r in recs]
    other_ends = [r["other_end"] for r in recs]
    min_other = min(other_starts)
    max_other = max(other_ends)
    other_chroms = set(r["other_chrom"] for r in recs)
    chrom_info = ", ".join(sorted(other_chroms))

    print(f"eRNA {eRNA_name}: {len(recs)} loops")
    print(f"  Other-end range: {chrom_info}:{min_other:,}-{max_other:,}")
    print(f"  Details:")
    for r in recs:
        print(f"    {r['other_anchor']}  cc={r['contact_count']}  q={r['q_value']:.4g}")

print("\nDone.")
