#!/usr/bin/env python3
"""
Phase 4: Motif parsing + subtype classification for AR-ESR1 co-occurrence analysis.
Reads ARE/ERE motif hits, H3K27ac counts, and existing ChIP-seq enrichment table.
Produces a unified summary CSV.
"""
import os, re, csv
import numpy as np

WORKDIR = ".."
ANALYSIS = os.path.join(WORKDIR, "analysis")

# ── 1. Parse motif hits ───────────────────────────────────────────────────
def parse_motif_hits(hit_file):
    """Parse HOMER scanMotifGenomeWide.pl BED output."""
    hits = {}
    with open(hit_file) as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) < 6:
                continue
            # Parse locus name from BED name field
            name_full = parts[0]  # e.g. eRNA_chr1_155158995::chr1:155184519-155188519
            locus = name_full.split('::')[0]
            if locus not in hits:
                hits[locus] = []
            score = float(parts[4])
            strand = parts[5]
            start = int(parts[1])
            end = int(parts[2])
            hits[locus].append({'start': start, 'end': end, 'score': score, 'strand': strand})
    return hits

are_hits = parse_motif_hits(os.path.join(ANALYSIS, "ARE_hits.bed"))
ere_hits = parse_motif_hits(os.path.join(ANALYSIS, "ERE_hits.bed"))
foxa1_hits = parse_motif_hits(os.path.join(ANALYSIS, "FOXA1_hits.bed"))
gata3_hits = parse_motif_hits(os.path.join(ANALYSIS, "GATA3_hits.bed"))

# All 10 loci
LOCI = [
    "eRNA_chr1_155158995", "eRNA_chr10_5528926", "eRNA_chr10_5531356",
    "eRNA_chr12_13371038", "eRNA_chr3_11236700", "eRNA_chr3_138070534",
    "eRNA_chr8_22624675", "eRNA_chr9_71398719", "eRNA_chr9_71398939",
    "eRNA_chr9_114689796"
]

# ── 2. Build motif summary per locus ──────────────────────────────────────
motif_rows = []
for locus in LOCI:
    are = sorted(are_hits.get(locus, []), key=lambda x: -x['score'])
    ere = sorted(ere_hits.get(locus, []), key=lambda x: -x['score'])
    foxa1 = sorted(foxa1_hits.get(locus, []), key=lambda x: -x['score'])
    gata3 = sorted(gata3_hits.get(locus, []), key=lambda x: -x['score'])
    row = {
        'eRNA': locus,
        'ARE_hits': len(are),
        'ARE_best_score': round(are[0]['score'], 3) if are else 0,
        'ERE_hits': len(ere),
        'ERE_best_score': round(ere[0]['score'], 3) if ere else 0,
        'FOXA1_hits': len(foxa1),
        'FOXA1_best_score': round(foxa1[0]['score'], 3) if foxa1 else 0,
        'GATA3_hits': len(gata3),
        'GATA3_best_score': round(gata3[0]['score'], 3) if gata3 else 0,
        'motif_cooccurrence': len(are) > 0 and len(ere) > 0,
    }
    # Print detailed hit info
    if are:
        print(f"  {locus}: ARE hits at " + ", ".join(f"{h['start']}-{h['end']}({h['strand']},score={h['score']:.1f})" for h in are))
    if ere:
        print(f"  {locus}: ERE hits at " + ", ".join(f"{h['start']}-{h['end']}({h['strand']},score={h['score']:.1f})" for h in ere))
    if foxa1:
        print(f"  {locus}: FOXA1 hits at " + ", ".join(f"{h['start']}-{h['end']}({h['strand']},score={h['score']:.1f})" for h in foxa1))
    if gata3:
        print(f"  {locus}: GATA3 hits at " + ", ".join(f"{h['start']}-{h['end']}({h['strand']},score={h['score']:.1f})" for h in gata3))
    motif_rows.append(row)

# ── Build mapping from H3K27ac naming to eRNA locus names ────────────────
BED_FILE = os.path.join(WORKDIR, "../shared_reference/prognostic_eRNA_10_2kb.bed")
locus_map = {}
with open(BED_FILE) as f:
    for line in f:
        parts = line.strip().split('\t')
        chrom = parts[0]
        start = int(parts[1])
        end = int(parts[2])
        name = parts[3]
        mid = (start + end) // 2
        h3_name = f"eRNA_{chrom}_{mid}"
        locus_map[h3_name] = name

# ── 3. Read existing ChIP-seq enrichment table ─────────────────────────────
# Uses midpoint-based naming (e.g., eRNA_chr10_5486963) - map to eRNA locus names
chip_data = {}
chip_path = os.path.join(ANALYSIS, "AR_ER_enrichment_table.csv")
with open(chip_path) as f:
    reader = csv.DictReader(f)
    for row in reader:
        chip_locus = row['eRNA']  # Uses H3K27ac-style naming (midpoint)
        # Map to canonical eRNA locus name using locus_map (h3_name → eRNA_name)
        mapped = locus_map.get(chip_locus, chip_locus)
        chip_data[mapped] = row

# ── 4. Read H3K27ac counts for subtype annotation ─────────────────────────
SUB_NAMES = ['MCF7_H3K27ac', 'SKBR3_H3K27ac', 'MB453_H3K27ac', 'MB231_H3K27ac', 'Hs578T_H3K27ac']

h3k27ac = {}
with open(os.path.join(ANALYSIS, "subtype_H3K27ac_counts.tab")) as f:
    header = f.readline().strip().split('\t')
    col_idx = {c.strip("'"): i for i, c in enumerate(header)}
    for line in f:
        parts = line.strip().split('\t')
        chrom = parts[0].strip("'")
        start = int(parts[1])
        end = int(parts[2])
        mid = (start + end) // 2
        h3_name = f"eRNA_{chrom}_{mid}"
        locus = locus_map.get(h3_name, h3_name)
        vals = {}
        for sn in SUB_NAMES:
            vals[sn] = float(parts[col_idx[sn]]) if col_idx[sn] < len(parts) else 0
        vals['TNBC_mean'] = np.mean([vals.get('MB231_H3K27ac',0), vals.get('Hs578T_H3K27ac',0)])
        h3k27ac[locus] = vals

# ── 5. Subtype classification ─────────────────────────────────────────────
# Uses z-score normalized H3K27ac across cell lines per locus
def classify_subtype(vals):
    """Classify a locus by its dominant H3K27ac activity pattern."""
    arr = np.array([vals[s] for s in SUB_NAMES])
    if np.max(arr) < 0.05:  # very low signal everywhere
        return "Inactive"

    # Find which cell line has the highest signal
    max_idx = np.argmax(arr)
    max_val = arr[max_idx]
    second_max = sorted(arr)[-2] if len(arr) > 1 else 0

    # Check if it's truly enriched in one subtype vs others
    ratio = max_val / second_max if second_max > 0 else float('inf')

    cell_map = {
        0: ("Luminal-enriched", "MCF7"),
        1: ("HER2-enriched", "SKBR3"),
        2: ("LAR-enriched", "MB453"),
        3: ("TNBC-enriched", "MB231"),
        4: ("TNBC-enriched", "Hs578T"),
    }

    if ratio > 1.5:
        label, _ = cell_map[max_idx]
        return label
    else:
        # Two or more subtypes have similar signal
        # Check if MCF7+LAR (AR-relevant) or broader
        if arr[0] > 0.05 and arr[2] > 0.05:  # MCF7 and MB453
            return "Pan-subtype (Luminal+LAR)"
        elif np.sum(arr > 0.05) >= 3:
            return "Pan-subtype"
        else:
            label, _ = cell_map[max_idx]
            return label

# ── 6. Merge everything into unified table ─────────────────────────────────
print("\n" + "="*90)
print("  UNIFIED SUMMARY: Motif co-occurrence + ChIP-seq enrichment + Subtype")
print("="*90)

header = (f"{'Locus':<30} {'ARE':>4} {'ERE':>4} {'Co-occur':>10} "
          f"{'AR/Input':>9} {'ER/Input':>9} "
          f"{'MCF7':>8} {'SKBR3':>8} {'MB453':>8} {'MB231':>8} {'Hs578T':>8} "
          f"{'Subtype':<25}")
print(header)
print("-"*len(header))

unified_rows = []
for mr in motif_rows:
    locus = mr['eRNA']

    # ChIP-seq data
    cd = chip_data.get(locus, {})
    ar_enr = float(cd.get('AR_enrichment', 0))
    er_enr = float(cd.get('ER_enrichment', 0))

    # H3K27ac data
    h3 = h3k27ac.get(locus, {})
    mcf7_v = h3.get('MCF7_H3K27ac', 0)
    skbr3_v = h3.get('SKBR3_H3K27ac', 0)
    mb453_v = h3.get('MB453_H3K27ac', 0)
    mb231_v = h3.get('MB231_H3K27ac', 0)
    hs578t_v = h3.get('Hs578T_H3K27ac', 0)

    subtype = classify_subtype(h3)

    cooccur = "✓ BOTH" if mr['motif_cooccurrence'] else ("ARE only" if mr['ARE_hits'] > 0 else ("ERE only" if mr['ERE_hits'] > 0 else "Neither"))

    foxa1_present = "✓" if mr['FOXA1_hits'] > 0 else ""
    gata3_present = "✓" if mr['GATA3_hits'] > 0 else ""

    print(f"{locus:<30} {mr['ARE_hits']:>4} {mr['ERE_hits']:>4} {cooccur:>10} "
          f"{ar_enr:>9.3f} {er_enr:>9.3f} "
          f"{mcf7_v:>8.3f} {skbr3_v:>8.3f} {mb453_v:>8.3f} {mb231_v:>8.3f} {hs578t_v:>8.3f} "
          f"{subtype:<25}")
    print(f"  {'':>30} FOXA1:{mr['FOXA1_hits']} {foxa1_present:>1}  GATA3:{mr['GATA3_hits']} {gata3_present:>1}")

    unified_rows.append({
        'eRNA': locus,
        'ARE_hits': mr['ARE_hits'],
        'ARE_best_score': mr['ARE_best_score'],
        'ERE_hits': mr['ERE_hits'],
        'ERE_best_score': mr['ERE_best_score'],
        'FOXA1_hits': mr['FOXA1_hits'],
        'FOXA1_best_score': mr['FOXA1_best_score'],
        'FOXA1_present': 1 if mr['FOXA1_hits'] > 0 else 0,
        'GATA3_hits': mr['GATA3_hits'],
        'GATA3_best_score': mr['GATA3_best_score'],
        'GATA3_present': 1 if mr['GATA3_hits'] > 0 else 0,
        'motif_cooccurrence': cooccur,
        'AR_enrichment': round(ar_enr, 3),
        'ER_enrichment': round(er_enr, 3),
        'ChIP_seq_co_bound': cd.get('co_bound', 'False'),
        'MCF7_H3K27ac': round(mcf7_v, 4),
        'SKBR3_H3K27ac': round(skbr3_v, 4),
        'MB453_H3K27ac': round(mb453_v, 4),
        'MB231_H3K27ac': round(mb231_v, 4),
        'Hs578T_H3K27ac': round(hs578t_v, 4),
        'subtype_classification': subtype,
    })

# ── 7. Write unified CSV ──────────────────────────────────────────────────
out_csv = os.path.join(ANALYSIS, "unified_summary.csv")
with open(out_csv, 'w', newline='') as f:
    w = csv.DictWriter(f, fieldnames=unified_rows[0].keys())
    w.writeheader()
    w.writerows(unified_rows)
print(f"\n  ✓ Written: {out_csv}")

# ── 8. Summary stats ──────────────────────────────────────────────────────
n_cooccur = sum(1 for r in unified_rows if r['motif_cooccurrence'] == '✓ BOTH')
n_are_only = sum(1 for r in unified_rows if r['motif_cooccurrence'] == 'ARE only')
n_ere_only = sum(1 for r in unified_rows if r['motif_cooccurrence'] == 'ERE only')
n_neither = sum(1 for r in unified_rows if r['motif_cooccurrence'] == 'Neither')
n_er_enriched = sum(1 for r in unified_rows if r['ER_enrichment'] > 1.5)
n_ere_positive_er_chip = sum(1 for r in unified_rows if r['ERE_hits'] > 0 and r['ER_enrichment'] > 1.5)

print("\n" + "="*60)
print("  KEY FINDINGS")
print("="*60)
print(f"  Motif co-occurrence (ARE+ERE): {n_cooccur}/10 loci")
print(f"  ARE only: {n_are_only}/10")
print(f"  ERE only: {n_ere_only}/10")
print(f"  Neither:  {n_neither}/10")
print(f"  ER ChIP-seq enriched (>1.5x): {n_er_enriched}/10")
print(f"  ERE motif + ER ChIP-seq positive: {n_ere_positive_er_chip}/10")
print()

# Highlight chr10 loci
print("  chr10 loci detail:")
for r in unified_rows:
    if 'chr10' in r['eRNA']:
        print(f"    {r['eRNA']}: motif={r['motif_cooccurrence']}, "
              f"AR/Input={r['AR_enrichment']}, ER/Input={r['ER_enrichment']}, "
              f"subtype={r['subtype_classification']}")

# Highlight subtype distribution
print("\n  Subtype distribution:")
subtypes = {}
for r in unified_rows:
    s = r['subtype_classification']
    subtypes[s] = subtypes.get(s, 0) + 1
for s, c in sorted(subtypes.items(), key=lambda x: -x[1]):
    print(f"    {s}: {c}/10")

print()