#!/bin/bash
# Suggestion 2.4 — AR-ESR1 Co-localization Analysis at Prognostic eRNA Loci
conda activate mytjw
set -e

WORKDIR=".."
REF="../../shared_reference"
DATA="$WORKDIR/data/GSE200300"
OUT="$WORKDIR/analysis"
mkdir -p "$OUT"

# eRNA coordinates (±2kb flanks to capture broader signal)
ERNA_BED="$REF/prognostic_eRNA_10_2kb.bed"

echo "================================================"
echo " 2.4 AR-ESR1 Colocalization Analysis"
echo "================================================"
echo "Started at: $(date)"

echo ""
echo "=== Step 2a: computeMatrix for AR, ER, Input at eRNA loci ==="
computeMatrix reference-point \
  --referencePoint center \
  -b 3000 -a 3000 \
  --binSize 50 \
  -S "$DATA/DMSO_AR.bw" "$DATA/DMSO_ER.bw" "$DATA/DMSO_Input.bw" \
  -R "$ERNA_BED" \
  --missingDataAsZero \
  -o "$OUT/AR_ER_matrix.gz" \
  --outFileSortedRegions "$OUT/AR_ER_sorted_regions.bed" \
  -p 8
echo "  ✓ computeMatrix done"

echo ""
echo "=== Step 2b: plotHeatmap — AR vs ER co-binding at eRNA loci ==="
plotHeatmap \
  -m "$OUT/AR_ER_matrix.gz" \
  -o "$OUT/Fig_2.4_AR_ER_heatmap.pdf" \
  --colorMap RdYlBu_r \
  --samplesLabel "AR (DMSO)" "ESR1 (DMSO)" "Input" \
  --regionsLabel "Prognostic eRNAs" \
  --heatmapHeight 12 \
  --heatmapWidth 6 \
  --zMin 0 --zMax 5 \
  --plotTitle "AR and ESR1 Co-binding at Prognostic eRNA Loci (MCF7)"
echo "  ✓ plotHeatmap done"

echo ""
echo "=== Step 2c: plotProfile — AR vs ER signal profile ==="
plotProfile \
  -m "$OUT/AR_ER_matrix.gz" \
  -o "$OUT/Fig_2.4_AR_ER_profile.pdf" \
  --samplesLabel "AR (DMSO)" "ESR1 (DMSO)" "Input" \
  --regionsLabel "Prognostic eRNAs" \
  --plotTitle "AR and ESR1 Binding Profile at Prognostic eRNA Loci" \
  --outFileNameData "$OUT/AR_ER_profile_data.tsv"
echo "  ✓ plotProfile done"

echo ""
echo "=== Step 2d: Per-eRNA signal quantification ==="
multiBigwigSummary BED-file \
  --BED "$ERNA_BED" \
  --bwfiles "$DATA/DMSO_AR.bw" "$DATA/DMSO_ER.bw" "$DATA/DMSO_Input.bw" \
  -o "$OUT/AR_ER_per_eRNA_summary.npz" \
  --outRawCounts "$OUT/AR_ER_per_eRNA_counts.tab" \
  -p 8
echo "  ✓ multiBigwigSummary done"

echo ""
echo "=== All deepTools steps completed at $(date) ==="