#!/bin/bash
# =========================================================
# 3.2_RunOnSeq - T47D PRO-seq with CORRECT hg38 coordinates
# =========================================================
set -e
START_TIME=$(date '+%Y-%m-%d %H:%M:%S')

BASE=".."
DATA="${BASE}/data/GSE227243_T47D"
ANALYSIS="${BASE}/analysis"
REF="../../shared_reference"
# conda environment needs to be activated manually before running

# eval "$CONDA_ACTIVATE"

ERNA_BED="${REF}/prognostic_eRNA_10.bed"

# T47D PRO-seq bigWig files
REP1_PLUS="${DATA}/GSM7093765_T47D_DMSO_Rep1_processed_plus.bw"
REP1_MINUS="${DATA}/GSM7093765_T47D_DMSO_Rep1_processed_minus.bw"
REP2_PLUS="${DATA}/GSM7093766_T47D_DMSO_Rep2_processed_plus.bw"
REP2_MINUS="${DATA}/GSM7093766_T47D_DMSO_Rep2_processed_minus.bw"

mkdir -p "$ANALYSIS"
cd "$ANALYSIS"

echo "=============================================================="
echo " 3.2_RunOnSeq T47D PRO-seq - Recompute with hg38"
echo " Started: $START_TIME"
echo "=============================================================="

echo ""
echo "=== Step 1: computeMatrix (reference-point, center ±2kb) ==="
computeMatrix reference-point \
  --referencePoint center \
  -R "$ERNA_BED" \
  -S ${REP1_PLUS} ${REP1_MINUS} ${REP2_PLUS} ${REP2_MINUS} \
  -b 2000 -a 2000 \
  --binSize 10 \
  --missingDataAsZero \
  -o "${ANALYSIS}/T47D_PROseq_matrix.gz" \
  -p 8
echo "  ✓ computeMatrix done"

echo ""
echo "=== Step 2: plotHeatmap final ==="
plotHeatmap \
  -m "${ANALYSIS}/T47D_PROseq_matrix.gz" \
  -o "${ANALYSIS}/T47D_PROseq_heatmap_final.pdf" \
  --colorMap YlOrRd \
  --zMin 0 \
  --samplesLabel "Rep1 +" "Rep1 -" "Rep2 +" "Rep2 -" \
  --regionsLabel "Prognostic eRNAs" \
  --heatmapHeight 8 \
  --heatmapWidth 6 \
  --plotTitle "T47D PRO-seq (hg38)"
echo "  ✓ heatmap done"

echo ""
echo "=== Step 3: plotProfile final ==="
plotProfile \
  -m "${ANALYSIS}/T47D_PROseq_matrix.gz" \
  -o "${ANALYSIS}/T47D_PROseq_profile_final.pdf" \
  --samplesLabel "Rep1 +" "Rep1 -" "Rep2 +" "Rep2 -" \
  --perGroup \
  --plotType se \
  --colors blue red cyan pink \
  --plotTitle ""
echo "  ✓ profile done"

echo ""
echo "=== Step 4: Plus strand only (aggregated) ==="
computeMatrix reference-point \
  --referencePoint center \
  -R "$ERNA_BED" \
  -S ${REP1_PLUS} ${REP2_PLUS} \
  -b 2000 -a 2000 \
  --binSize 10 \
  --missingDataAsZero \
  -o "${ANALYSIS}/T47D_PROseq_matrix_plus_final.gz" \
  -p 8

plotHeatmap \
  -m "${ANALYSIS}/T47D_PROseq_matrix_plus_final.gz" \
  -o "${ANALYSIS}/T47D_PROseq_heatmap_plus_v5.pdf" \
  --colorMap YlOrRd \
  --zMin 0 \
  --samplesLabel "Rep1 +" "Rep2 +" \
  --heatmapHeight 8 --heatmapWidth 4 \
  --plotTitle ""

echo ""
echo "=== Step 5: Minus strand only ==="
computeMatrix reference-point \
  --referencePoint center \
  -R "$ERNA_BED" \
  -S ${REP1_MINUS} ${REP2_MINUS} \
  -b 2000 -a 2000 \
  --binSize 10 \
  --missingDataAsZero \
  -o "${ANALYSIS}/T47D_PROseq_matrix_minus_final.gz" \
  -p 8

# Summary
END_TIME=$(date '+%Y-%m-%d %H:%M:%S')
echo ""
echo "=============================================================="
echo " 3.2_RunOnSeq RECOMPUTE COMPLETE"
echo " Start: $START_TIME"
echo " End:   $END_TIME"
echo "=============================================================="
echo ""
echo "Generated files:"
ls -lh "${ANALYSIS}/T47D_PROseq_heatmap_final.pdf" "${ANALYSIS}/T47D_PROseq_profile_final.pdf" 2>/dev/null
