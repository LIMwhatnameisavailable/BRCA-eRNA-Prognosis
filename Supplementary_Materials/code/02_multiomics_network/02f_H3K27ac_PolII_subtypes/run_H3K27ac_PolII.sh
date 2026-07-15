#!/bin/bash
# =========================================================
# 3.3_H3K27ac_PolII_subtypes - Re-run with CORRECT hg38 coords
# =========================================================
set -e
START_TIME=$(date '+%Y-%m-%d %H:%M:%S')

BASE=".."
REF="../../shared_reference"
ANALYSIS="${BASE}/analysis"
# conda environment needs to be activated manually before running
# eval "$CONDA_ACTIVATE"

ERNA_BED="${REF}/prognostic_eRNA_10_2kb.bed"
ERNA_BED_1kb="${REF}/prognostic_eRNA_10.bed"
mkdir -p "$ANALYSIS"
cd "$ANALYSIS"

echo "=============================================================="
echo " 3.3 Multi-subtype H3K27ac - Recompute with hg38"
echo " Started: $START_TIME"
echo "=============================================================="

# MCF7 (GSE245868) - 5 tracks
MCF7_K27="${BASE}/GSM7850113_Sorted_E_H3K27acNormalised.bw"
MCF7_K4me1="${BASE}/GSM7850111_Sorted_ctrl_H3K4meNormalised.bw"
MCF7_K4me3="${BASE}/GSM7850112_Sorted_ctrl_H3K4me3Normalised.bw"
MCF7_POLII="${BASE}/GSM7850116_Sorted_ctrl_PolIINormalised.bw"
MCF7_P300="${BASE}/GSM7850133_Sorted_p300Normalised.bw"

# Subtype H3K27ac - use first replicate (UT1 for SKBR3/MB453, WO1 for others)
SKBR3_K27="${BASE}/GSM7989319_SKBR3_K27_UT1.bw"
MB453_K27="${BASE}/MB453_H3K27ac_mean.bw"
MB231_K27="${BASE}/GSM7989345_MB231_K27_WO1.bw"
HS578T_K27="${BASE}/Hs578T_H3K27ac_mean.bw"

# =========================================================
# v1: 8-sample scale-regions
# =========================================================
echo ""
echo "=============================================================="
echo " v1: 8-sample scale-regions"
echo "=============================================================="

computeMatrix scale-regions \
  -R "$ERNA_BED" \
  -S ${MCF7_K27} ${MCF7_K4me1} ${MCF7_K4me3} ${MCF7_POLII} ${MCF7_P300} ${SKBR3_K27} ${MB231_K27} ${HS578T_K27} \
  --regionBodyLength 1000 \
  -b 2000 -a 2000 \
  --binSize 100 \
  --missingDataAsZero \
  -o "${ANALYSIS}/multisubtype_matrix.gz" \
  -p 8
echo "  ✓ v1 computeMatrix done"

plotHeatmap \
  -m "${ANALYSIS}/multisubtype_matrix.gz" \
  -o "${ANALYSIS}/Fig_3.3_multisubtype_heatmap.pdf" \
  --colorMap YlOrRd \
  --zMin 0 \
  --samplesLabel "MCF7 H3K27ac" "MCF7 H3K4me1" "MCF7 H3K4me3" "MCF7 PolII" "MCF7 p300" "SKBR3 H3K27ac" "MB231 H3K27ac" "Hs578T H3K27ac" \
  --regionsLabel "Prognostic eRNAs" \
  --heatmapHeight 8 \
  --heatmapWidth 10 \
  --plotTitle ""
echo "  ✓ v1 heatmap done"

plotProfile \
  -m "${ANALYSIS}/multisubtype_matrix.gz" \
  -o "${ANALYSIS}/Fig_3.3_multisubtype_profile.pdf" \
  --samplesLabel "MCF7 H3K27ac" "MCF7 H3K4me1" "MCF7 H3K4me3" "MCF7 PolII" "MCF7 p300" "SKBR3 H3K27ac" "MB231 H3K27ac" "Hs578T H3K27ac" \
  --perGroup \
  --plotType se \
  --plotTitle ""
echo "  ✓ v1 profile done"

# =========================================================
# v2: 9-sample reference-point (center ±3kb)
# =========================================================
echo ""
echo "=============================================================="
echo " v2/v3: 9-sample reference-point"
echo "=============================================================="

computeMatrix reference-point \
  --referencePoint center \
  -R "$ERNA_BED_1kb" \
  -S ${MCF7_K27} ${MCF7_K4me1} ${MCF7_K4me3} ${MCF7_POLII} ${MCF7_P300} ${SKBR3_K27} ${MB453_K27} ${MB231_K27} ${HS578T_K27} \
  -b 3000 -a 3000 \
  --binSize 100 \
  --missingDataAsZero \
  -o "${ANALYSIS}/multisubtype_matrix_v2.gz" \
  -p 8
echo "  ✓ v2 computeMatrix done"

# v2: uniform zMax across all columns
plotHeatmap \
  -m "${ANALYSIS}/multisubtype_matrix_v2.gz" \
  -o "${ANALYSIS}/Fig_3.3_multisubtype_heatmap_v2.pdf" \
  --colorMap YlOrRd \
  --zMin 0 \
  --samplesLabel "MCF7 H3K27ac" "MCF7 H3K4me1" "MCF7 H3K4me3" "MCF7 PolII" "MCF7 p300" "SKBR3 H3K27ac" "MB453 H3K27ac" "MB231 H3K27ac" "Hs578T H3K27ac" \
  --heatmapHeight 8 \
  --heatmapWidth 12 \
  --plotTitle ""
echo "  ✓ v2 heatmap done"

# v3: per-column independent zMax
plotHeatmap \
  -m "${ANALYSIS}/multisubtype_matrix_v2.gz" \
  -o "${ANALYSIS}/Fig_3.3_multisubtype_heatmap_v3.pdf" \
  --colorMap YlOrRd \
  --zMin 0 \
  --samplesLabel "MCF7 H3K27ac" "MCF7 H3K4me1" "MCF7 H3K4me3" "MCF7 PolII" "MCF7 p300" "SKBR3 H3K27ac" "MB453 H3K27ac" "MB231 H3K27ac" "Hs578T H3K27ac" \
  --heatmapHeight 8 \
  --heatmapWidth 12 \
  --zMax 0.5 0.5 0.1 0.2 0.5 0.5 0.1 0.2 0.5 \
  --plotTitle ""
echo "  ✓ v3 heatmap done (per-column zMax)"

# Profile v2 and v3
plotProfile \
  -m "${ANALYSIS}/multisubtype_matrix_v2.gz" \
  -o "${ANALYSIS}/Fig_3.3_multisubtype_profile_v2.pdf" \
  --samplesLabel "MCF7 H3K27ac" "MCF7 H3K4me1" "MCF7 H3K4me3" "MCF7 PolII" "MCF7 p300" "SKBR3 H3K27ac" "MB453 H3K27ac" "MB231 H3K27ac" "Hs578T H3K27ac" \
  --perGroup \
  --plotType se \
  --plotTitle ""
echo "  ✓ v2 profile done"

plotProfile \
  -m "${ANALYSIS}/multisubtype_matrix_v2.gz" \
  -o "${ANALYSIS}/Fig_3.3_multisubtype_profile_v3.pdf" \
  --samplesLabel "MCF7 H3K27ac" "MCF7 H3K4me1" "MCF7 H3K4me3" "MCF7 PolII" "MCF7 p300" "SKBR3 H3K27ac" "MB453 H3K27ac" "MB231 H3K27ac" "Hs578T H3K27ac" \
  --perGroup \
  --plotType se \
  --colors red green blue purple orange cyan pink brown gray \
  --yMax 0.5 \
  --plotTitle ""
echo "  ✓ v3 profile done (yMax=0.5)"

# =========================================================
# SUMMARY
# =========================================================
END_TIME=$(date '+%Y-%m-%d %H:%M:%S')
echo ""
echo "=============================================================="
echo " 3.3 RECOMPUTE COMPLETE"
echo " Start: $START_TIME"
echo " End:   $END_TIME"
echo "=============================================================="
echo ""
echo "Generated files:"
ls -lh "${ANALYSIS}/Fig_3.3_"*.pdf 2>/dev/null
