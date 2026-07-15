#!/bin/bash
# =========================================================
# 3.2_RNAseq_eRNA_signal - Re-run with CORRECT hg38 coordinates
# =========================================================
set -eo pipefail

START_TIME=$(date '+%Y-%m-%d %H:%M:%S')

BASE=".."
REF="../../shared_reference"
ANALYSIS="${BASE}/analysis"
# conda environment needs to be activated manually before running

# eval "$CONDA_ACTIVATE"
mkdir -p "$ANALYSIS"
cd "$ANALYSIS"

ERNA_BED="${REF}/prognostic_eRNA_10_2kb.bed"
ERNA_BED_1kb="${REF}/prognostic_eRNA_10.bed"

# BigWig files
RNA_FWD="${BASE}/GSE298771_Veh_Forward.bw"
RNA_REV="${BASE}/GSE298771_Veh_Reverse.bw"
K27="${BASE}/../2.2_ChIPseq_TF_tracks/GSE298767_Veh_H3K27ac_merge.bw"
PRO_PLUS="${BASE}/GSE298770_MCF7_Veh_plus.bw"
PRO_MINUS="${BASE}/GSE298770_MCF7_Veh_minus.bw"
E2_PLUS="${BASE}/GSE298770_MCF7_10nM_E2_plus.bw"
E2_MINUS="${BASE}/GSE298770_MCF7_10nM_E2_minus.bw"

echo "=============================================================="
echo " 3.2 RNA-seq eRNA Signal - Recompute with hg38"
echo " Started: $START_TIME"
echo "=============================================================="

# =========================================================
# PART A: H3K27ac + RNA-seq initial analysis (scale-regions)
# =========================================================
echo ""
echo "=== A1: H3K27ac only heatmap (scale-regions) ==="
computeMatrix scale-regions \
  -R "$ERNA_BED" \
  -S "$K27" \
  --regionBodyLength 1000 \
  -b 2000 -a 2000 \
  --binSize 50 \
  --missingDataAsZero \
  -o "${ANALYSIS}/matrix_eRNA_H3K27ac.gz" \
  -p 8
echo "  ✓ computeMatrix A1"

plotHeatmap \
  -m "${ANALYSIS}/matrix_eRNA_H3K27ac.gz" \
  -o "${ANALYSIS}/Fig_3.2_H3K27ac_only_heatmap.pdf" \
  --colorMap YlOrRd \
  --zMin 0 \
  --samplesLabel "H3K27ac" \
  --heatmapHeight 8 \
  --heatmapWidth 4 \
  --plotTitle ""
echo "  ✓ H3K27ac heatmap"

echo ""
echo "=== A2: H3K27ac vs background ==="
computeMatrix scale-regions \
  -R "$ERNA_BED" "${ANALYSIS}/random_bg_2kb.bed" \
  -S "$K27" \
  --regionBodyLength 1000 \
  -b 2000 -a 2000 \
  --binSize 50 \
  --missingDataAsZero \
  -o "${ANALYSIS}/matrix_bg.gz" \
  -p 8
echo "  ✓ computeMatrix A2"

plotProfile \
  -m "${ANALYSIS}/matrix_bg.gz" \
  -o "${ANALYSIS}/Fig_3.2_H3K27ac_vs_bg_profile.pdf" \
  --perGroup \
  --plotType se \
  --colors red gray \
  --plotTitle ""
echo "  ✓ H3K27ac vs bg profile"

echo ""
echo "=== A3: Multi-track (RNA-seq + H3K27ac) ==="
computeMatrix scale-regions \
  -R "$ERNA_BED" \
  -S "$RNA_FWD" "$RNA_REV" "$K27" \
  --regionBodyLength 1000 \
  -b 2000 -a 2000 \
  --binSize 50 \
  --missingDataAsZero \
  -o "${ANALYSIS}/matrix_eRNA.gz" \
  -p 8
echo "  ✓ computeMatrix A3"

plotHeatmap \
  -m "${ANALYSIS}/matrix_eRNA.gz" \
  -o "${ANALYSIS}/Fig_3.2_heatmap.pdf" \
  --colorMap RdYlBu_r \
  --samplesLabel "RNA-seq Fwd" "RNA-seq Rev" "H3K27ac" \
  --heatmapHeight 8 \
  --heatmapWidth 6 \
  --plotTitle ""
echo "  ✓ multi-track heatmap"

plotProfile \
  -m "${ANALYSIS}/matrix_eRNA.gz" \
  -o "${ANALYSIS}/Fig_3.2_profile_vs_bg.pdf" \
  --samplesLabel "RNA-seq Fwd" "RNA-seq Rev" "H3K27ac" \
  --perGroup \
  --plotType se \
  --plotTitle ""
echo "  ✓ multi-track profile"

# =========================================================
# PART B: PRO-seq analysis (reference-point, center ±3kb)
# =========================================================
echo ""
echo "=== B1: PRO-seq DMSO ==="
computeMatrix reference-point \
  --referencePoint center \
  -R "$ERNA_BED_1kb" \
  -S "$PRO_PLUS" "$PRO_MINUS" \
  -b 3000 -a 3000 \
  --binSize 50 \
  --missingDataAsZero \
  -o "${ANALYSIS}/GSE298770_PROseq_matrix.gz" \
  -p 8
echo "  ✓ computeMatrix B1"

plotHeatmap \
  -m "${ANALYSIS}/GSE298770_PROseq_matrix.gz" \
  -o "${ANALYSIS}/Fig_3.2_GSE298770_PROseq_heatmap.pdf" \
  --colorMap YlOrRd \
  --zMin 0 \
  --samplesLabel "PRO-seq +" "PRO-seq -" \
  --heatmapHeight 8 \
  --heatmapWidth 4 \
  --plotTitle ""
echo "  ✓ PRO-seq heatmap"

plotProfile \
  -m "${ANALYSIS}/GSE298770_PROseq_matrix.gz" \
  -o "${ANALYSIS}/Fig_3.2_GSE298770_PROseq_profile.pdf" \
  --samplesLabel "PRO-seq +" "PRO-seq -" \
  --perGroup \
  --plotType se \
  --plotTitle ""
echo "  ✓ PRO-seq profile"

echo ""
echo "=== B2: PRO-seq DMSO vs E2 ==="
computeMatrix reference-point \
  --referencePoint center \
  -R "$ERNA_BED_1kb" \
  -S "$PRO_PLUS" "$PRO_MINUS" "$E2_PLUS" "$E2_MINUS" \
  -b 3000 -a 3000 \
  --binSize 50 \
  --missingDataAsZero \
  -o "${ANALYSIS}/GSE298770_PROseq_DMSOvsE2_matrix.gz" \
  -p 8

plotHeatmap \
  -m "${ANALYSIS}/GSE298770_PROseq_DMSOvsE2_matrix.gz" \
  -o "${ANALYSIS}/Fig_3.2_GSE298770_PROseq_DMSOvsE2_heatmap.pdf" \
  --colorMap YlOrRd \
  --zMin 0 \
  --samplesLabel "DMSO +" "DMSO -" "E2 +" "E2 -" \
  --heatmapHeight 8 \
  --heatmapWidth 6 \
  --plotTitle ""

plotProfile \
  -m "${ANALYSIS}/GSE298770_PROseq_DMSOvsE2_matrix.gz" \
  -o "${ANALYSIS}/Fig_3.2_GSE298770_PROseq_DMSOvsE2_profile.pdf" \
  --samplesLabel "DMSO +" "DMSO -" "E2 +" "E2 -" \
  --perGroup \
  --plotType se \
  --plotTitle ""
echo "  ✓ DMSO vs E2 done"

echo ""
echo "=== B3: Final multi-omics integrated ==="
computeMatrix reference-point \
  --referencePoint center \
  -R "$ERNA_BED_1kb" \
  -S "$RNA_FWD" "$RNA_REV" "$K27" "$PRO_PLUS" "$PRO_MINUS" "$E2_PLUS" "$E2_MINUS" \
  -b 3000 -a 3000 \
  --binSize 50 \
  --missingDataAsZero \
  -o "${ANALYSIS}/GSE298770_final_combined_matrix.gz" \
  -p 8

plotHeatmap \
  -m "${ANALYSIS}/GSE298770_final_combined_matrix.gz" \
  -o "${ANALYSIS}/Fig_3.2_final_multiomics_heatmap.pdf" \
  --colorMap YlOrRd \
  --zMin 0 \
  --samplesLabel "RNA-Fwd" "RNA-Rev" "H3K27ac" "PRO+" "PRO-" "E2+" "E2-" \
  --heatmapHeight 8 \
  --heatmapWidth 8 \
  --plotTitle ""
echo "  ✓ Multi-omics heatmap"

# =========================================================
# PART C: v4 style PRO-seq (±strand separate heatmaps)
# =========================================================
echo ""
echo "=== C1: PRO-seq plus strand (v4) ==="
computeMatrix reference-point \
  --referencePoint center \
  -R "$ERNA_BED_1kb" \
  -S "$PRO_PLUS" \
  -b 3000 -a 3000 \
  --binSize 50 \
  --missingDataAsZero \
  -o "${ANALYSIS}/matrix_plus_v2.gz" \
  -p 8

plotHeatmap \
  -m "${ANALYSIS}/matrix_plus_v2.gz" \
  -o "${ANALYSIS}/Fig_3.1_PROseq_heatmap_plus_v4.pdf" \
  --colorMap YlOrRd \
  --zMin 0 --zMax 2 \
  --samplesLabel "PRO-seq +" \
  --heatmapHeight 8 --heatmapWidth 3 \
  --sortRegions descend --sortUsing mean \
  --whatToShow 'heatmap only' \
  --plotTitle ""
echo "  ✓ PRO-seq plus heatmap"

echo ""
echo "=== C2: PRO-seq minus strand (v4) ==="
computeMatrix reference-point \
  --referencePoint center \
  -R "$ERNA_BED_1kb" \
  -S "$PRO_MINUS" \
  -b 3000 -a 3000 \
  --binSize 50 \
  --missingDataAsZero \
  --scale -1 \
  -o "${ANALYSIS}/matrix_minus_v2.gz" \
  -p 8

plotHeatmap \
  -m "${ANALYSIS}/matrix_minus_v2.gz" \
  -o "${ANALYSIS}/Fig_3.1_PROseq_heatmap_minus_v4.pdf" \
  --colorMap YlOrRd \
  --zMin 0 --zMax 2 \
  --samplesLabel "PRO-seq -" \
  --heatmapHeight 8 --heatmapWidth 3 \
  --sortRegions descend --sortUsing mean \
  --whatToShow 'heatmap only' \
  --plotTitle ""
echo "  ✓ PRO-seq minus heatmap"

echo ""
echo "=== C3: Both strands combined v4 profile ==="
computeMatrix reference-point \
  --referencePoint center \
  -R "$ERNA_BED_1kb" \
  -S "$PRO_PLUS" "$PRO_MINUS" \
  -b 3000 -a 3000 \
  --binSize 50 \
  --missingDataAsZero \
  -o "${ANALYSIS}/matrix_bothstrands_v3.gz" \
  -p 8

plotProfile \
  -m "${ANALYSIS}/matrix_bothstrands_v3.gz" \
  -o "${ANALYSIS}/Fig_3.1_PROseq_profile_v4.pdf" \
  --samplesLabel "PRO-seq +" "PRO-seq -" \
  --perGroup \
  --plotType se \
  --plotTitle "" \
  --plotHeight 4 --plotWidth 6
echo "  ✓ v4 profile done"

# =========================================================
# SUMMARY
# =========================================================
END_TIME=$(date '+%Y-%m-%d %H:%M:%S')
echo ""
echo "=============================================================="
echo " 3.2 RECOMPUTE COMPLETE"
echo " Start: $START_TIME"
echo " End:   $END_TIME"
echo "=============================================================="
echo ""
echo "Generated files:"
ls -lh "${ANALYSIS}/Fig_3.2_"*.pdf "${ANALYSIS}/Fig_3.1_"*.pdf 2>/dev/null