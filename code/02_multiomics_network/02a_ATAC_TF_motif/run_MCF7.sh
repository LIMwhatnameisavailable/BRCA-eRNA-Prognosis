#!/bin/bash
# =========================================================
# MCF7 ATAC-seq Motif Analysis Pipeline (GSE298769)
# Rebuild peaks from bigWig and perform TF motif enrichment analysis
# =========================================================
set -e
START_TIME=$(date '+%Y-%m-%d %H:%M:%S')

BASE=".."
DATA="${BASE}/data"
ANALYSIS="${BASE}/analysis"
REF="../../shared_reference"
BW="${DATA}/GSE298769_Veh_merged.bw"
ERNA_BED="${REF}/prognostic_eRNA_10_2kb.bed"
HOMER_GENOME="hg38"
HOMER_DIR="../tools/homer"
FASTA="${HOMER_DIR}/data/genomes/hg38/genome.fa"
LOG="${BASE}/analysis/run_MCF7.log"
# conda environment needs to be activated manually before running

exec > >(tee -a "$LOG") 2>&1

echo "=============================================================="
echo " MCF7 ATAC-seq Motif Analysis Pipeline (GSE298769)"
echo " Started: $START_TIME"
echo "=============================================================="

# Activate conda environment (adjust path as needed)
# eval "$CONDA_ACTIVATE"

# Verify tools
echo "[CHECK] Checking required tools..."
for cmd in bigWigToBedGraph macs2 bedtools perl; do
    if which $cmd &>/dev/null; then
        echo "  ✅ $cmd found at $(which $cmd)"
    else
        echo "  ❌ $cmd NOT found"
        exit 1
    fi
done

echo "[CHECK] HOMER hg38 genome..."
if [ -f "$FASTA" ]; then
    echo "  ✅ hg38 FASTA found ($(du -h "$FASTA" | cut -f1))"
else
    echo "  ❌ hg38 FASTA not found at $FASTA"
    exit 1
fi

echo "[CHECK] eRNA BED..."
if [ -f "$ERNA_BED" ]; then
    echo "  ✅ eRNA BED found: $(wc -l < "$ERNA_BED") regions"
else
    echo "  ❌ eRNA BED not found"
    exit 1
fi

# =========================================================
# Step 1: bigWig → bedGraph
# =========================================================
echo ""
echo "=============================================================="
echo " Step 1: bigWig → bedGraph"
echo "=============================================================="
BEDGRAPH="${ANALYSIS}/GSE298769_Veh.bedGraph"

if [ -f "$BEDGRAPH" ]; then
    echo "  ⏩ bedGraph already exists, skipping"
else
    echo "  Running bigWigToBedGraph..."
    bigWigToBedGraph "$BW" "$BEDGRAPH"
    echo "  ✅ Output: $BEDGRAPH"
    echo "  Lines: $(wc -l < "$BEDGRAPH")"
    echo "  Size: $(du -h "$BEDGRAPH" | cut -f1)"
fi

# =========================================================
# Step 2: macs2 bdgpeakcall — bedGraph → peaks
# =========================================================
echo ""
echo "=============================================================="
echo " Step 2: macs2 bdgpeakcall (cutoff analysis + peak calling)"
echo "=============================================================="

# First run cutoff analysis to determine thresholds
CUTOFF_ANALYSIS="${ANALYSIS}/bdgpeak_cutoff_analysis.txt"
PEAKS_NARROW="${ANALYSIS}/ATAC_peaks.narrowPeak"
PEAKS_BED="${ANALYSIS}/ATAC_peaks.bed"

# Check signal distribution
echo "  Analyzing signal distribution..."
SORTED_BG="${ANALYSIS}/GSE298769_Veh_sorted.bedGraph"
awk '{print $4}' "$BEDGRAPH" | sort -n > "${ANALYSIS}/scores_sorted.txt"
TOTAL=$(wc -l < "${ANALYSIS}/scores_sorted.txt")
P90=$(awk -v t="$TOTAL" 'NR==int(t*0.9) {print $1}' "${ANALYSIS}/scores_sorted.txt")
P95=$(awk -v t="$TOTAL" 'NR==int(t*0.95) {print $1}' "${ANALYSIS}/scores_sorted.txt")
P99=$(awk -v t="$TOTAL" 'NR==int(t*0.99) {print $1}' "${ANALYSIS}/scores_sorted.txt")
MEAN=$(awk '{sum+=$1} END{printf "%.2f", sum/NR}' "${ANALYSIS}/scores_sorted.txt")
echo "  Score distribution:"
echo "    Mean: $MEAN"
echo "    P90:  $P90"
echo "    P95:  $P95"
echo "    P99:  $P99"

# Choose cutoff: start at P90 (rounded up to nearest integer)
if (( $(echo "$P90 > 2" | bc -l) )); then
    CUTOFF=$(echo "$P90" | awk '{printf "%.1f", $1}')
elif (( $(echo "$P95 > 1" | bc -l) )); then
    CUTOFF=2
else
    CUTOFF=2
fi
echo "  Selected cutoff: $CUTOFF"

if [ -f "$PEAKS_NARROW" ]; then
    echo "  ⏩ Peaks already exist, skipping"
else
    echo "  Running macs2 bdgpeakcall with cutoff=${CUTOFF}..."
    macs2 bdgpeakcall \
        -i "$BEDGRAPH" \
        -c "$CUTOFF" \
        -l 100 \
        -g 100 \
        --cutoff-analysis \
        -o "$PEAKS_NARROW" \
        2>&1 | tee "${ANALYSIS}/bdgpeakcall.log"
fi

if [ -f "$PEAKS_NARROW" ] && [ -s "$PEAKS_NARROW" ]; then
    PEAK_COUNT=$(wc -l < "$PEAKS_NARROW")
    echo "  ✅ Peaks called: $PEAK_COUNT"
    # Convert to BED3
    cut -f1-3 "$PEAKS_NARROW" > "$PEAKS_BED"
else
    echo "  ⚠ No peaks with cutoff=$CUTOFF, retrying with lower cutoff..."
    # Retry with cutoff 1
    macs2 bdgpeakcall \
        -i "$BEDGRAPH" \
        -c 1 \
        -l 100 \
        -g 100 \
        --cutoff-analysis \
        -o "$PEAKS_NARROW" \
        2>&1 | tee "${ANALYSIS}/bdgpeakcall_retry.log"
    if [ -f "$PEAKS_NARROW" ] && [ -s "$PEAKS_NARROW" ]; then
        PEAK_COUNT=$(wc -l < "$PEAKS_NARROW")
        echo "  ✅ Peaks called with cutoff=1: $PEAK_COUNT"
        cut -f1-3 "$PEAKS_NARROW" > "$PEAKS_BED"
    else
        echo "  ❌ No peaks called even with cutoff=1"
        echo "  Creating test peaks via simple threshold..." 
        awk -v cutoff=1 '$4 >= cutoff' "$BEDGRAPH" | mergeBed -i stdin > "$PEAKS_BED" 2>/dev/null || \
        awk -v cutoff=1 '$4 >= cutoff' "$BEDGRAPH" > "${ANALYSIS}/peaks_raw.bed"
        PEAK_COUNT=$(wc -l < "$PEAKS_BED" 2>/dev/null || echo "0")
        echo "  ✅ Merged regions above cutoff=1: $PEAK_COUNT"
    fi
fi

# =========================================================
# Step 3: bedtools intersect — ATAC peaks ∩ eRNA loci
# =========================================================
echo ""
echo "=============================================================="
echo " Step 3: bedtools intersect (ATAC peaks ∩ eRNA loci)"
echo "=============================================================="
INTERSECT_BED="${ANALYSIS}/ATAC_eRNA_intersect.bed"

if [ -f "$INTERSECT_BED" ] && [ -s "$INTERSECT_BED" ]; then
    echo "  ⏩ Intersect BED already exists, skipping"
else
    echo "  Running bedtools intersect..."
    bedtools intersect -a "$PEAKS_BED" -b "$ERNA_BED" -wa -wb > "$INTERSECT_BED" 2>&1
    echo "  ✅ Intersect regions: $(wc -l < "$INTERSECT_BED")"
fi

# Also create a clean BED (no overlap info) for HOMER
INTERSECT_CLEAN="${ANALYSIS}/ATAC_eRNA_intersect_clean.bed"
cut -f1-3 "$INTERSECT_BED" | sort -u > "$INTERSECT_CLEAN"
echo "  Clean intersect regions (unique): $(wc -l < "$INTERSECT_CLEAN")"

# =========================================================
# Step 4: Extract FASTA from intersect regions
# =========================================================
echo ""
echo "=============================================================="
echo " Step 4: bedtools getfasta (extract sequences)"
echo "=============================================================="
INTERSECT_FASTA="${ANALYSIS}/eRNA_intersect.fa"

if [ -f "$INTERSECT_FASTA" ] && [ -s "$INTERSECT_FASTA" ]; then
    echo "  ⏩ FASTA already exists, skipping"
else
    if [ -s "$INTERSECT_CLEAN" ]; then
        echo "  Running bedtools getfasta..."
        bedtools getfasta -fi "$FASTA" -bed "$INTERSECT_CLEAN" -fo "$INTERSECT_FASTA" -name 2>&1
        echo "  ✅ FASTA generated: $(grep -c '^>' "$INTERSECT_FASTA") sequences"
    else
        echo "  ⚠ No intersect regions to extract FASTA from"
    fi
fi

# =========================================================
# Step 5: HOMER findMotifsGenome.pl (intersect regions)
# =========================================================
echo ""
echo "=============================================================="
echo " Step 5: HOMER known motif enrichment (ATAC ∩ eRNA)"
echo "=============================================================="
HOMER_INTERSECT="${ANALYSIS}/homer_ATAC_eRNA"

if [ -f "${HOMER_INTERSECT}/knownResults.txt" ]; then
    echo "  ⏩ HOMER results already exist, skipping"
else
    if [ -s "$INTERSECT_CLEAN" ]; then
        echo "  Running findMotifsGenome.pl on ATAC∩eRNA regions..."
        findMotifsGenome.pl "$INTERSECT_CLEAN" "$HOMER_GENOME" "$HOMER_INTERSECT" -size given -mask 2>&1
        echo "  ✅ HOMER (ATAC∩eRNA) completed"
    else
        echo "  ⚠ No intersect regions, skipping HOMER for intersect"
    fi
fi

# =========================================================
# Step 6: HOMER findMotifsGenome.pl (all ATAC peaks, background)
# =========================================================
echo ""
echo "=============================================================="
echo " Step 6: HOMER known motif enrichment (all ATAC peaks)"
echo "=============================================================="
HOMER_ALL="${ANALYSIS}/homer_ATAC_all"

if [ -f "${HOMER_ALL}/knownResults.txt" ]; then
    echo "  ⏩ HOMER results already exist, skipping"
else
    if [ -s "$PEAKS_BED" ]; then
        echo "  Running findMotifsGenome.pl on all ATAC peaks..."
        findMotifsGenome.pl "$PEAKS_BED" "$HOMER_GENOME" "$HOMER_ALL" -size given -mask 2>&1
        echo "  ✅ HOMER (all ATAC peaks) completed"
    else
        echo "  ⚠ No peaks, skipping HOMER for all peaks"
    fi
fi

# =========================================================
# Results summary
# =========================================================
echo ""
echo "=============================================================="
echo " RESULTS SUMMARY"
echo "=============================================================="
END_TIME=$(date '+%Y-%m-%d %H:%M:%S')
echo "Start: $START_TIME"
echo "End:   $END_TIME"

echo ""
echo "--- ATAC peaks ---"
if [ -f "$PEAKS_NARROW" ] && [ -s "$PEAKS_NARROW" ]; then
    echo "Total ATAC peaks: $(wc -l < "$PEAKS_NARROW")"
else
    echo "ATAC peaks: NONE"
fi

echo ""
echo "--- ATAC ∩ eRNA ---"
if [ -f "$INTERSECT_CLEAN" ] && [ -s "$INTERSECT_CLEAN" ]; then
    echo "Intersect regions: $(wc -l < "$INTERSECT_CLEAN")"
    echo "Intersect regions (with eRNA info): $(wc -l < "$INTERSECT_BED")"
else
    echo "Intersect regions: NONE"
fi

echo ""
echo "--- Top 10 known motifs (ATAC ∩ eRNA) ---"
if [ -f "${HOMER_INTERSECT}/knownResults.txt" ]; then
    head -12 "${HOMER_INTERSECT}/knownResults.txt" | column -t
else
    echo "No HOMER results for ATAC∩eRNA"
fi

echo ""
echo "--- Top 10 known motifs (all ATAC peaks) ---"
if [ -f "${HOMER_ALL}/knownResults.txt" ]; then
    head -12 "${HOMER_ALL}/knownResults.txt" | column -t
else
    echo "No HOMER results for all ATAC peaks"
fi

echo ""
echo "=============================================================="
echo " PIPELINE COMPLETE"
echo "=============================================================="
