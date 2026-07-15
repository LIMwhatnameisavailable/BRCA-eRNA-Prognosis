#!/bin/bash
set -e
BASE="../.."
GEO_BASE="https://ftp.ncbi.nlm.nih.gov/geo/series"
LOG="${BASE}/download.log"

mkdir -p "${BASE}/2.1_ATAC_TF_motif"
mkdir -p "${BASE}/2.2_ChIPseq_TF_tracks"
mkdir -p "${BASE}/3.2_RNAseq_eRNA_signal"
mkdir -p "${BASE}/3.3_H3K27ac_PolII_subtypes"
mkdir -p "${BASE}/3.7_HiChIP_MCF7"

echo "========================================" | tee -a "$LOG"
echo "Download started: $(date)" | tee -a "$LOG"
echo "========================================" | tee -a "$LOG"

echo "[TEST] Testing network connectivity..." | tee -a "$LOG"
curl -s --max-time 15 -o /dev/null -w "HTTP %{http_code} | Speed: %{speed_download} bytes/s\n" \
  "${GEO_BASE}/GSE157nnn/GSE157381/suppl/GSE157381_RAW.tar" | tee -a "$LOG"

echo "[3.7] Downloading GSE157381 HiChIP loop BED (2.1 Mb)..." | tee -a "$LOG"
wget -c --tries=5 --timeout=60 \
  -O "${BASE}/3.7_HiChIP_MCF7/GSE157381_RAW.tar" \
  "${GEO_BASE}/GSE157nnn/GSE157381/suppl/GSE157381_RAW.tar" \
  2>&1 | tee -a "$LOG"
cd "${BASE}/3.7_HiChIP_MCF7" && tar -xf GSE157381_RAW.tar
echo "[3.7] Done. Contents:" | tee -a "$LOG"
ls -lh "${BASE}/3.7_HiChIP_MCF7/" | tee -a "$LOG"

echo "[2.1] Downloading GSE298769 ATAC-seq BW (~86 Mb)..." | tee -a "$LOG"
wget -c --tries=5 --timeout=120 \
  -O "${BASE}/2.1_ATAC_TF_motif/GSE298769_Veh_merged.bw" \
  "${GEO_BASE}/GSE298nnn/GSE298769/suppl/GSE298769_Veh_merged.bw" \
  2>&1 | tee -a "$LOG"
echo "[2.1] Done." | tee -a "$LOG"

echo "[2.2] Downloading GSE298767 ChIP-seq BW files (5 parallel)..." | tee -a "$LOG"
wget -c --tries=5 --timeout=120 \
  -O "${BASE}/2.2_ChIPseq_TF_tracks/GSE298767_ER_Veh.bw" \
  "${GEO_BASE}/GSE298nnn/GSE298767/suppl/GSE298767_ER_Veh.bw" \
  2>&1 | tee -a "$LOG" &
PID_ER=$!

wget -c --tries=5 --timeout=120 \
  -O "${BASE}/2.2_ChIPseq_TF_tracks/GSE298767_FOXA1_veh_R1_R2.bw" \
  "${GEO_BASE}/GSE298nnn/GSE298767/suppl/GSE298767_FOXA1_veh_R1_R2.bw" \
  2>&1 | tee -a "$LOG" &
PID_FOXA1=$!

wget -c --tries=5 --timeout=120 \
  -O "${BASE}/2.2_ChIPseq_TF_tracks/GSE298767_GATA3_veh_R1_R2.bw" \
  "${GEO_BASE}/GSE298nnn/GSE298767/suppl/GSE298767_GATA3_veh_R1_R2.bw" \
  2>&1 | tee -a "$LOG" &
PID_GATA3=$!

wget -c --tries=5 --timeout=120 \
  -O "${BASE}/2.2_ChIPseq_TF_tracks/GSE298767_Veh_H3K27ac_merge.bw" \
  "${GEO_BASE}/GSE298nnn/GSE298767/suppl/GSE298767_Veh_H3K27ac_merge.bw" \
  2>&1 | tee -a "$LOG" &
PID_K27=$!

wget -c --tries=5 --timeout=120 \
  -O "${BASE}/2.2_ChIPseq_TF_tracks/GSE298767_IN_Veh.bw" \
  "${GEO_BASE}/GSE298nnn/GSE298767/suppl/GSE298767_IN_Veh.bw" \
  2>&1 | tee -a "$LOG" &
PID_IN=$!

wait $PID_ER $PID_FOXA1 $PID_GATA3 $PID_K27 $PID_IN
echo "[2.2] All ChIP-seq BW done." | tee -a "$LOG"

echo "[3.2] Downloading GSE298771 RNA-seq BW (Veh, 2 files parallel)..." | tee -a "$LOG"
wget -c --tries=5 --timeout=120 \
  -O "${BASE}/3.2_RNAseq_eRNA_signal/GSE298771_Veh_Forward.bw" \
  "${GEO_BASE}/GSE298nnn/GSE298771/suppl/GSE298771_E2dosage_Veh_2H.sorted_strandSp.Forward.bw" \
  2>&1 | tee -a "$LOG" &
PID_FWD=$!

wget -c --tries=5 --timeout=120 \
  -O "${BASE}/3.2_RNAseq_eRNA_signal/GSE298771_Veh_Reverse.bw" \
  "${GEO_BASE}/GSE298nnn/GSE298771/suppl/GSE298771_E2dosage_Veh_2H.sorted_strandSp.Reverse.bw" \
  2>&1 | tee -a "$LOG" &
PID_REV=$!

wait $PID_FWD $PID_REV
echo "[3.2] RNA-seq BW done." | tee -a "$LOG"

echo "[3.3] Downloading 3 large TAR files in parallel..." | tee -a "$LOG"
wget -c --tries=5 --timeout=600 \
  -O "${BASE}/3.3_H3K27ac_PolII_subtypes/GSE251871_RAW.tar" \
  "${GEO_BASE}/GSE251nnn/GSE251871/suppl/GSE251871_RAW.tar" \
  2>&1 | tee -a "$LOG" &
PID_TNBC=$!

wget -c --tries=5 --timeout=600 \
  -O "${BASE}/3.3_H3K27ac_PolII_subtypes/GSE251868_RAW.tar" \
  "${GEO_BASE}/GSE251nnn/GSE251868/suppl/GSE251868_RAW.tar" \
  2>&1 | tee -a "$LOG" &
PID_HER2=$!

wget -c --tries=5 --timeout=600 \
  -O "${BASE}/3.3_H3K27ac_PolII_subtypes/GSE245868_RAW.tar" \
  "${GEO_BASE}/GSE245nnn/GSE245868/suppl/GSE245868_RAW.tar" \
  2>&1 | tee -a "$LOG" &
PID_POLII=$!

wait $PID_TNBC $PID_HER2 $PID_POLII
echo "[3.3] All TAR files downloaded. Extracting..." | tee -a "$LOG"

cd "${BASE}/3.3_H3K27ac_PolII_subtypes"
tar -xf GSE251871_RAW.tar && echo "[EXTRACT] GSE251871 done" | tee -a "$LOG"
tar -xf GSE251868_RAW.tar && echo "[EXTRACT] GSE251868 done" | tee -a "$LOG"
tar -xf GSE245868_RAW.tar && echo "[EXTRACT] GSE245868 done" | tee -a "$LOG"

echo "========================================" | tee -a "$LOG"
echo "ALL DONE: $(date)" | tee -a "$LOG"
echo "========================================" | tee -a "$LOG"
du -sh "${BASE}"/*/  | tee -a "$LOG"
df -h / | tee -a "$LOG"
