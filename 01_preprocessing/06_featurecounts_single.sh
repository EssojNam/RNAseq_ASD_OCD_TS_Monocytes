#!/usr/bin/env bash
set -euo pipefail

# ============================================================
# Step 6: Gene-level quantification with featureCounts (single-end)
# Input:  results/alignment/*.sorted.bam
# Output: results/counts/gene_counts.txt (+ .summary)
# Params (as used in the project): -O -t exon -g gene_id
# ============================================================

BAM_DIR="results/alignment"
COUNT_DIR="results/counts"
GTF="data/genome/Homo_sapiens.GRCh38.115.gtf"
LOG_DIR="results/logs/counts"
THREADS=4

mkdir -p "${COUNT_DIR}"
mkdir -p "${LOG_DIR}"

LOG_FILE="${LOG_DIR}/featurecounts.log"
OUTPUT="${COUNT_DIR}/gene_counts.txt"

# Collect BAM files
BAMS=( "${BAM_DIR}"/*.sorted.bam )

if [[ ${#BAMS[@]} -eq 0 ]]; then
  echo "ERROR: No BAM files found in ${BAM_DIR} (*.sorted.bam)." | tee -a "${LOG_FILE}"
  exit 1
fi

if [[ ! -f "${GTF}" ]]; then
  echo "ERROR: GTF file not found: ${GTF}" | tee -a "${LOG_FILE}"
  exit 1
fi

echo "Running featureCounts..." | tee -a "${LOG_FILE}"
echo "BAM dir: ${BAM_DIR}" | tee -a "${LOG_FILE}"
echo "GTF:     ${GTF}" | tee -a "${LOG_FILE}"
echo "Output:  ${OUTPUT}" | tee -a "${LOG_FILE}"
echo "" | tee -a "${LOG_FILE}"

featureCounts \
  -T "${THREADS}" \
  -O \
  -s 0 \
  -t exon \
  -g gene_id \
  -a "${GTF}" \
  -o "${COUNT_DIR}/gene_counts.txt" \
  "${BAMS[@]}"\
  2>&1 | tee -a "${LOG_FILE}"

echo "" | tee -a "${LOG_FILE}"
echo "featureCounts completed." | tee -a "${LOG_FILE}"
echo "Counts:   ${COUNT_DIR}/gene_counts.txt" | tee -a "${LOG_FILE}"
echo "Summary:  ${COUNT_DIR}/gene_counts.txt.summary" | tee -a "${LOG_FILE}"
echo "Log file: ${LOG_FILE}" | tee -a "${LOG_FILE}"
echo "Finished at: $(date)" | tee -a "${LOG_FILE}"

