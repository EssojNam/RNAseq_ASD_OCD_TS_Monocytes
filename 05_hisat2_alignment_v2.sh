#!/usr/bin/env bash
set -euo pipefail
# ============================================================
# Step 5: Single-end alignment with HISAT2 (parallelized)
# Input:  data/trimmed_fastq/*_R1_trimmed.fastq.gz
# Output: results/alignment/*.sorted.bam (+ .bai + .flagstat.txt)
# Requires: hisat2 (conda), samtools (conda)
# ============================================================

# ── Paths ─────────────────────────────────────────────────
INDEX_PERMANENT="data/genome_index_v2/hg38"   # pool01 HDD (master copy)
INDEX_TMP="/tmp/genome_index_v2/hg38"          # NVMe (fast working copy)
INPUT_DIR="data/trimmed_fastq"
OUTPUT_DIR="results/alignment"
LOG_DIR="results/logs/alignment"
PROGRESS_LOG="results/logs/alignment_progress.log"

# ── Concurrency settings ─────────────────────────────────────
PARALLEL_JOBS=3        # samples processed simultaneously
THREADS=8              # threads per sample (3 × 8 = 24 total CPUs used)

# ── Setup ────────────────────────────────────────────────────
mkdir -p "${OUTPUT_DIR}" "${LOG_DIR}" "$(dirname ${PROGRESS_LOG})"

# ── Verify conda tools are being used ────────────────────────
HISAT2_PATH=$(which hisat2)
SAMTOOLS_PATH=$(which samtools)

if [[ "${HISAT2_PATH}" != *"miniconda"* ]]; then
    echo "ERROR: hisat2 is not from conda: ${HISAT2_PATH}"
    echo "Run: export PATH=/opt/miniconda3/bin:\$PATH"
    exit 1
fi
echo "hisat2:   ${HISAT2_PATH}"
echo "samtools: ${SAMTOOLS_PATH}"
echo ""

# ── Copy index to NVMe if not already there ──────────────────
if [ ! -f "${INDEX_TMP}.1.ht2" ]; then
    echo "[$(date '+%H:%M:%S')] Copying genome index to NVMe /tmp..."
    mkdir -p /tmp/genome_index_v2
    cp /home/jose/Desktop/RNAseq/data/genome_index_v2/*.ht2 /tmp/genome_index_v2/
    echo "[$(date '+%H:%M:%S')] Index ready."
else
    echo "[$(date '+%H:%M:%S')] Index already in /tmp, skipping copy."
fi
echo ""

# ── Count total samples ──────────────────────────────────────
TOTAL=$(ls "${INPUT_DIR}"/*_R1_trimmed.fastq.gz | wc -l)
echo "Total samples to align: ${TOTAL}"
echo "Parallel jobs:          ${PARALLEL_JOBS}"
echo "Threads per sample:     ${THREADS}"
echo "Estimated total time:   ~$(( (TOTAL / PARALLEL_JOBS) * 3 )) minutes"
echo ""
echo "Started at: $(date)" | tee -a "${PROGRESS_LOG}"
echo "────────────────────────────────────────────" | tee -a "${PROGRESS_LOG}"

# ── Per-sample alignment function ────────────────────────────
align_sample() {
    local R1="$1"
    local INDEX_TMP="$2"
    local OUTPUT_DIR="$3"
    local LOG_DIR="$4"
    local THREADS="$5"
    local PROGRESS_LOG="$6"

    local SAMPLE
    SAMPLE=$(basename "${R1}" _R1_trimmed.fastq.gz)
    local BAM_OUT="${OUTPUT_DIR}/${SAMPLE}.sorted.bam"

    # ── Resume: skip if BAM + index already exist ─────────────
    if [ -f "${BAM_OUT}" ] && [ -f "${BAM_OUT}.bai" ]; then
        echo "[$(date '+%H:%M:%S')] SKIP (already done): ${SAMPLE}" | tee -a "${PROGRESS_LOG}"
        return 0
    fi

    echo "[$(date '+%H:%M:%S')] START: ${SAMPLE}" | tee -a "${PROGRESS_LOG}"

    # Copy FASTQ to NVMe for fast I/O
    local TMP_FASTQ="/tmp/${SAMPLE}_R1.fastq.gz"
    cp "${R1}" "${TMP_FASTQ}"

    # Align
    hisat2 \
        -p "${THREADS}" \
        --dta \
        -x "${INDEX_TMP}" \
        -U "${TMP_FASTQ}" \
        2>"${LOG_DIR}/${SAMPLE}.hisat2.log" \
    | samtools view -bS - \
    | samtools sort -@ "${THREADS}" \
        -o "${BAM_OUT}"

    # Index BAM
    samtools index "${BAM_OUT}"

    # Alignment summary
    samtools flagstat "${BAM_OUT}" \
        > "${OUTPUT_DIR}/${SAMPLE}.flagstat.txt"

    # Clean up NVMe temp FASTQ
    rm -f "${TMP_FASTQ}"

    # Log alignment rate from hisat2 log
    local ALIGN_RATE
    ALIGN_RATE=$(grep "overall alignment rate" "${LOG_DIR}/${SAMPLE}.hisat2.log" || echo "unknown")
    echo "[$(date '+%H:%M:%S')] DONE:  ${SAMPLE} — ${ALIGN_RATE}" | tee -a "${PROGRESS_LOG}"
}

export -f align_sample

# ── Dispatch ─────────────────────────────────────────────────
if command -v parallel &>/dev/null; then
    echo "Using GNU parallel"
    find "${INPUT_DIR}" -name '*_R1_trimmed.fastq.gz' | sort \
    | parallel \
        --jobs "${PARALLEL_JOBS}" \
        --eta \
        align_sample {} \
            "${INDEX_TMP}" \
            "${OUTPUT_DIR}" \
            "${LOG_DIR}" \
            "${THREADS}" \
            "${PROGRESS_LOG}"
else
    echo "GNU parallel not found — using bash background jobs"
    active=0
    pids=()

    for R1 in "${INPUT_DIR}"/*_R1_trimmed.fastq.gz; do
        align_sample \
            "${R1}" \
            "${INDEX_TMP}" \
            "${OUTPUT_DIR}" \
            "${LOG_DIR}" \
            "${THREADS}" \
            "${PROGRESS_LOG}" &
        pids+=($!)
        (( active++ ))

        if (( active >= PARALLEL_JOBS )); then
            wait "${pids[0]}"
            pids=("${pids[@]:1}")
            (( active-- ))
        fi
    done

    wait
fi

# ── Summary ──────────────────────────────────────────────────
echo "" | tee -a "${PROGRESS_LOG}"
echo "════════════════════════════════════════════" | tee -a "${PROGRESS_LOG}"
echo "Finished at: $(date)" | tee -a "${PROGRESS_LOG}"

DONE=$(ls "${OUTPUT_DIR}"/*.sorted.bam 2>/dev/null | wc -l)
echo "BAMs completed: ${DONE} / ${TOTAL}" | tee -a "${PROGRESS_LOG}"

if (( DONE < TOTAL )); then
    echo "WARNING: ${DONE} of ${TOTAL} samples completed. Check logs in ${LOG_DIR}" | tee -a "${PROGRESS_LOG}"
else
    echo "All samples aligned successfully." | tee -a "${PROGRESS_LOG}"
fi

echo ""
echo "BAM files : ${OUTPUT_DIR}"
echo "Logs      : ${LOG_DIR}"
echo "Progress  : ${PROGRESS_LOG}"

