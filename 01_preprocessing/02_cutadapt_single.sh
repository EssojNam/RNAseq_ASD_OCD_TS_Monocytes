#!/usr/bin/env bash
set -euo pipefail

# ============================================================
# Step 2: Adapter trimming with cutadapt (single-end RNA-seq)
# Library prep: QIAseq FastSelect RNA Library Kit
# Adapters: Illumina TruSeq
# Input:  data/raw_fastq/*_R1_001.fastq.gz
# Output: data/trimmed_fastq/
# ============================================================

RAW_FASTQ_DIR="data/raw_fastq"
TRIMMED_FASTQ_DIR="data/trimmed_fastq"
THREADS=8

# Illumina TruSeq adapter sequences
ADAPTER_R1="AGATCGGAAGAGCACACGTCTGAACTCCAGTCA"

# Create output directory
mkdir -p "${TRIMMED_FASTQ_DIR}"

echo "Running cutadapt for single-end FASTQ files..."

for R1 in "${RAW_FASTQ_DIR}"/*_R1_001.fastq.gz; do

  SAMPLE_NAME=$(basename "${R1}" _R1_001.fastq.gz)

  OUT_R1="${TRIMMED_FASTQ_DIR}/${SAMPLE_NAME}_R1_trimmed.fastq.gz"

  echo "Trimming adapters for sample: ${SAMPLE_NAME}"

  cutadapt \
    -j "${THREADS}" \
    -a "${ADAPTER_R1}" \
    -m 20 \
    -o "${OUT_R1}" \
    "${R1}"

done

echo "Cutadapt trimming completed successfully."

