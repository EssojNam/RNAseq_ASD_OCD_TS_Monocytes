#!/usr/bin/env bash
set -euo pipefail

# ============================================================
# Step 3.2: MultiQC summary for all FastQC and Cutadapt outputs
# Library prep: QIAseq FastSelect RNA Library Kit
# Input:  results/fastqc/trimmed    /data/trimmed_fastq/
# Output: results/multiqc
# ============================================================

RAW_TRIM_DIR="data/trimmed_fastq"
RAW_MULTIQC_DIR="results/fastqc/trimmed"
MULTIQC_DIR="results/multiqc"

# Create output directory
mkdir -p "${MULTIQC_DIR}"

echo "Running multiqc on all results..."

multiqc "${RAW_MULTIQC_DIR}" "${RAW_TRIM_DIR}" -o "${MULTIQC_DIR}" 

echo "MultiQC report generated succesfully at: ${MULTIQC_DIR}"
