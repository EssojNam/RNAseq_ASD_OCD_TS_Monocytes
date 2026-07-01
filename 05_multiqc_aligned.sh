#!/usr/bin/env bash
set -euo pipefail

# ============================================================
# Step 3.2: MultiQC summary for all FastQC and Cutadapt outputs
# Library prep: QIAseq FastSelect RNA Library Kit
# Input:  results/fastqc/trimmed    /data/trimmed_fastq/
# Output: results/multiqc
# ============================================================

ALIG_DIR="results/alignment"
LOG_DIR="results/logs/alignment"
MULTIQC_DIR="results/multiqc/alignment"

# Create output directory
mkdir -p "${MULTIQC_DIR}"

echo "Running multiqc on all results..."

multiqc "${ALIG_DIR}" "${LOG_DIR}" -o "${MULTIQC_DIR}" 

echo "MultiQC report generated succesfully at: ${MULTIQC_DIR}"
