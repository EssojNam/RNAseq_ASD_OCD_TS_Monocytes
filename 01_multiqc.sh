#!/usr/bin/env bash
set -euo pipefail

# ============================================================
# Step 1.2: MultiQC summary for initial FastQC output
# Library prep: QIAseq FastSelect RNA Library Kit
# Input:  results/fastqc/raw   
# Output: results/multiqc/raw
# ============================================================

RAW_MULTIQC_DIR="results/fastqc/raw"
MULTIQC_DIR="results/multiqc/raw"

# Create output directory
mkdir -p "${MULTIQC_DIR}"

echo "Running multiqc on all results..."

multiqc "${RAW_MULTIQC_DIR}" -o "${MULTIQC_DIR}" 

echo "MultiQC report generated succesfully at: ${MULTIQC_DIR}"
