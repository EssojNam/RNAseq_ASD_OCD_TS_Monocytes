#!/usr/bin/env bash
set -euo pipefail

# ============================================================
# Step 4: Download human reference genome and build HISAT2 index
# Genome: Homo Sapiens GRCh38 (hg38)
# ============================================================

GENOME_DIR="data/genome"
INDEX_DIR="data/genome_index"
THREADS=12

mkdir -p "${GENOME_DIR}"
mkdir -p "${INDEX_DIR}"

cd "${GENOME_DIR}"

echo "Downloading human reference genome (GRCh38)..."
wget -c https://ftp.ensembl.org/pub/release-115/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
gunzip -f Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz 

echo "Downloading gene annotation (GTF)..."
wget -c https://ftp.ensembl.org/pub/release-115/gtf/homo_sapiens/Homo_sapiens.GRCh38.115.gtf.gz
gunzip -f Homo_sapiens.GRCh38.115.gtf.gz

echo "Building HISAT2 index..."
hisat2-build -p "${THREADS}" \
  Homo_sapiens.GRCh38.dna.primary_assembly.fa \
  "../genome_index/hg38"

echo "Genome indexing completed."

