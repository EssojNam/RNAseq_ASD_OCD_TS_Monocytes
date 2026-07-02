# RNAseq_ASD_OCD_TS_Monocytes

**Monocytes transcriptional responses in Basal or LPS-stimulated conditions in three neuropsychiatric disorders (ASD, OCD, TS)**

This repository contains the full analysis pipeline to characterize gene expression profiles of monocytes in basal or LPS conditions in patients with Autism Spectrum Disorder (ASD), Obsessive-Compulsive Disorder (OCD), and Tourette Syndrome (TS). The study involves bulk RNA-seq data  and applies quality control and alignment of raw data, followed by differential gene expression and gene set enrichment analyses.

---

## 📚 Project Overview

The goal of this project is to investigate the biological mechanisms underlying immune stimulation in neuropsychiatric disorders. To achieve this, we characterized monocytes transcriptional responses after an immune stimulation in neurodevelopmental and neuropsychiatric disorders. The study uses bulk RNA-seq data from immuno-stilated monocytes populations and applies:
- **Differential gene expression (DEG) analysis**: to identify genes uniquely regulated after an before treatment for each diagnosis and pairwise comparisons.
- **Gene Set Enrichment Analysis (GSEA)**: to identify molecular programs associated with diagnosis-specific basal and LPS-stimulated profiles.

**Study Design**:
- **DisordersCell types**: Autism Spectrum Disorder (ASD), Obsessive-Compulsive Disorder (OCD), and Tourette Syndrome (TS).
- **Treatment**: Basal or LPS-stimulated conditions.
- **Analyses**: DEG analysis and GSEA/Hallmark, KEGG, Reactome, and GO.

---

## 🧪 Analysis Pipeline

The analysis is structured into sequential modules:

| Folder | Description |
|--------|-------------|
| `0_setup/` | Installation of required command-line tools, version checks, and environment setup|
| `00_Descriptive_data/`| Statistical analysis of clinical data and monocytes population
| `1_preprocessing/` | Raw FASTQ files preprocessing: quality control (QC), adapter trimming, alignment to the human genome (GRCh38), and gene-level quantification|
| `2_bulk_matrix/` | Generation of a clean bulk gene expression matrix, handling batch effects (e.g., ComBat-seq) and quality control (PCA) before downstream analyses|
| `3_dataAnalysis/` | Core statistical analysis: DEG (DESeq2), GSEA (clusterProfiler) and KEGG and Reactome categorization|

---

## 🛠 Setup (required tools and packages)

### Command-line tools: 
The preprocessing pipeline requires the following tools (installable via Linux):
- FastQC
- cutadapt
- HISAT2
- samtools
- subread (featureCounts)

### R packages:

- **Core RNA-seq analysis**:  
  `DESeq2`, `sva`, `BiocGenerics`

- **Functional enrichment and GO analysis**:  
  `clusterProfiler`, `org.Mm.eg.db`, `AnnotationDbi`, `GO.db`
  
- **Data manipulation**:  
  `dplyr`, `tidyr`, `readr`, `tibble`, `stringr`, `purrr`, `reshape2`, `tidyverse`

- **Visualization**:  
  `ggplot2`, `pheatmap`, `RColorBrewer`, `ggrepel`, `VennDiagram`, `ggvenn`, `grid`

- **Input / Output utilities**:  
  `openxlsx`
  
The analysis depends on several Bioconductor and CRAN packages. You can install all dependencies by running:

```r
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
# Bioconductor packages
bioc_pkgs <- c("DESeq2", "sva", "BiocGenerics", "clusterProfiler", "org.Mm.eg.db", 
               "AnnotationDbi", "GO.db", 
               "rhdf5", "biomaRt")
# CRAN packages
cran_pkgs <- c(tidyverse", "pheatmap", "RColorBrewer", "ggrepel", 
               "VennDiagram", "ggvenn", "openxlsx", "reshape2")
BiocManager::install(bioc_pkgs)
install.packages(cran_pkgs)
```
---

## ⚙️ Reproducibility

To reproduce the analysis:

1. Clone the repository:
   ```bash
   git clone https://github.com/EssojNam/RNAseq_ASD_OCD_TS_Monocytes.git
   ```

2. Set up the environment: Run the scripts in 0_setup/ to ensure all command-line tools are available.


3. Run the scripts in the order specified above.

