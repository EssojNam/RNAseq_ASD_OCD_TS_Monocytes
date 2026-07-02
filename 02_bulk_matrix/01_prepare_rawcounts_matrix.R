############################################
# 01_prepare_rawcounts_matrix.R
# Bulk RNA-seq raw count matrix preparation (featureCounts -> merged gene-level matrix)
# - Merge featureCounts outputs across samples
# - Remove Ensembl version suffix
# - Map Ensembl IDs to gene symbols
# - Collapse duplicated symbols by summing counts
############################################

############################################
## 0. Configuration and libraries
############################################
counts_dir  <- "/home/jose/Desktop/RNAseq/results/counts"
out_dir     <- "/home/jose/Desktop/RNAseq/results/counts/bulk"
out_file   <- file.path(out_dir, "rawcounts.csv")
pattern_fc  <- "\\.txt$"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

library(dplyr)
library(org.Hs.eg.db)
library(AnnotationDbi)
library(tibble) 

############################################
## 1. Read and merge featureCounts outputs
############################################
# Function to read the file and return sampple columns
read_featurecounts <- function(path) {
  tmp <- read.table(
    path,
    header = TRUE,
    sep = "\t",
    comment.char = "#",
    stringsAsFactors = FALSE
  )
  

  counts   <- tmp[, 7:ncol(tmp)]
  rownames(counts) <- tmp$Geneid
  
  #clean sample names
  colnames(counts) <- gsub(".sorted.bam$", "", colnames(counts))
  colnames(counts) <- gsub("results.alignment.", "", colnames(counts))
  
  return (counts)
}

# Read files and merge
fc_files <- list.files(counts_dir, pattern = pattern_fc, full.names = TRUE)
stopifnot(length(fc_files) > 0)

count_list <- lapply(fc_files, read_featurecounts)

# Keep common genes across files
common_genes <- Reduce(intersect, lapply(count_list, rownames)) 
raw_counts <- do.call(cbind, lapply(count_list, `[`, common_genes, , drop = FALSE))

# Remove Ensembl version numbers
ens_ids_nover <- gsub("\\.\\d+$", "", rownames(raw_counts))
rownames(raw_counts) <- ens_ids_nover             

############################################
## 2. Map Ensembl IDs to gene symbols and collapse to gene-level counts
############################################
symbols <- AnnotationDbi::mapIds(
  x = org.Hs.eg.db,
  keys = ens_ids_nover,
  column = "SYMBOL",
  keytype = "ENSEMBL",
  multiVals = "first"
) 

# Use symbol if not keep Ensembl ID
gene_labels <- ifelse(is.na(symbols), ens_ids_nover, symbols) 

df <- as.data.frame(raw_counts, check.names = FALSE) 
df$symbol <- gene_labels

df_collapsed <- df %>%
  dplyr::group_by(symbol) %>%
  dplyr::summarise(
    dplyr::across(where(is.numeric), ~ sum(.x, na.rm = TRUE)),
    .groups = "drop"
  )

raw_counts_gene <- df_collapsed %>%
  tibble::column_to_rownames(var = "symbol") %>%
  as.matrix()
dim(raw_counts_gene)
head(raw_counts_gene)

############################################
## 3. Save outputs
############################################
write.csv(raw_counts_gene, file = out_file, row.names = TRUE)
