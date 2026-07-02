# ============================================================
# Visualización por categorías — Barplot + Boxplot
# ============================================================

library(readODS)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggbeeswarm)   # install.packages("ggbeeswarm") si no lo tienes

# ── Paths ─────────────────────────────────────────────────────
table_path <- "/home/jose/Desktop/RNAseq/results/Enrichment/GSEA_inputs/Reports/CrossModules/REACTOME/M3/Table_REACTOME_M3_comparison.ods"
plot_dir   <- "/home/jose/Desktop/RNAseq/results/plots/Enrichment/Categorization/M3/M3.REACTOME"
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

FDR_THR <- 0.05

# ── Cargar tabla ──────────────────────────────────────────────
df <- read_ods(table_path, sheet = "Table_REACTOME_M3_comparison")

# ============================================================
# 1. Pasar a formato largo por diagnóstico
# ============================================================
df_long <- df %>%
  # Solo pathways categorizados
  filter(!is.na(M3_category), M3_category != "") %>%
  pivot_longer(
    cols = c(M3_TEA_vs_TOC_NES, M3_TEA_vs_Tourette_NES, M3_TOC_vs_Tourette_NES),
    names_to  = "diagnosis_raw",
    values_to = "NES"
  ) %>%
  mutate(
    diagnosis = gsub("^M3_|_NES$", "", diagnosis_raw),
    FDR = case_when(
      diagnosis == "TEA_vs_TOC"      ~ M3_TEA_vs_TOC_FDR,
      diagnosis == "TEA_vs_Tourette"      ~ M3_TEA_vs_Tourette_FDR,
      diagnosis == "TOC_vs_Tourette" ~ M3_TOC_vs_Tourette_FDR
    ),
    sig = case_when(
      diagnosis == "TEA_vs_TOC"      ~ M3_TEA_vs_TOC_sig,
      diagnosis == "TEA_vs_Tourette"      ~ M3_TEA_vs_Tourette_sig,
      diagnosis == "TOC_vs_Tourette" ~ M3_TOC_vs_Tourette_sig
    ),
    significant = sig == "sig" & FDR <= FDR_THR,
    diagnosis   = factor(diagnosis, levels = c("TEA_vs_TOC", "TEA_vs_Tourette", "TOC_vs_Tourette"))
  ) %>%
  select(pathway, M3_category, M3_sharing, diagnosis,
         NES, FDR, significant)

# ============================================================
# 2. Tabla resumen por categoría y diagnóstico
# ============================================================
# Media de NES — solo pathways significativos
category_summary <- df_long %>%
  filter(significant) %>%
  group_by(M3_category, diagnosis) %>%
  summarise(
    mean_NES  = mean(NES, na.rm = TRUE),
    sd_NES    = sd(NES, na.rm = TRUE),
    n_sig     = n(),
    n_up      = sum(NES > 0),
    n_down    = sum(NES < 0),
    .groups   = "drop"
  )

write.csv(category_summary,
          file.path(plot_dir, "Category_summary.csv"),
          row.names = FALSE)

# ============================================================
# 3. Barplot — media NES por categoría y diagnóstico
# ============================================================
# Orden de categorías por media NES absoluta global
cat_order <- category_summary %>%
  group_by(M3_category) %>%
  summarise(mean_abs = mean(abs(mean_NES))) %>%
  arrange(desc(mean_abs)) %>%
  pull(M3_category)

p_bar <- category_summary %>%
  mutate(M3_category = factor(M3_category, levels = rev(cat_order))) %>%
  ggplot(aes(x = M3_category, y = mean_NES, fill = diagnosis)) +
  geom_col(position = position_dodge(width = 0.7),
           width = 0.65, alpha = 0.85) +
  geom_errorbar(
    aes(ymin = mean_NES - sd_NES,
        ymax = mean_NES + sd_NES),
    position = position_dodge(width = 0.7),
    width = 0.25, linewidth = 0.4, color = "grey30"
  ) +
  geom_hline(yintercept = 0, color = "grey40", linewidth = 0.4) +
  geom_text(
    aes(label = paste0("n=", n_sig),
         y = mean_NES + sign(mean_NES) * (sd_NES + 0.2)),
    position = position_dodge(width = 0.7),
    size = 2.5, color = "grey30"
  ) +
  coord_flip() +
  scale_fill_manual(values = c(
    TEA_vs_TOC      = "#9b59b6",
    TEA_vs_Tourette      = "#2ecc71",
    TOC_vs_Tourette = "#e67e22"
  )) +
  scale_y_continuous(limits = c(-3.5, 3.5)) +
  labs(
    title    = "Mean NES per category by diagnosis",
    subtitle = paste0("Only significant pathways (FDR ≤ ", FDR_THR,
                      ")  |  Error bars = SD  |  n = number of pathways"),
    x = NULL, y = "Mean NES",
    fill = "Diagnosis"
  ) +
  theme_bw(base_size = 12) +
  theme(
    plot.title      = element_text(face = "bold"),
    legend.position = "bottom"
  )

ggsave(file.path(plot_dir, "Barplot_category_NES.png"),
       p_bar, width = 10, height = 7, dpi = 300)
message("Barplot guardado.")

# ============================================================
# 4. Boxplot — distribución de NES por categoría y diagnóstico
# ============================================================
# Todos los pathways de la tabla, significativos marcados distinto

p_box <- df_long %>%
  mutate(M3_category = factor(M3_category, levels = rev(cat_order))) %>%
  ggplot(aes(x = M3_category, y = NES, color = diagnosis)) +
  geom_hline(yintercept = 0, color = "grey60",
             linewidth = 0.4, linetype = "dashed") +
  geom_boxplot(
    aes(fill = diagnosis),
    alpha    = 0.15,
    outlier.shape = NA,
    position = position_dodge(width = 0.75),
    width    = 0.6
  ) +
  geom_beeswarm(
    aes(shape = significant, size = significant),
    dodge.width = 0.75,
    alpha = 0.8,
    cex   = 2
  ) +
  scale_color_manual(values = c(
    TEA_vs_TOC      = "#9b59b6",
    TEA_vs_Tourette      = "#2ecc71",
    TOC_vs_Tourette = "#e67e22"
  )) +
  scale_fill_manual(values = c(
    TEA_vs_TOC      = "#9b59b6",
    TEA_vs_Tourette      = "#2ecc71",
    TOC_vs_Tourette = "#e67e22"
  )) +
  scale_shape_manual(values = c("TRUE" = 19, "FALSE" = 1),
                     labels = c("TRUE" = "Significant", "FALSE" = "ns")) +
  scale_size_manual(values  = c("TRUE" = 2.5, "FALSE" = 1.5),
                    labels  = c("TRUE" = "Significant", "FALSE" = "ns")) +
  coord_flip() +
  scale_y_continuous(limits = c(-3.5, 3.5)) +
  labs(
    title    = "NES distribution per category by diagnosis",
    subtitle = paste0("All selected pathways  |  Filled points = FDR ≤ ",
                      FDR_THR, "  |  Open points = ns"),
    x = NULL, y = "NES",
    color = "Diagnosis", fill = "Diagnosis",
    shape = NULL, size = NULL
  ) +
  theme_bw(base_size = 12) +
  theme(
    plot.title      = element_text(face = "bold"),
    legend.position = "bottom"
  )

ggsave(file.path(plot_dir, "Boxplot_category_NES.png"),
       p_box, width = 10, height = 7, dpi = 300)
message("Boxplot guardado.")

# ============================================================
# 5. Plot por categoría individual — detalle de pathways
# ============================================================
# Para categorías con subcategorías (ej. Cellular stress)
# muestra cada pathway individualmente

plot_category_detail <- function(category_name, df_long, plot_dir) {
  
  df_cat <- df_long %>%
    filter(M3_category == category_name)
  
  if (nrow(df_cat) == 0) return(invisible(NULL))
  
  # Ordenar pathways por NES medio
  path_order <- df_cat %>%
    group_by(pathway) %>%
    summarise(mean_NES = mean(NES, na.rm = TRUE)) %>%
    arrange(mean_NES) %>%
    pull(pathway)
  
  # Limpiar nombre del pathway para mostrar
  df_cat <- df_cat %>%
    mutate(
      pathway_clean = gsub("^REACTOME_|^KEGG_", "", pathway),
      pathway_clean = gsub("_", " ", pathway_clean),
      pathway_clean = factor(pathway_clean,
                             levels = gsub("^REACTOME_|^KEGG_", "",
                                           path_order) %>%
                               gsub("_", " ", .))
    )
  
  p <- ggplot(df_cat,
              aes(x = pathway_clean, y = NES,
                  color = diagnosis, shape = significant)) +
    geom_hline(yintercept = 0, color = "grey60",
               linewidth = 0.4, linetype = "dashed") +
    geom_point(aes(size = significant),
               position = position_dodge(width = 0.6),
               alpha = 0.85) +
    geom_linerange(
      aes(ymin = 0, ymax = NES),
      position = position_dodge(width = 0.6),
      linewidth = 0.4, alpha = 0.5
    ) +
    scale_color_manual(values = c(
      TEA_vs_TOC      = "#9b59b6",
      TEA_vs_Tourette      = "#2ecc71",
      TOC_vs_Tourette = "#e67e22"
    )) +
    scale_shape_manual(values = c("TRUE" = 19, "FALSE" = 1),
                       labels = c("TRUE" = "Significant", "FALSE" = "ns")) +
    scale_size_manual(values  = c("TRUE" = 3, "FALSE" = 1.8),
                      labels  = c("TRUE" = "Significant", "FALSE" = "ns")) +
    coord_flip() +
    labs(
      title    = category_name,
      subtitle = paste0("FDR ≤ ", FDR_THR, " = significant"),
      x = NULL, y = "NES",
      color = "Diagnosis", shape = NULL, size = NULL
    ) +
    theme_bw(base_size = 11) +
    theme(
      plot.title      = element_text(face = "bold"),
      legend.position = "bottom",
      axis.text.y     = element_text(size = 8)
    )
  
  filename <- paste0("Detail_",
                     gsub("[^A-Za-z0-9]", "_", category_name), ".png")
  ggsave(file.path(plot_dir, filename),
         p,
         width  = 10,
         height = min(20, max(4, length(unique(df_cat$pathway)) * 0.5 + 2)),
         dpi    = 300)
  message("Detail plot: ", filename)
}

# Generar un plot de detalle por cada categoría
for (cat in unique(df_long$M3_category)) {
  plot_category_detail(cat, df_long, plot_dir)
}

message("\nFin. Plots en: ", plot_dir)

# ============================================================
# 6. Barplot por diagnóstico (uno por figura)
# ============================================================

for (diag in levels(df_long$diagnosis)) {
  
  df_diag <- category_summary %>%
    filter(diagnosis == diag)
  
  p_bar_diag <- df_diag %>%
    mutate(M3_category = factor(M3_category, levels = rev(cat_order))) %>%
    ggplot(aes(x = M3_category, y = mean_NES)) +
    geom_col(fill = "#34495e", alpha = 0.85, width = 0.65) +
    geom_errorbar(
      aes(ymin = mean_NES - sd_NES,
          ymax = mean_NES + sd_NES),
      width = 0.25, linewidth = 0.4, color = "grey30"
    ) +
    geom_hline(yintercept = 0, color = "grey40", linewidth = 0.4) +
    coord_flip() +
    scale_y_continuous(limits = c(-3.5, 3.5)) +
    labs(
      title = paste0("Mean NES per category — ", diag),
      subtitle = paste0("FDR ≤ ", FDR_THR),
      x = NULL, y = "Mean NES"
    ) +
    theme_bw(base_size = 12)
  
  ggsave(file.path(plot_dir, paste0("Barplot_", diag, ".png")),
         p_bar_diag, width = 8, height = 6, dpi = 300)
}

# ============================================================
# 7. Boxplot por diagnóstico (uno por figura)
# ============================================================

for (diag in levels(df_long$diagnosis)) {
  
  df_diag <- df_long %>%
    filter(diagnosis == diag) %>%
    mutate(M3_category = factor(M3_category, levels = rev(cat_order)))
  
  p_box_diag <- ggplot(df_diag, aes(x = M3_category, y = NES)) +
    geom_hline(yintercept = 0, color = "grey60",
               linewidth = 0.4, linetype = "dashed") +
    geom_boxplot(fill = "#3498db", alpha = 0.2, outlier.shape = NA) +
    geom_beeswarm(
      aes(shape = significant, size = significant),
      alpha = 0.8
    ) +
    scale_shape_manual(values = c(19, 1)) +
    scale_size_manual(values = c(2.5, 1.5)) +
    coord_flip() +
    scale_y_continuous(limits = c(-3.5, 3.5)) +
    labs(
      title = paste0("NES distribution — ", diag),
      subtitle = paste0("FDR ≤ ", FDR_THR),
      x = NULL, y = "NES"
    ) +
    theme_bw(base_size = 12)
  
  ggsave(file.path(plot_dir, paste0("Boxplot_", diag, ".png")),
         p_box_diag, width = 8, height = 6, dpi = 300)
}

#############
#Per category_summary
#######################

plot_category_detail <- function(category_name, df_long, plot_dir,
                                 n_top = 7) {
  
  df_cat <- df_long %>%
    filter(M3_category == category_name)
  
  if (nrow(df_cat) == 0) return(invisible(NULL))
  
  # media NES por pathway
  df_summary <- df_cat %>%
    group_by(pathway) %>%
    summarise(mean_NES = mean(NES, na.rm = TRUE),
              .groups = "drop")
  
  # seleccionar TOP UP y DOWN
  top_up <- df_summary %>%
    arrange(desc(mean_NES)) %>%
    head(n_top)
  
  top_down <- df_summary %>%
    arrange(mean_NES) %>%
    head(n_top)
  
  selected_pathways <- unique(c(top_up$pathway, top_down$pathway))
  
  df_cat <- df_cat %>%
    filter(pathway %in% selected_pathways)
  
  # orden
  path_order <- df_summary %>%
    filter(pathway %in% selected_pathways) %>%
    arrange(mean_NES) %>%
    pull(pathway)
  
  df_cat <- df_cat %>%
    mutate(
      pathway_clean = gsub("^REACTOME_|^KEGG_", "", pathway),
      pathway_clean = gsub("_", " ", pathway_clean),
      pathway_clean = stringr::str_wrap(pathway_clean, width = 30),
      pathway_clean = factor(pathway_clean,
                             levels = gsub("^REACTOME_|^KEGG_", "", path_order) %>%
                               gsub("_", " ", .) %>%
                             stringr::str_wrap(width = 30))
    )
  
  p <- ggplot(df_cat,
              aes(x = pathway_clean, y = NES,
                  color = diagnosis, shape = significant)) +
    geom_hline(yintercept = 0, color = "grey60",
               linewidth = 0.4, linetype = "dashed") +
    geom_point(aes(size = significant),
               position = position_dodge(width = 0.6),
               alpha = 0.85) +
    geom_linerange(
      aes(ymin = 0, ymax = NES),
      position = position_dodge(width = 0.6),
      linewidth = 0.4, alpha = 0.5
    ) +
    coord_flip() +
    labs(
      title = paste0(category_name, " (Top ", n_top, " up/down)"),
      x = NULL, y = "NES"
    ) +
    theme_bw(base_size = 11) +
    theme(
      axis.text.y = element_text(size = 8),
      plot.margin = margin(10, 10, 10, 100)
    )
  
  filename <- paste0("Detail_TOP_",
                     gsub("[^A-Za-z0-9]", "_", category_name), ".png")
  
  ggsave(file.path(plot_dir, filename),
         p, width = 9, height = 10, dpi = 300)
}

for (cat in unique(df_long$M3_category)) {
  plot_category_detail(cat, df_long, plot_dir, n_top = 10)
}



################################################################################

# ============================================================
# 8. Signature matrix (categoría × diagnóstico)
# ============================================================

# Puedes cambiar esto:
USE_ONLY_SIG <- TRUE

df_sig <- if (USE_ONLY_SIG) {
  df_long %>% filter(significant)
} else {
  df_long
}

signature_mat <- df_sig %>%
  group_by(M3_category, diagnosis) %>%
  summarise(
    mean_NES = mean(NES, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  tidyr::pivot_wider(
    names_from = diagnosis,
    values_from = mean_NES
  )

# Guardar tabla
write.csv(signature_mat,
          file.path(plot_dir, "Signature_matrix.csv"),
          row.names = FALSE)

library(pheatmap)

mat <- signature_mat %>%
  as.data.frame()

rownames(mat) <- mat$M3_category
mat$M3_category <- NULL

mat <- as.matrix(mat)

mat <- mat[complete.cases(mat), ]
mat <- mat[apply(mat, 1, sd, na.rm = TRUE) > 0, ]

p_heat <- pheatmap(
  mat,
  scale = "row",
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  fontsize_row = 10,
  main = "Category signature (NES)"
)

png(file.path(plot_dir, "Signature_heatmap.png"),
    width = 1800, height = 1800, res = 300)

grid::grid.newpage()
grid::grid.draw(p_heat$gtable)

dev.off()

df_sig_plot <- signature_mat %>%
  pivot_longer(
    cols = -M3_category,
    names_to = "diagnosis",
    values_to = "mean_NES"
  )

p_sig <- ggplot(df_sig_plot,
                aes(x = diagnosis,
                    y = reorder(M3_category, mean_NES),
                    color = mean_NES,
                    size = abs(mean_NES))) +
  geom_point() +
  scale_color_gradient2(
    low = "#3498db",
    mid = "white",
    high = "#e74c3c",
    midpoint = 0
  ) +
  labs(
    title = "Category signature per diagnosis",
    x = NULL,
    y = NULL,
    color = "Mean NES",
    size = "|NES|"
  ) +
  theme_bw(base_size = 12)

ggsave(file.path(plot_dir, "Signature_dotplot.png"),
       p_sig, width = 7, height = 6, dpi = 300)

signature_weighted <- df_sig %>%
  group_by(M3_category, diagnosis) %>%
  summarise(
    mean_NES = mean(NES),
    n = n(),
    weighted = mean_NES * log1p(n),
    .groups = "drop"
  )









