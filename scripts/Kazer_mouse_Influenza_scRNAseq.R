library(Seurat)
library(ggplot2)
library(dplyr)
library(patchwork)

# Load data
load("Kazer.RData")  # Contains object: epi_cells

# Generate UMAP
p_umap <- DimPlot(
  epi_cells,
  reduction = "umap",
  group.by = "cell.type2",   # or use "cluster_label"
  label = FALSE
) +
  ggtitle("Mouse (Kazer et al. 2024)") +
  theme_void() +
  theme(
    plot.title = element_text(size = 24, hjust = 0.5, face = "bold"),
    legend.text = element_text(size = 18),
    legend.title = element_blank(),
    legend.key.size = unit(1, "cm")
  )

ggsave("figures/umap_celltypes.png", p_umap, width = 8, height = 6, dpi = 300)


# FeaturePlot for Nos2 and Scgb1c1
feature_colors <- c("lightgrey", "red")

p_nos2 <- FeaturePlot(
  epi_cells,
  features = "Nos2",
  reduction = "umap",
  cols = feature_colors
) +
  ggtitle("Nos2") +
  theme_void() +
  theme(plot.title = element_text(hjust = 0.5, size = 20))

p_scgb1c1 <- FeaturePlot(
  epi_cells,
  features = "Scgb1c1",
  reduction = "umap",
  cols = feature_colors
) +
  ggtitle("Scgb1c1") +
  theme_void() +
  theme(plot.title = element_text(hjust = 0.5, size = 20))

p_combined <- p_nos2 / p_scgb1c1
ggsave("figures/featureplot_nos2_scgb1c1.png", p_combined, width = 8, height = 10, dpi = 300)


# Boxplot for Nos2 in Krt13+Il1a+ cells
# Create infection status label
epi_cells$InfectionStatus <- ifelse(epi_cells$timepoint == "Naive", "Naive", "Infected")

# Boxplot function
plot_krt13_il1a_boxplot <- function(seurat_obj,
                                    gene = "Nos2",
                                    cluster_label = "Krt13+Il1a+ Epi",
                                    group_col = "InfectionStatus",
                                    group_levels = c("Naive", "Infected"),
                                    cluster_col = "cluster_label") {
  
  expr_data <- GetAssayData(seurat_obj, layer = "data")
  df <- data.frame(
    Cell = colnames(seurat_obj),
    Group = seurat_obj@meta.data[[group_col]],
    Cluster = seurat_obj@meta.data[[cluster_col]],
    Expression = expr_data[gene, colnames(seurat_obj)],
    stringsAsFactors = FALSE
  )
  
  df_filtered <- df %>%
    filter(Cluster == cluster_label, Group %in% group_levels) %>%
    mutate(Group = factor(Group, levels = group_levels), Gene = gene)
  
  p <- ggplot(df_filtered, aes(x = Group, y = Expression)) +
    geom_boxplot(
      width = 0.6,
      fill = c("white", "firebrick"),
      color = c("firebrick", "black")
    ) +
    theme_minimal(base_size = 16) +
    theme(
      plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
      axis.title.x = element_blank(),
      axis.title.y = element_text(face = "bold"),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
      axis.text.x = element_text(size = 14, face = "bold")
    ) +
    labs(y = "Normalized Expression", title = gene)
  
  return(p)
}

p_box <- plot_krt13_il1a_boxplot(epi_cells)
ggsave("figures/boxplot_Nos2_Krt13Il1aEpi.png", p_box, width = 5, height = 6, dpi = 300)

# Export session info to a text file
sink("Kazer_session_info.txt")
sessionInfo()
sink()