library(Seurat)
library(ggplot2)
library(dplyr)
library(patchwork)

# Set relative paths
expression_file <- "data/expression_matrix.txt"
metadata_file   <- "data/metadata.tsv"
umap_file       <- "data/umap_layout.txt"

# Load raw count matrix
raw_counts <- read.delim(expression_file, row.names = 1, check.names = FALSE)
seurat_obj <- CreateSeuratObject(counts = raw_counts, project = "SCP1289")

# Load and integrate metadata
metadata <- read.delim(metadata_file, row.names = 1, check.names = FALSE)
metadata <- metadata[colnames(seurat_obj), ]
seurat_obj <- AddMetaData(seurat_obj, metadata)

# Load UMAP coordinates and cell types
umap_data <- read.delim(umap_file, row.names = 1, skip = 1)
colnames(umap_data) <- c("UMAP_1", "UMAP_2", "cell_type")
umap_embeddings <- apply(as.matrix(umap_data[, c("UMAP_1", "UMAP_2")]), 2, as.numeric)
rownames(umap_embeddings) <- rownames(umap_data)
seurat_obj[["umap"]] <- CreateDimReducObject(embeddings = umap_embeddings, key = "UMAP_", assay = "RNA")
seurat_obj$cell_type <- umap_data$cell_type

# Normalize expression data
seurat_obj <- NormalizeData(seurat_obj, assay = "RNA")

# Standardized UMAP Plot
standard_colors <- c(
  "T Cells" = "#d73027", "B Cells" = "#fc8d59", "Plasma Cells" = "#fdae61",
  "Macrophages" = "#4575b4", "Monocytes" = "#74add1", "Dendritic Cells" = "#2c7bb6",
  "Mast Cells" = "#91bfdb", "Neutrophils" = "#313695",
  "Ciliated Cells" = "#bdbdbd", "Basal Cells" = "#969696", "Squamous Cells" = "#cccccc",
  "Club Cells" = "#e0e0e0", "Ionocytes" = "#f0f0f0", "Goblet Cells" = "#1a1a1a",
  "Erythrocytes" = "#f7f7f7", "Unknown / low signal" = "#bbbbbb"
)

plot_umap_standardized <- function(seurat_obj, title_text = "Human (Ziegler et al. 2021)") {
  DimPlot(seurat_obj, reduction = "umap", group.by = "cell_type", label = FALSE, raster = FALSE) +
    scale_color_manual(values = standard_colors, na.translate = FALSE) +
    ggtitle(title_text) +
    theme_void() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 24),
      legend.position = "right",
      legend.text = element_text(size = 20),
      legend.title = element_blank(),
      plot.margin = margin(10, 10, 10, 10, "pt")
    )
}

ggsave("plots/UMAP_overview.png", plot_umap_standardized(seurat_obj), width = 9, height = 7)

# DotPlot
plot_dot_standardized <- function(seurat_obj, features = c("CXCL10", "IFNG", "NOS2"),
                                  title_text = "Human BAL (COVID-19)",
                                  scale_colors = c("grey90", "purple", "blue")) {
  DotPlot(seurat_obj, features = features, group.by = "cell_type") +
    scale_color_gradientn(colors = scale_colors, name = "Avg Expression") +
    theme_minimal(base_size = 14) +
    labs(x = NULL, y = NULL, title = title_text) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 24),
      axis.text.x = element_text(angle = 45, hjust = 1, face = "bold", size = 14),
      axis.text.y = element_text(face = "bold", size = 14),
      legend.text = element_text(size = 14),
      legend.title = element_text(size = 14, face = "bold")
    )
}

ggsave("plots/DotPlot_CXCL10_IFNG_NOS2.png", plot_dot_standardized(seurat_obj), width = 10, height = 6)


# Boxplots (Macrophages & Secretory Cells)

# Define helper
plot_expression_boxplot <- function(seurat_obj, gene, celltype_filter, condition_col = "disease__ontology_label",
                                    conditions = c("normal", "COVID-19"), title_text) {
  expr_data <- GetAssayData(seurat_obj, assay = "RNA", slot = "data")
  df <- data.frame(
    Cell = colnames(seurat_obj),
    Expression = expr_data[gene, ],
    Condition = seurat_obj[[condition_col]][,1],
    Celltype = seurat_obj$cell_type
  )
  
  df <- df %>%
    filter(Celltype %in% celltype_filter, Condition %in% conditions) %>%
    mutate(Condition = factor(Condition, levels = conditions))
  
  ggplot(df, aes(x = Condition, y = Expression)) +
    geom_boxplot(width = 0.6, fill = c("white", "firebrick"), color = c("black", "firebrick")) +
    theme_minimal(base_size = 14) +
    theme(
      plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
      axis.title.y = element_text(face = "bold"),
      panel.border = element_rect(color = "black", fill = NA),
      axis.text.x = element_text(face = "bold", size = 14)
    ) +
    labs(y = "Expression", title = title_text)
}

# Create plots
macrophage_box <- plot_expression_boxplot(seurat_obj, gene = "NOS2", celltype_filter = "Macrophages",
                                          title_text = "NOS2 (Macrophages)")
secretory_types <- c("Secretory Cells", "Goblet Cells", "Club Cells")
secretory_box <- plot_expression_boxplot(seurat_obj, gene = "NOS2", celltype_filter = secretory_types,
                                         title_text = "NOS2 (Secretory Cells)")

# Save combined figure
combined_box <- macrophage_box + secretory_box + plot_layout(ncol = 2)
ggsave("plots/NOS2_boxplot_macrophages_secretory.png", combined_box, width = 10, height = 5)

# Save Session
save(seurat_obj, file = "output/seurat_obj_SCP1289.RData")

# Export session info to a text file
sink("Ziegler_human_session_info.txt")
sessionInfo()
sink()