options(future.globals.maxSize= 891289600)
library(tidyverse)
library(Seurat)
library(patchwork)
library(cowplot)
library(Azimuth)

# Load published data and create Seurat object
Hillman.data <- Read10X("L:/RD-Corleis/External_scRNAseq/8_GSE214237_TB_human_PBMCs_Hillman/8_download")
meta.data <- read.delim("./GSE214237_HaHi02_metadata.tsv")
rownames(meta.data) <- meta.data$Cell_Barcode
Hillman <- CreateSeuratObject(Hillman.data$`Gene Expression`, meta.data = meta.data)

# Standard Seurat workflow
Hillman <- NormalizeData(Hillman)
Hillman <- FindVariableFeatures(Hillman)
Hillman <- ScaleData(Hillman)
Hillman <- RunPCA(Hillman)
Hillman <- FindNeighbors(Hillman, dims = 1:20)
Hillman <- FindClusters(Hillman)
Hillman <- RunUMAP(Hillman, dims = 1:20)
DimPlot(Hillman, reduction = "umap")

# Identify cell types with Azimuth
Hillman <- RunAzimuth(Hillman, reference = "pbmcref")
Idents(Hillman) <- "predicted.celltype.l2"
DimPlot(Hillman, reduction = "umap")

# Reorder identities alphabetically
current_levels <- levels(Idents(Hillman))
new_levels <- sort(current_levels)
Idents(Hillman) <- factor(Idents(Hillman), levels = new_levels)

# Define consistent cell type colors (T cells vs. myeloid clearly different)
standard_colors <- c(
  # T cells
  "CD4 TCM" = "#d73027",
  "CD8 TCM" = "#fc8d59",
  "CD8 TEM" = "#FF3F34",
  "Treg" = "#b2182b",
  "gdT" = "#f46d43",
  "CD4 CTL" = "#ef6548",
  "NK" = "#e34a33",
  "NK Proliferating" = "#D1001C",
  
  # Myeloid
  "pDC" = "#4575b4",
  "CD14 Mono" = "#91bfdb",
  "CD16 Mono" = "#74add1",
  "cDC2" = "#2c7bb6",
  
  # Other
  "B intermediate" = "#999999",
  "B memory" = "#bbbbbb",
  "B naive" = "#D9D9D9",
  "Eryth" = "#f0f0f0",
  "HSPC" = "#dddddd",
  "Plasmablast" = "#dddddd",
  "Platelet" = "#cccccc"
)

# Function to plot a standardized UMAP
plot_umap_standardized <- function(seurat_obj, title_text = "UMAP", label_clusters = FALSE) {
  DimPlot(seurat_obj, reduction = "umap", label = label_clusters, raster = FALSE) +
    scale_color_manual(values = standard_colors, na.translate = FALSE) +
    ggtitle(title_text) +
    theme_void() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 24),
      plot.title.position = "plot",
      legend.text = element_text(size = 20),
      legend.title =  element_text(size = 20, face = "bold"),
      plot.margin = margin(t = 10, r = 10, b = 10, l = 10, unit = "pt")
    )
}

# Generate the UMAP
plot_umap_standardized(Hillman, title_text = "Human (Hillman et al. 2023)")

# Function for standardized dot plot
plot_dot_standardized <- function(seurat_obj, features = c("CXCL10", "IFNG", "NOS2"), 
                                  title_text = NULL, scale_colors = c("grey90", "purple", "blue")) {
  
  DotPlot(seurat_obj, features = features) +
    scale_color_gradientn(colors = scale_colors, name = "Average Expression") +
    theme_minimal(base_size = 14) +
    labs(x = NULL, y = NULL, title = title_text) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 24),
      axis.text.x = element_text(angle = 0, hjust = 0.5, face = "bold", size = 20),
      axis.text.y = element_text(face = "bold", size = 18),
      legend.text = element_text(size = 16),
      legend.title = element_text(size = 16, face = "bold"),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    ) +
    guides(
      color = guide_colorbar(order = 1),
      size = guide_legend(order = 2, title = "Percent Expressed")
    )
}

# Generate dot plot
plot_dot_standardized(Hillman, title_text = "Human")

# Export session info to a text file
sink("Hillman_session_info.txt")
sessionInfo()
sink()