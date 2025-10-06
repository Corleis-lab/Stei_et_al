library(tidyverse)
library(Seurat)
library(tidyseurat)
library(cowplot)
library(patchwork)
library(sctransform)
theme_set(theme_bw())

# Load published data
Pisu <- readRDS("./GSE167232_mtb_integrated.RDS")
Idents(Pisu) <- "seurat_clusters"
DimPlot(Pisu, reduction = "umap", label = TRUE)

# Rename cluster as shown in the paper
Pisu <- RenameIdents(Pisu, '0'="AM_1",
                     '1'="IM_3",
                     '2'="AM_2",
                     '3'="IM_2",
                     '4'="IM_1",
                     '5'="AM_3",
                     '6'="AM_4",
                     '7'="IM_4",
                     '8'="Monocytes",
                     '9'="Neutrophils",
                     '10'="DC",
                     '11'="Macrophages/T cells",
                     '12'="DC.103+11B-",
                     '13'="Macrophages/B cells")
Pisu$cell_types <- Idents(Pisu)
Idents(Pisu) <- "cell_types"

# Reorder identities alphabetically
current_levels <- levels(Idents(Pisu))
new_levels <- sort(current_levels)
Idents(Pisu) <- factor(Idents(Pisu), levels = new_levels)

# Define consistent cell type colors (T cells vs. myeloid clearly different)
standard_colors <- c(
  "Macrophages/T cells" = "#d73027",
  
  # Myeloid
  "AM_1" = "#4575b4",
  "Monocytes" = "#91bfdb",
  "DC" = "#A569BD",
  "DC.103+11B-" = "#313695",
  "Neutrophils" = "#2c7bb6",
  "AM_2" = "#74a9cf",
  "AM_3" = "#1F77B4",
  "AM_4" = "#2E86C1",
  "IM_1" = "#5DADE2",
  "IM_2" = "#85C1E9",
  "IM_3" = "#2980B9",
  "IM_4" = "#74add1",
  
  # Other
  "Macrophages/B cells" = "#999999"
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
plot_umap_standardized(Pisu, title_text = "Mouse (Pisu et al. 2021)")

# Function for standardized dot plot
plot_dot_standardized <- function(seurat_obj, features = c("Cxcl10", "Ifng", "Nos2"), 
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
plot_dot_standardized(Pisu, title_text = "Mouse")

# Define myeloid cell types
myeloid_cell_types <- levels(Pisu)

# Universal boxplot function
plot_myeloid_bin_boxplot <- function(seurat_obj,
                                     gene = "NOS2",
                                     group_col = "HL_group",
                                     group_levels = c("Low", "High"),
                                     celltype_col = "cell_type_collapsed",
                                     myeloid_types = myeloid_cell_types) {
  
  # Check required metadata
  missing_cols <- setdiff(c(group_col, celltype_col), colnames(seurat_obj@meta.data))
  if (length(missing_cols) > 0) {
    stop(paste("Missing metadata column(s):", paste(missing_cols, collapse = ", ")))
  }
  
  # Get expression data (Seurat v5 compatible)
  expr_data <- GetAssayData(seurat_obj, layer = "data")
  if (!(gene %in% rownames(expr_data))) {
    stop(paste("Gene", gene, "not found in expression matrix."))
  }
  
  # Build base data frame
  df <- data.frame(
    Cell = colnames(seurat_obj),
    Group = seurat_obj@meta.data[[group_col]],
    Celltype = seurat_obj@meta.data[[celltype_col]],
    Expression = expr_data[gene, colnames(seurat_obj)],
    stringsAsFactors = FALSE
  )
  
  # Filter to desired cells
  df_filtered <- df %>%
    filter(Group %in% group_levels, Celltype %in% myeloid_types) %>%
    mutate(
      Group = factor(Group, levels = group_levels),
      Gene = gene
    )
  
  if (nrow(df_filtered) == 0) {
    stop("No matching cells found after filtering. Check group names and cell types.")
  }
  
  # Create boxplot
  p <- ggplot(df_filtered, aes(x = Group, y = Expression)) +
    geom_boxplot(
      width = 0.6,                     # Medium box width (not too narrow)
      fill = c("white", "lightskyblue"),  # Fill colors for each group
      color = c("lightskyblue", "black")  # Border colors
    ) +
    theme_minimal(base_size = 16) +
    theme(
      plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
      axis.title.x = element_blank(),
      axis.title.y = element_text(face = "bold"),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
      plot.margin = margin(5, 5, 5, 5),
      axis.text.x = element_text(size = 14, face = "bold")
    ) +
    labs(y = "Expression", title = gene)
  
  return(p)
}


# Generate box plot
plot_myeloid_bin_boxplot(
  seurat_obj = Pisu,
  gene = "Nos2",
  group_col = "Treatment",
  group_levels = c("Uninfected", "Infected"),
  celltype_col = "cell_types"
)

# Export session info to a text file
sink("Pisu_session_info.txt")
sessionInfo()
sink()