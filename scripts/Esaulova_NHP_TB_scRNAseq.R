library(tidyverse)
library(Seurat)
library(tidyseurat)
library(cowplot)
library(patchwork)
theme_set(theme_bw())

# Load published data and generate Seurat object
load(file="./GSE149758_CD45_expData.Rda")
meta.data <- read_delim(file = "./GSE149758_CD45_meta_data.txt")
Esaulova <- CreateSeuratObject(expData, meta.data = meta.data)

# Standard Seurat workflow
Esaulova <- NormalizeData(object = Esaulova)
Esaulova <- FindVariableFeatures(object = Esaulova)
Esaulova <- ScaleData(object = Esaulova)
Esaulova <- RunPCA(object = Esaulova)
Esaulova <- FindNeighbors(object = Esaulova, dims = 1:20)
Esaulova <- FindClusters(object = Esaulova, resolution = 0.25)
Esaulova <- RunUMAP(object = Esaulova, dims = 1:20)
DimPlot(Esaulova, reduction = "umap")

# Keep only clusters of interest and rename clusters
Idents(Esaulova) <- "Cluster_new"
Esaulova <- subset(Esaulova, idents = c("B cells",
                                        "Myeloid 6",
                                        "Mast cells",
                                        "T cells 2",
                                        "proliferating",
                                        "Non-immune 4",
                                        "Myeloid 8",
                                        "pDC",
                                        "T cells 1",
                                        "Plasmablasts",
                                        "Non-immune 3",
                                        "Non-immune 5",
                                        "Non-immune 1",
                                        "Non-immune 2"))

Esaulova <- RenameIdents(Esaulova, 'Myeloid 6'="Myeloid",
                         'T cells 2'="T cells",
                         'Non-immune 4'="Non-immune",
                         'Myeloid 8'="Myeloid",
                         'T cells 1'="T cells",
                         'Non-immune 3'="Non-immune",
                         'Non-immune 5'="Non-immune",
                         'Non-immune 1'="Non-immune",
                         'Non-immune 2'="Non-immune")

Esaulova$cell_type <- Idents(Esaulova)
Idents(Esaulova) <- "cell_type"

# Reorder identities alphabetically
current_levels <- levels(Idents(Esaulova))
new_levels <- sort(current_levels)
Idents(Esaulova) <- factor(Idents(Esaulova), levels = new_levels)

# Define consistent cell type colors (T cells vs. myeloid clearly different)
standard_colors <- c(
  # T cells
  "T cells" = "#d73027",
  
  # Myeloid
  "Myeloid" = "#4575b4",
  "Mast cells" = "#74a9cf",
  "pDC" = "#2c7bb6",
  
  # Other
  "B cells" = "#999999",
  "Non-immune" = "#f0f0f0",
  "proliferating" = "#cccccc",
  "Plasmablasts" = "#bbbbbb"
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
plot_umap_standardized(Esaulova, title_text = "Non-human primates (Esaulova et al. 2020)")

# Function to plot a standardized dot plot 
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
plot_dot_standardized(Esaulova, title_text = "Non-human primates")

# Define myeloid cell types
myeloid_cell_types <- c(
  "Myeloid", "pDC", "Mast cells"
)

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
      fill = c("white", "darkorange"),  # Fill colors for each group
      color = c("darkorange", "black")  # Border colors
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
  seurat_obj = Esaulova,
  gene = "NOS2",
  group_col = "Condition",
  group_levels = c("control", "active"),
  celltype_col = "cell_type"
)

# Export session info to a text file
sink("Esaulova_session_info.txt")
sessionInfo()
sink()