library(tidyverse)
library(Seurat)
library(tidyseurat)
library(cowplot)
library(patchwork)
theme_set(theme_bw())

# Load published data and create Seurat object
UMI_counts <- read.delim(file = "GSE167650_Whole_Processed_umi_counts.tsv")
meta.data <- read.delim(file = "GSE167650_Whole_Processed_annotations.tsv")
UMI_counts <- as(UMI_counts, "sparseMatrix")
Kang <- CreateSeuratObject(UMI_counts, meta.data = meta.data)

# Standard Seurat workflow
Kang <- NormalizeData(Kang)
Kang <- FindVariableFeatures(Kang)
Kang <- ScaleData(Kang)
Kang <- RunPCA(Kang)
Kang <- FindNeighbors(Kang, dims = 1:30, reduction = "pca")
Kang <- FindClusters(Kang, resolution = 1)
Kang <- RunUMAP(Kang, dims = 1:30, reduction = "pca")
DimPlot(Kang, reduction = "umap", group.by = "orig.ident")

# Rename and Reorder identities alphabetically
Idents(Kang) <- "Celltype"
DimPlot(Kang, reduction = "umap")
Kang <- RenameIdents(Kang, 'B_Macrophages'="Mixed",
                         'Neutrophil_Macrophages'="Mixed",
                         'Monocytes_Ly6c1+'="Monocytes",
                         'B_Neutrophil'="Mixed",
                         'Bcell'="B cells",
                         'T_CD8'="T cells",
                         'T_CD4'="T cells")
Kang$celltype_new <- Idents(Kang)
current_levels <- levels(Idents(Kang))
new_levels <- sort(current_levels)
Idents(Kang) <- factor(Idents(Kang), levels = new_levels)

# Define consistent cell type colors (T cells vs. myeloid clearly different)
standard_colors <- c(
  # T cells
  "T cells" = "#d73027",
  "NK" = "#e34a33",
  
  # Myeloid
  "Macrophages" = "#4575b4",
  "Monocytes" = "#91bfdb",
  "DC" = "#74a9cf",
  "Neutrophils" = "#2c7bb6",
  
  
  # Other
  "B cells" = "#999999",
  "Mixed" = "#f0f0f0",
  "ILC" = "#dddddd"
  
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
plot_umap_standardized(Kang, title_text = "Mouse (Kang et al. 2021)")

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
plot_dot_standardized(Kang, title_text = "Mouse")

# Define myeloid cell types
myeloid_cell_types <- c(
  "Neutrophils", "Monocytes", "Macrophages", "DC"
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
  seurat_obj = Kang,
  gene = "Nos2",
  group_col = "Group",
  group_levels = c("Naive", "MtbK"),
  celltype_col = "celltype_new"
)

# Export session info to a text file
sink("Kang_session_info.txt")
sessionInfo()
sink()