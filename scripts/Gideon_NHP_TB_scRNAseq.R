library(tidyverse)
library(Seurat)
library(tidyseurat)
library(cowplot)
library(patchwork)
library(sctransform)
theme_set(theme_bw())

# Set up Seurat object and add published metadata
Gideon.10wk.data <- Read10X(data.dir = "./data/10_weeks", gene.column = 1)
Gideon_10wk <- CreateSeuratObject(counts = Gideon.10wk.data, 
                                  project = "Gideon_10wk", 
                                  min.cells = 3, 
                                  min.features = 200)
wk10_metadata <- read.delim("./data/10_weeks/GSE200151_Updated10wk_alexandria_structured_metadata10.txt")
rownames(wk10_metadata) <- wk10_metadata$NAME
Gideon_10wk <- AddMetaData(Gideon_10wk, wk10_metadata)

# Standard Seurat workflow
Gideon_10wk <- NormalizeData(object = Gideon_10wk)
Gideon_10wk <- FindVariableFeatures(object = Gideon_10wk)
Gideon_10wk <- ScaleData(object = Gideon_10wk)
Gideon_10wk <- RunPCA(object = Gideon_10wk)
Gideon_10wk <- FindNeighbors(object = Gideon_10wk, dims = 1:30)
Gideon_10wk <- FindClusters(object = Gideon_10wk)
Gideon_10wk <- RunUMAP(object = Gideon_10wk, dims = 1:30)
Idents(Gideon_10wk) <- "GenericFinal"
DimPlot(object = Gideon_10wk, reduction = "umap")

# Set up Seurat object and add published metadata
Gideon.4wk.data <- Read10X(data.dir = "./data/4_weeks", gene.column = 1)
Gideon_4wk <- CreateSeuratObject(counts = Gideon.4wk.data, 
                                 project = "Gideon_4wk", 
                                 min.cells = 3, 
                                 min.features = 200)
wk4_metadata <- read.delim("./data/4_weeks/GSE200151_Updated4wk_alexandria_structured_metadata3.txt")
rownames(wk4_metadata) <- wk4_metadata$NAME
Gideon_4wk <- AddMetaData(Gideon_4wk, wk4_metadata)

# Standard Seurat workflow
Gideon_4wk <- NormalizeData(object = Gideon_4wk)
Gideon_4wk <- FindVariableFeatures(object = Gideon_4wk)
Gideon_4wk <- ScaleData(object = Gideon_4wk)
Gideon_4wk <- RunPCA(object = Gideon_4wk)
Gideon_4wk <- FindNeighbors(object = Gideon_4wk, dims = 1:30)
Gideon_4wk <- FindClusters(object = Gideon_4wk)
Gideon_4wk <- RunUMAP(object = Gideon_4wk, dims = 1:30)
Idents(Gideon_4wk) <- "CellTypeAnnotations"
DimPlot(object = Gideon_4wk, reduction = "umap")

# Merge Seurat objects
Gideon <- merge(Gideon_10wk, Gideon_4wk, add.cell.ids = c("wk10", "wk4"), merge.dr = TRUE)
Gideon <- Gideon %>% mutate(celltype_new = coalesce(CellTypeAnnotations, GenericFinal))
Gideon$week <- sub("\\_.*", "", rownames(Gideon@meta.data))
Gideon <- NormalizeData(Gideon)
Gideon <- FindVariableFeatures(Gideon)
Gideon <- ScaleData(Gideon)
Gideon <- RunPCA(Gideon)
Gideon <- FindNeighbors(Gideon, dims = 1:30, reduction = "pca")
Gideon <- FindClusters(Gideon, cluster.name = "unintegrated_clusters")
Gideon <- RunUMAP(Gideon, dims = 1:30, reduction = "pca", reduction.name = "umap.unintegrated")
Gideon <- IntegrateLayers(object = Gideon, method = CCAIntegration, orig.reduction = "pca", new.reduction = "integrated.cca", verbose = FALSE)
Gideon[["RNA"]] <- JoinLayers(Gideon[["RNA"]])
Gideon <- FindNeighbors(Gideon, reduction = "integrated.cca", dims = 1:30)
Gideon <- FindClusters(Gideon)
Gideon <- RunUMAP(Gideon, dims = 1:30, reduction = "integrated.cca")
DimPlot(Gideon, reduction = "umap", group.by = "week")

# Rename and Reorder identities alphabetically
Idents(Gideon) <- "celltype_new"
Gideon <- RenameIdents(Gideon, 'B' = "B cells",
                       'Endo' = "Endothelial cells",
                       'Endothelial' = "Endothelial cells",
                       'Fibro' = "Fibroblasts",
                       'Fibroblast' = "Fibroblasts",
                       'Mast' = "Mast cells",
                       'Macrophage' = "Macrophages",
                       'Mphage' = "Macrophages",
                       'Neutrophil' = "Neutrophils",
                       'Plasma' = "Plasma cells",
                       'T' = "T cells",
                       'T1P' = "Type 1 pneumocytes",
                       'T2P' = "Type 2 pneumocytes")
Gideon$celltype_new <- Idents(Gideon)
current_levels <- levels(Idents(Gideon))
new_levels <- sort(current_levels)
Idents(Gideon) <- factor(Idents(Gideon), levels = new_levels)

# Define consistent cell type colors (T cells vs. myeloid clearly different)
standard_colors <- c(
  # T cells
  "T cells" = "#d73027",
  
  # Myeloid
  "Macrophages" = "#4575b4",
  "cDC" = "#74add1",
  "pDC" = "#2c7bb6",
  "Neutrophils" = "#313695",
  "Mast cells" = "#74a9cf",
  
  # Other
  "Type 1 pneumocytes" = "#D9D9D9",
  "Type 2 pneumocytes" = "#BFBFBF",
  "B cells" = "#999999",
  "RBC" = "#f0f0f0",
  "Fibroblasts" = "#dddddd",
  "Endothelial cells" = "#cccccc",
  "Plasma cells" = "#bbbbbb"
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
plot_umap_standardized(Gideon, title_text = "Non-human primates (Gideon et al. 2022)")


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
plot_dot_standardized(Gideon, title_text = "Non-human primates")

# Define myeloid cell types
myeloid_cell_types <- c(
  "cDC", "pCD", "Mast cells", "Neutrophils", "Macrophages"
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
  seurat_obj = Gideon,
  gene = "NOS2",
  group_col = "week",
  group_levels = c("wk4", "wk10"),
  celltype_col = "celltype_new"
)

# Export session info to a text file
sink("Gideon_session_info.txt")
sessionInfo()
sink()