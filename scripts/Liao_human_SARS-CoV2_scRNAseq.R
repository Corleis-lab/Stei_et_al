options(future.globals.maxSize= 891289600)
library(tidyverse)
library(Seurat)
library(tidyseurat)
library(cowplot)
library(patchwork)
library(Azimuth)
theme_set(theme_bw())

# Load published data
C141.data <- Read10X_h5("./H5/GSM4339769_C141_filtered_feature_bc_matrix.h5")
C142.data <- Read10X_h5("./H5/GSM4339770_C142_filtered_feature_bc_matrix.h5")
C143.data <- Read10X_h5("./H5/GSM4339771_C143_filtered_feature_bc_matrix.h5")
C144.data <- Read10X_h5("./H5/GSM4339772_C144_filtered_feature_bc_matrix.h5")
C145.data <- Read10X_h5("./H5/GSM4339773_C145_filtered_feature_bc_matrix.h5")
C146.data <- Read10X_h5("./H5/GSM4339774_C146_filtered_feature_bc_matrix.h5")
C51.data <- Read10X_h5("./H5/GSM4475048_C51_filtered_feature_bc_matrix.h5")
C52.data <- Read10X_h5("./H5/GSM4475049_C52_filtered_feature_bc_matrix.h5")
C100.data <- Read10X_h5("./H5/GSM4475050_C100_filtered_feature_bc_matrix.h5")
C148.data <- Read10X_h5("./H5/GSM4475051_C148_filtered_feature_bc_matrix.h5")
C149.data <- Read10X_h5("./H5/GSM4475052_C149_filtered_feature_bc_matrix.h5")
C152.data <- Read10X_h5("./H5/GSM4475053_C152_filtered_feature_bc_matrix.h5")

# Create Seurat objects
C141 <- CreateSeuratObject(C141.data, project = "mild")
C142 <- CreateSeuratObject(C142.data, project = "mild")
C143 <- CreateSeuratObject(C143.data, project = "severe")
C144 <- CreateSeuratObject(C144.data, project = "mild")
C145 <- CreateSeuratObject(C145.data, project = "severe")
C146 <- CreateSeuratObject(C146.data, project = "severe")
C51 <- CreateSeuratObject(C51.data, project = "CTR")
C52 <- CreateSeuratObject(C52.data, project = "CTR")
C100 <- CreateSeuratObject(C100.data, project = "CTR")
C148 <- CreateSeuratObject(C148.data, project = "severe")
C149 <- CreateSeuratObject(C149.data, project = "severe")
C152 <- CreateSeuratObject(C152.data, project = "severe")

# Merge Seurat objects
mild <- merge(C141, y = c(C142, C144), add.cell.ids = c("C141", "C142", "C144"), project = "mild")
severe <- merge(C143, y = c(C145, C146, C148, C149, C152), add.cell.ids = c("C143", "C145", "C146", "C148", "C149", "C152"), project = "severe")
CTR <- merge(C51, y = c(C52, C100), add.cell.ids = c("C51", "C52", "C100"), project = "CTR")
Liao <- merge(CTR, y = c(mild, severe), add.cell.ids = c("CTR", "mild", "severe"), project = "Liao")
Liao$barcode <- rownames(Liao@meta.data)
Liao$Donor <- str_match(Liao$barcode, "_\\s*(.*?)\\s*_")[,2]
Liao <- JoinLayers(Liao)
Liao[["RNA"]] <- split(Liao[["RNA"]], f = Liao$orig.ident)
Liao <- NormalizeData(Liao)
Liao <- FindVariableFeatures(Liao)
Liao <- ScaleData(Liao)
Liao <- RunPCA(Liao)
Liao <- FindNeighbors(Liao, dims = 1:30, reduction = "pca")
Liao <- FindClusters(Liao, resolution = 1.2, cluster.name = "unintegrated_clusters")
Liao <- RunUMAP(Liao, dims = 1:30, reduction = "pca", reduction.name = "umap.unintegrated")
DimPlot(Liao, reduction = "umap.unintegrated", group.by = "orig.ident")
Liao <- IntegrateLayers(
  object = Liao, method = CCAIntegration,
  orig.reduction = "pca", new.reduction = "integrated.cca",
  verbose = FALSE
)
Liao <- FindNeighbors(Liao, reduction = "integrated.cca", dims = 1:30)
Liao <- FindClusters(Liao, resolution = 1.2, cluster.name = "cca_clusters")
Liao <- RunUMAP(Liao, reduction = "integrated.cca", dims = 1:30, reduction.name = "umap.cca")
DimPlot(Liao, reduction = "umap.cca", group.by = "orig.ident")

# Identify cell types with Azimuth
Liao <- JoinLayers(Liao)
Liao <- RunAzimuth(Liao, reference = "lungref", assay = "RNA")
DimPlot(Liao, reduction = "umap.cca", group.by = "predicted.ann_level_4", raster = FALSE) +
  theme_void()
Idents(Liao) <- "predicted.ann_level_4"

# Rename and Reorder identities alphabetically
Liao <- RenameIdents(Liao, 'Basal resting' = "Non-immun cells",
                    'Club' = "Non-immun cells",
                    'Deuterosomal' = "Non-immun cells",
                    'EC aerocyte capillary' = "Non-immun cells",
                    'EC general capillary' = "Non-immun cells",
                    'Goblet' = "Non-immun cells",
                    'Ionocyte' = "Non-immun cells",
                    'Multiciliated' = "Non-immun cells",
                    'None' = "Non-immun cells",
                    'SMG duct' = "Non-immun cells",
                    'SMG serous' = "Non-immun cells",
                    'Suprabasal' = "Non-immun cells",
                    'Transitional Club-AT2' = "Non-immun cells")
Liao$celltype <- Idents(Liao)
current_levels <- levels(Idents(Liao))
new_levels <- sort(current_levels)
Idents(Liao) <- factor(Idents(Liao), levels = new_levels)

# Define consistent cell type colors (T cells vs. myeloid clearly different)
standard_colors <- c(
  # T cells
  "CD4 T cells" = "#d73027",
  "CD8 T cells" = "#fc8d59",
  "T cells proliferating" = "#ef6548",
  "NK cells" = "#e34a33",
  
  # Myeloid
  "Alveolar macrophages" = "#4575b4",
  "Classical monocytes" = "#91bfdb",
  "Non-classical monocytes" = "#74add1",
  "Migratory DCs" = "#2c7bb6",
  "DC2" = "#313695",
  "Plasmacytoid DCs" = "#74a9cf",
  "Interstitial macrophages" = "#1F77B4",
  
  # Other
  "B cells" = "#999999",
  "Non-immun cells" = "#D9D9D9",
  "Plasma cells" = "#595959"
)

# Function to plot a standardized UMAP
plot_umap_standardized <- function(seurat_obj, title_text = "UMAP", label_clusters = FALSE) {
  DimPlot(seurat_obj, reduction = "umap.cca", label = label_clusters, raster = FALSE) +
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
plot_umap_standardized(Liao, title_text = "Human (Liao et al. 2020)")

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
plot_dot_standardized(Liao, title_text = "Human")

# Define myeloid cell types
myeloid_cell_types <- c(
  "Alveolar macrophages", "Classical monocytes", "DC2", "Interstitial macrophages",
  "Migratory DCs", "Non-classical monocytes", "Plasmacytoid DCs"
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
      fill = c("white", "firebrick"),  # Fill colors for each group
      color = c("firebrick", "black")  # Border colors
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
  seurat_obj = Liao,
  gene = "NOS2",
  group_col = "orig.ident",
  group_levels = c("CTR", "severe"),
  celltype_col = "predicted.ann_level_4"
)

# Export session info to a text file
sink("Liao_session_info.txt")
sessionInfo()
sink()