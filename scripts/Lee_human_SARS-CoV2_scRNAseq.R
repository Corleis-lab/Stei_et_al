options(future.globals.maxSize= 891289600)
library(tidyverse)
library(Seurat)
library(tidyseurat)
library(cowplot)
library(patchwork)
library(Azimuth)
theme_set(theme_bw())

# Load published data and create Seurat object
Lee.data <- Read10X(data.dir = "./data")
Lee <- CreateSeuratObject(counts = Lee.data, 
                                  project = "Lee_SARS-CoV2", 
                                  min.cells = 3, 
                                  min.features = 200)

# Define samples, severity and groups
Lee$Barcode <- rownames(Lee@meta.data)
Lee <- Lee %>% mutate(Sample = case_when(
  endsWith(Barcode, "1") ~ "nCoV1",
  endsWith(Barcode, "2") ~ "nCoV2",
  endsWith(Barcode, "3") ~ "Flu1",
  endsWith(Barcode, "4") ~ "Flu2",
  endsWith(Barcode, "5") ~ "CTR1",
  endsWith(Barcode, "6") ~ "Flu3",
  endsWith(Barcode, "7") ~ "Flu4",
  endsWith(Barcode, "8") ~ "Flu5",
  endsWith(Barcode, "9") ~ "nCoV3",
  endsWith(Barcode, "10") ~ "nCoV4",
  endsWith(Barcode, "20") ~ "nCoV11"
))

Lee <- Lee %>% mutate(Severity = case_when(
  endsWith(Barcode, "1") ~ "severe",
  endsWith(Barcode, "2") ~ "mild",
  endsWith(Barcode, "3") ~ "NA",
  endsWith(Barcode, "4") ~ "NA",
  endsWith(Barcode, "5") ~ "NA",
  endsWith(Barcode, "6") ~ "NA",
  endsWith(Barcode, "7") ~ "NA",
  endsWith(Barcode, "8") ~ "NA",
  endsWith(Barcode, "9") ~ "severe",
  endsWith(Barcode, "10") ~ "severe",
  endsWith(Barcode, "20") ~ "asymptomatic"
))

Lee <- Lee %>% mutate(Group = case_when(
  startsWith(Sample, "n") ~ "nCoV",
  startsWith(Sample, "F") ~ "Flu",
  startsWith(Sample, "C") ~ "CTR"
))

# Standard Seurat workflow
Lee <- NormalizeData(object = Lee)
Lee <- FindVariableFeatures(object = Lee)
Lee <- ScaleData(object = Lee)
Lee <- RunPCA(object = Lee)
Lee <- FindNeighbors(object = Lee, dims = 1:30)
Lee <- FindClusters(object = Lee, resolution = 0.5)
Lee <- RunUMAP(object = Lee, dims = 1:30)
DimPlot(Lee, reduction = "umap")

# Identify cell types with Azimuth
Lee <- RunAzimuth(Lee, reference = "lungref")
Idents(Lee) <- "predicted.ann_level_4"
DimPlot(Lee, reduction = "umap") +
  theme_void(base_size = 25)

# Rename and Reorder identities alphabetically
Lee <- RenameIdents(Lee, 'Alveolar fibroblasts' = "Non-immune cells",
                    'AT2 proliferating' = "Non-immune cells",
                    'Club' = "Non-immune cells",
                    'EC aerocyte capillary' = "Non-immune cells",
                    'EC general capillary' = "Non-immune cells",
                    'EC venous systemic' = "Non-immune cells",
                    'Ionocyte' = "Non-immune cells",
                    'Multiciliated' = "Non-immune cells",
                    'None' = "Non-immune cells",
                    'Pericytes' = "Non-immune cells",
                    'Suprabasal' = "Non-immune cells")
Lee$celltype <- Idents(Lee)
current_levels <- levels(Idents(Lee))
new_levels <- sort(current_levels)
Idents(Lee) <- factor(Idents(Lee), levels = new_levels)

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
  "DC1" = "#2c7bb6",
  "DC2" = "#313695",
  "Plasmacytoid DCs" = "#74a9cf",
  "Interstitial macrophages" = "#1F77B4",
  
  # Other
  "B cells" = "#999999",
  "Non-immune cells" = "#D9D9D9",
  "Plasma cells" = "#595959"
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
plot_umap_standardized(Lee, title_text = "Human (Lee et al. 2020)")

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
plot_dot_standardized(Lee, title_text = "Human")

# Define myeloid cell types
myeloid_cell_types <- c(
  "Alveolar macrophages", "Classical monocytes", "DC1", "DC2",
  "Interstitial macrophages", "Non-classical monocytes", "Plasmacytoid DCs"
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
  seurat_obj = Lee,
  gene = "NOS2",
  group_col = "Group",
  group_levels = c("CTR", "nCoV"),
  celltype_col = "predicted.ann_level_4"
)

# Export session info to a text file
sink("Lee_session_info.txt")
sessionInfo()
sink()