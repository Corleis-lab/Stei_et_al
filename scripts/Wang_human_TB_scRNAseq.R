# reanalyze_Wang_GSE192483.R
# ------------------------------------------------------------------------------

####################################################################################################################
## Sort published data into subfolders and rename files (Windows Powershell)
#cd "C:\path\to\your\files"
#
#$prefixes = Get-ChildItem -File | ForEach-Object {
#  ($_ -replace "(\.barcodes\.tsv\.gz|\.features\.tsv\.gz|\.matrix\.mtx\.gz|\.filtered_feature_bc_matrix\.h5)$","")
#} | Sort-Object -Unique
#
#foreach ($prefix in $prefixes) {
#  # Create subfolder with the prefix name
#  New-Item -ItemType Directory -Name $prefix -Force | Out-Null
#  
#  # Get all files with this prefix
#  $files = Get-ChildItem -File | Where-Object { $_.Name -like "$prefix*" }
#  
#  foreach ($file in $files) {
#    # Strip the "prefix." part (everything up to and including the first dot)
#    $newName = $file.Name -replace "^$prefix\.",""
#    
#    # Move with cleaned name
#    Move-Item $file.FullName (Join-Path $prefix $newName)
#  }
#}
####################################################################################################################

library(tidyverse)
library(Seurat)
library(tidyseurat)
library(cowplot)
library(patchwork)
library(sctransform)
theme_set(theme_bw())

# =========================
# 1) Load data + metadata
# =========================

# Adjust if needed
low_dir  <- "./data/Wang_GSE192483/Low"
high_dir <- "./data/Wang_GSE192483/High"

# Helper: read one 10x sample directory into a Seurat object
.read_one_10x <- function(sample_dir, add_group) {
  candidates <- c(
    file.path(sample_dir, "filtered_feature_bc_matrix"),
    file.path(sample_dir, "filtered_gene_bc_matrices", "GRCh38"),
    sample_dir
  )
  tenx_dir <- candidates[file.exists(candidates)][1]
  if (is.na(tenx_dir)) stop("No 10x folder found under: ", sample_dir)
  
  dat <- Read10X(data.dir = tenx_dir, gene.column = 1)
  obj <- CreateSeuratObject(counts = dat,
                            project = basename(sample_dir),
                            min.cells = 3,
                            min.features = 200)
  obj$sample <- basename(sample_dir)
  obj$Group  <- add_group  # "Low" or "High"
  obj
}

# Read all samples within a group folder (each subdir is one sample)
.read_group <- function(group_dir, group_label) {
  subs <- list.dirs(group_dir, full.names = TRUE, recursive = FALSE)
  if (length(subs) == 0) subs <- group_dir
  
  objs <- lapply(subs, .read_one_10x, add_group = group_label)
  
  if (length(objs) == 1) {
    return(objs[[1]])
  } else {
    merge(x = objs[[1]], 
          y = objs[-1], 
          add.cell.ids = sapply(objs, function(o) o$sample[1]), 
          merge.dr = TRUE)
  }
}

Wang_Low  <- .read_group(low_dir,  "Low")
Wang_High <- .read_group(high_dir, "High")

# =========================
# 2) Standard workflows
# =========================

# Low
Wang_Low <- NormalizeData(Wang_Low)
Wang_Low <- FindVariableFeatures(Wang_Low)
Wang_Low <- ScaleData(Wang_Low)
Wang_Low <- RunPCA(Wang_Low)
Wang_Low <- FindNeighbors(Wang_Low, dims = 1:30)
Wang_Low <- FindClusters(Wang_Low)
Wang_Low <- RunUMAP(Wang_Low, dims = 1:30)

# High
Wang_High <- NormalizeData(Wang_High)
Wang_High <- FindVariableFeatures(Wang_High)
Wang_High <- ScaleData(Wang_High)
Wang_High <- RunPCA(Wang_High)
Wang_High <- FindNeighbors(Wang_High, dims = 1:30)
Wang_High <- FindClusters(Wang_High)
Wang_High <- RunUMAP(Wang_High, dims = 1:30)

# =========================
# 3) Merge + integrate
# =========================

Wang <- merge(Wang_Low, Wang_High, add.cell.ids = c("Low","High"), merge.dr = TRUE)

# Keep the imported Group ("Low"/"High")
Wang$Group <- Wang$Group %>% factor(levels = c("Low","High"))

# Convert Ensembl IDs to gene symbols
ConvertEnsemblToSymbol <- function(seurat_obj) {
  if (!requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
    BiocManager::install("org.Hs.eg.db")
  }
  library(org.Hs.eg.db)
  library(Matrix)
  
  # Map Ensembl → Symbol
  gene_map <- AnnotationDbi::mapIds(
    org.Hs.eg.db,
    keys = rownames(seurat_obj),
    column = "SYMBOL",
    keytype = "ENSEMBL",
    multiVals = "first"   # if multiple symbols, pick the first
  )
  
  # Replace rownames, keep Ensembl if no symbol found
  new_names <- ifelse(is.na(gene_map), rownames(seurat_obj), gene_map)
  
  # Extract counts explicitly from "counts" layer (Seurat v5 safe)
  counts <- GetAssayData(seurat_obj, assay = DefaultAssay(seurat_obj), layer = "counts")
  counts <- as.matrix(counts)
  rownames(counts) <- new_names
  
  # Collapse duplicates by summing rows with the same symbol
  counts_collapsed <- rowsum(counts, group = new_names)
  
  # Rebuild Seurat object
  new_obj <- CreateSeuratObject(
    counts = counts_collapsed,
    meta.data = seurat_obj@meta.data,
    assay = DefaultAssay(seurat_obj)
  )
  
  return(new_obj)
}


Wang <- ConvertEnsemblToSymbol(Wang)

# Recompute on merged
Wang <- NormalizeData(Wang)
Wang <- FindVariableFeatures(Wang)
Wang <- ScaleData(Wang)
Wang <- RunPCA(Wang)
Wang <- FindNeighbors(Wang, dims = 1:30, reduction = "pca")
Wang <- FindClusters(Wang, resolution = 0.25)
Wang <- RunUMAP(Wang, dims = 1:30)
Wang <- JoinLayers(Wang)

# Find cluster markers
Wang.markers <- FindAllMarkers(Wang, only.pos = TRUE)
Wang.markers <- Wang.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC>1)

top50 <- Wang.markers %>%
  mutate(score = avg_log2FC * log10(1 / (p_val + 1e-300)) * pct.1) %>%
  group_by(cluster) %>%
  slice_max(order_by = score, n = 50, with_ties = FALSE) %>%
  ungroup() %>%
  arrange(cluster, desc(score))

gene_table <- top50 %>%
  dplyr::select(cluster, gene) %>%
  group_by(cluster) %>%
  mutate(rank = row_number()) %>%  # rank 1..50 inside each cluster
  ungroup() %>%
  pivot_wider(
    names_from = cluster,
    values_from = gene
  )

gene_table <- gene_table %>% dplyr::select(-rank)

# =========================
# 4) Harmonize labels
# =========================

# If needed, map to consistent names (edit as appropriate for Wang)
Wang <- RenameIdents(Wang,
                     '0' = "CD4 T cells",
                     '1' = "CD8 T cells",
                     '2' = "Macrophages",
                     '3' = "Macrophages",
                     '4' = "NK / Cytotoxic T cells",
                     '5' = "B cells",
                     '6' = "Gamma delta T cells",
                     '7' = "Monocytes (non-classical)",
                     '8' = "B cells",
                     '9' = "non-immune cells",
                     '10' = "Regulatory T cells (Tregs)",
                     '11' = "Monocytes (classical)",
                     '12' = "non-immune cells",
                     '13' = "B cells",
                     '14' = "Cycling T cells",
                     '15' = "Dendritic cells",
                     '16' = "Mast Cells",
                     '17' = "non-immune cells",
                     '18' = "Cycling T cells",
                     '19' = "Dendritic cells",
                     '20' = "B cells"
)
Wang$celltype_new <- Idents(Wang)
Idents(Wang) <- factor(Idents(Wang), levels = sort(levels(Idents(Wang))))
DimPlot(Wang, reduction = "umap")

# =========================
# 5) Colors + plot helpers
# =========================

standard_colors <- c(
  # T/innate-like
  "CD4 T cells" = "#d73027",
  "CD8 T cells" = "#f4a582",
  "Regulatory T cells (Tregs)" = "#ca0020",
  "Cycling T cells" = "#fc8d59",
  "NK / Cytotoxic T cells" = "#b2182b",
  "Gamma delta T cells" = "#ef8a62",
  
  # Myeloid
  "Macrophages" = "#4575b4",
  "Monocytes (classical)" = "#74add1",
  "Monocytes (non-classical)" = "#2c7bb6",
  "Dendritic cells" = "#92c5de",
  "cDC" = "#74add1",
  "pDC" = "#2c7bb6",
  "Neutrophils" = "#313695",
  "Mast cells" = "#74a9cf",
  
  # Other
  "B cells" = "#999999",
  "Endothelial cells" = "#cccccc",
  "Fibroblasts" = "#dddddd",
  "Plasma cells" = "#bbbbbb",
  "non-immune cells" = "#D9D9D9"
)

plot_umap_standardized <- function(Wang, title_text = "Human (Wang et al. 2024)", label_clusters = FALSE) {
  DimPlot(Wang, reduction = "umap", label = label_clusters, raster = FALSE) +
    scale_color_manual(values = standard_colors, na.translate = FALSE) +
    ggtitle(title_text) +
    theme_void() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 24),
      legend.text = element_text(size = 20),
      legend.title = element_text(size = 20, face = "bold"),
      plot.margin = margin(10, 10, 10, 10)
    )
}

plot_dot_standardized <- function(Wang,
                                  features = c("CXCL10","IFNG","NOS2"),
                                  title_text = "Human") {
  DotPlot(Wang, features = features) +
    scale_color_gradientn(colors = c("grey90","purple","blue"), name = "Average Expression") +
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
      size   = guide_legend(order = 2, title = "Percent Expressed")
    )
}

# Myeloid set used for the boxplot
myeloid_cell_types <- c("Macrophages","Monocytes (classical)","Monocytes (non-classical)",
                        "Dendritic cells","cDC","pDC","Neutrophils","Mast cells")

plot_myeloid_bin_boxplot <- function(Wang,
                                     gene = "NOS2",
                                     group_col = "Group",
                                     group_levels = c("Low","High"),
                                     celltype_col = "celltype_new",
                                     myeloid_types = myeloid_cell_types) {
  
  missing_cols <- setdiff(c(group_col, celltype_col), colnames(Wang@meta.data))
  if (length(missing_cols) > 0) stop("Missing metadata column(s): ", paste(missing_cols, collapse = ", "))
  
  expr <- GetAssayData(Wang, layer = "data")
  if (!(gene %in% rownames(expr))) stop("Gene not found in expression matrix: ", gene)
  
  df <- tibble(
    Cell       = colnames(Wang),
    Group      = Wang@meta.data[[group_col]],
    Celltype   = Wang@meta.data[[celltype_col]],
    Expression = as.numeric(expr[gene, colnames(Wang)])
  ) %>%
    filter(Group %in% group_levels, Celltype %in% myeloid_types) %>%
    mutate(Group = factor(Group, levels = group_levels))
  
  if (nrow(df) == 0) stop("No cells after filtering; check group names and cell types.")
  
  ggplot(df, aes(x = Group, y = Expression, fill = Group, color = Group)) +
    geom_boxplot(width = 0.6) +
    scale_fill_manual(values = c("Low" = "white", "High" = "darkorange")) +
    scale_color_manual(values = c("Low" = "darkorange", "High" = "black")) +
    theme_minimal(base_size = 16) +
    theme(
      plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
      axis.title.x = element_blank(),
      axis.title.y = element_text(face = "bold"),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
      plot.margin = margin(5,5,5,5),
      axis.text.x = element_text(size = 14, face = "bold")
    ) +
    labs(y = "Expression", title = gene)
}

# =========================
# 6) Generate plots
# =========================

Idents(Wang) <- "celltype_new"

# UMAP
p_umap <- plot_umap_standardized(Wang, title_text = "Human (Wang et al. 2024)")

# Dot plot
p_dot  <- plot_dot_standardized(Wang, features = c("CXCL10","IFNG","NOS2"), title_text = "Human")

# NOS2 boxplot for myeloid; Low vs High
p_box  <- plot_myeloid_bin_boxplot(Wang, gene = "NOS2", group_col = "Group", group_levels = c("Low","High"))

# =========================
# 7) Save outputs + session
# =========================

dir.create("plots_Wang_GSE192483", showWarnings = FALSE)
ggsave("plots_Wang_GSE192483/UMAP_Wang_GSE192483.png", p_umap, width = 7, height = 6, dpi = 300)
ggsave("plots_Wang_GSE192483/Dotplot_CXCL10_IFNG_NOS2_Wang_GSE192483.png", p_dot, width = 7, height = 6, dpi = 300)
ggsave("plots_Wang_GSE192483/Boxplot_NOS2_Low_vs_High_Myeloid_Wang_GSE192483.png", p_box, width = 4.5, height = 5.5, dpi = 300)

sink("plots_Wang_GSE192483/Wang_GSE192483_session_info.txt")
sessionInfo()
sink()