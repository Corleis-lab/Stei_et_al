library(Seurat)
library(ggplot2)
library(dplyr)
library(patchwork)
library(forcats)
library(biomaRt)

# Custom plotting theme
theme_pub <- theme_minimal(base_size = 16) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.text = element_text(face = "bold"),
    panel.grid.minor = element_blank()
  )

# Load preprocessed Seurat object
load("path/to/lung_filtered.rds")  # Replace with actual path

# Optional: load gene map for gene name conversion
gene_map <- read.csv("path/to/gene_map.csv")  # Columns: ensembl_gene_id, hgnc_symbol
rev_lookup <- setNames(gene_map$ensembl_gene_id, gene_map$hgnc_symbol)

# Cell type simplification
epi_map <- c(
  "pulmonary alveolar type 1 cell" = "AT1",
  "pulmonary alveolar type 2 cell" = "AT2",
  "bronchial goblet cell" = "secretory-like",
  "tracheobronchial goblet cell" = "secretory-like",
  "club cell" = "secretory-like",
  "respiratory basal cell" = "basal-like",
  "multiciliated columnar cell of tracheobronchial tree" = "ciliated-like",
  "multiciliated epithelial cell" = "ciliated-like"
)

lung_filtered$cell_type_simple <- epi_map[lung_filtered$cell_type]

# UMAP with simplified types
DimPlot(lung_filtered, group.by = "cell_type_simple", label = FALSE, repel = TRUE) +
  labs(title = "Human Lung Cell Atlas") +
  theme_pub +
  theme(legend.text = element_text(size = 16))

# Dotplot of top markers + NOS2
top_genes <- c("MUC5AC", "BPIFB1", "KRT5", "TP63", "NOS2")
ens_ids <- rev_lookup[top_genes]
ens_ids <- ens_ids[!is.na(ens_ids) & ens_ids %in% rownames(lung_filtered)]
gene_names <- names(ens_ids)

DotPlot(lung_filtered, features = ens_ids, group.by = "cell_type_simple") +
  scale_y_discrete(limits = rev) +
  scale_size(range = c(1, 6)) +
  labs(title = "Human Lung Cell Atlas", x = "Gene", y = "Cell Type") +
  theme_pub +
  theme(axis.text.x = element_text(size = 12, angle = 45, hjust = 1))

# Featureplots
plots <- lapply(ens_ids, function(gene_id) {
  FeaturePlot(lung_filtered, features = gene_id, order = TRUE, pt.size = 0.3) +
    ggtitle(names(ens_ids)[which(ens_ids == gene_id)]) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"))
})
wrap_plots(plots, ncol = 3)

# Boxplot epithelial
epi_simplified <- c("AT1", "AT2", "secretory-like", "ciliated-like", "basal-like")
ordered_tissues <- c(
  "inferior turbinate", "trachea", "left or right main bronchus",
  "segmental bronchi", "lobular bronchi", "distal lobular airways",
  "parenchyma upper lobe", "parenchyma right middle lobe", "parenchyma lower lobe"
)

df_epi <- lung_filtered@meta.data %>%
  mutate(NOS2 = FetchData(lung_filtered, vars = "ENSG00000007171")[, 1]) %>%
  filter(cell_type_simple %in% epi_simplified, tissue_level_2 %in% ordered_tissues) %>%
  group_by(tissue_level_2, cell_type_simple) %>%
  mutate(label = paste0(sum(NOS2 > 0), "/", n())) %>%
  ungroup() %>%
  mutate(
    tissue_level_2 = factor(tissue_level_2, levels = ordered_tissues),
    cell_type_simple = factor(cell_type_simple, levels = epi_simplified)
  )

ggplot(df_epi, aes(x = cell_type_simple, y = NOS2)) +
  geom_boxplot(outlier.shape = NA, fill = "grey80") +
  geom_jitter(width = 0.25, size = 0.7, alpha = 0.5, color = "black") +
  geom_text(aes(label = label, y = 3.1), size = 4, fontface = "plain", color = "black") +
  facet_wrap(~ tissue_level_2, ncol = 3) +
  scale_y_continuous(limits = c(0, 3.2)) +
  labs(y = "NOS2 expression (raw)", x = NULL, title = "Human Lung Cell Atlas") +
  theme_pub +
  theme(strip.text = element_text(size = 16), panel.spacing = unit(1.2, "lines"))

# Boxplot myeloid
myeloid_types <- c("alveolar macrophage", "classical monocyte", "non-classical monocyte",
                   "elicited macrophage", "lung macrophage", "dendritic cell",
                   "plasmacytoid dendritic cell", "CD1c-positive myeloid dendritic cell",
                   "conventional dendritic cell")

df_mye <- lung_filtered@meta.data %>%
  mutate(NOS2 = FetchData(lung_filtered, vars = "ENSG00000007171")[, 1]) %>%
  filter(cell_type %in% myeloid_types, tissue_level_2 %in% ordered_tissues) %>%
  group_by(tissue_level_2, cell_type) %>%
  mutate(label = paste0(sum(NOS2 > 0), "/", n())) %>%
  ungroup() %>%
  mutate(
    tissue_level_2 = factor(tissue_level_2, levels = ordered_tissues),
    cell_type = factor(cell_type, levels = myeloid_types)
  )

ggplot(df_mye, aes(x = cell_type, y = NOS2)) +
  geom_boxplot(outlier.shape = NA, fill = "grey80") +
  geom_jitter(width = 0.25, size = 0.7, alpha = 0.5, color = "black") +
  geom_text(aes(label = label, y = 3.1), size = 4, fontface = "plain", color = "black") +
  facet_wrap(~ tissue_level_2, ncol = 3) +
  scale_y_continuous(limits = c(0, 3.2)) +
  labs(y = "NOS2 expression (raw)", x = NULL, title = "Human Lung Cell Atlas") +
  theme_pub +
  theme(strip.text = element_text(size = 16), panel.spacing = unit(1.2, "lines"))

# Export session info to a text file
sink("Lung_cell_atlas_session_info.txt")
sessionInfo()
sink()