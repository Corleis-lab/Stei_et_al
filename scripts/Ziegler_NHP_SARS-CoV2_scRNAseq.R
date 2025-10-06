library(ggplot2)
library(dplyr)
library(readr)
library(tidyr)
library(grid)  # for unit()

# Set file paths
base_dir <- "/home/RNAadmin/SCP806"
umap_file <- file.path(base_dir, "cluster/epi_cells_umap.txt")
metadata_file <- file.path(base_dir, "metadata/epithelial_metadata_alexandria.txt")
expression_file <- file.path(base_dir, "expression/epithelial_cells_for_scp.txt.gz")

# Load UMAP and metadata
umap_data <- read_tsv(umap_file, show_col_types = FALSE)
colnames(umap_data) <- c("NAME", "UMAP_1", "UMAP_2")

metadata <- read_tsv(metadata_file, show_col_types = FALSE)

# Merge metadata and UMAP, remove invalid label
merged_data <- inner_join(umap_data, metadata, by = "NAME") %>%
  filter(cell_type__ontology_label != "group")

# UMAP Plot
umap_plot <- ggplot(merged_data, aes(x = UMAP_1, y = UMAP_2, color = cell_type__ontology_label)) +
  geom_point(size = 1.2, alpha = 0.85) +
  ggtitle("UMAP of Epithelial Cells") +
  labs(color = "Cell Type") +
  theme_void() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 18),
    legend.position = "right",
    legend.box.margin = margin(0, 20, 0, 0),
    legend.title = element_text(size = 20),
    legend.text = element_text(size = 18),
    legend.key.size = unit(0.8, "cm"),
    plot.margin = unit(c(1, 1, 1, 1), "cm")
  ) +
  guides(color = guide_legend(override.aes = list(size = 5)))

print(umap_plot)


# NOS2 Feature Plot
expr_matrix <- read_tsv(expression_file, show_col_types = FALSE)
colnames(expr_matrix)[1] <- "GENE"

# Extract NOS2 and reshape
nos2_expr <- expr_matrix %>%
  filter(GENE == "NOS2") %>%
  pivot_longer(-GENE, names_to = "NAME", values_to = "NOS2")

# Merge expression
feature_data <- merged_data %>%
  left_join(nos2_expr, by = "NAME") %>%
  mutate(NOS2 = replace_na(NOS2, 0))

# Plot NOS2 on UMAP
nos2_plot <- ggplot(feature_data, aes(x = UMAP_1, y = UMAP_2, color = NOS2)) +
  geom_point(size = 1.2, alpha = 0.9) +
  scale_color_gradient(low = "gray90", high = "red", trans = "log1p", name = "NOS2") +
  ggtitle("NOS2 Expression (log scale)") +
  theme_void() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 14),
    legend.key.height = unit(0.6, "cm"),
    legend.key.width = unit(0.4, "cm"),
    plot.margin = unit(c(1, 1, 1, 1), "cm")
  )

print(nos2_plot)


# Bar Plot: NOS2+ Secretory Cells

secretory_data <- feature_data %>%
  filter(cell_type__ontology_label == "secretory cell") %>%
  filter(Granuloma %in% c("Granuloma", "Uninvolved lung")) %>%
  mutate(NOS2_positive = NOS2 > 0)

summary_df <- secretory_data %>%
  group_by(Granuloma) %>%
  summarise(
    total_cells = n(),
    NOS2_positive_cells = sum(NOS2_positive),
    .groups = "drop"
  ) %>%
  mutate(label_text = paste0(NOS2_positive_cells, " / ", total_cells))

bar_plot <- ggplot(summary_df, aes(x = Granuloma, y = NOS2_positive_cells)) +
  geom_col(fill = c("gray80", "firebrick"), color = "black", width = 0.6) +
  geom_text(aes(label = label_text), vjust = -0.5, size = 6, fontface = "bold") +
  theme_minimal(base_size = 16) +
  theme(
    plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
    axis.title.x = element_blank(),
    axis.title.y = element_text(face = "bold", size = 16),
    axis.text.x = element_text(size = 14, face = "bold"),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
    plot.margin = margin(10, 10, 10, 10)
  ) +
  labs(y = "NOS2+ Cell Count", title = "NOS2+ Secretory Cells")

print(bar_plot)

# Export session info to a text file
sink("Ziegler_NHP_session_info.txt")
sessionInfo()
sink()