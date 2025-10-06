library(tidyverse)
library(DESeq2)
library(EnhancedVolcano)
library(org.Mm.eg.db)
library(AnnotationDbi)

# Load published data
lung <- read.delim("./mouse_lung_raw_data.txt")
lung.metadata <- read.delim("./lung_metadata.txt")

# Set up DESeq2 data set
lung.dds <- DESeqDataSetFromMatrix(countData = lung,
                                    colData = lung.metadata,
                                    design = ~condition, tidy = TRUE)
lung.dds$condition <- factor(lung.dds$condition, levels = c("CTR", "low", "high"))
lung.dds$condition <- relevel(lung.dds$condition, ref = "CTR")

# Standard DESeq2 workflow
lung.dds <- DESeq(lung.dds)
lung.res <- results(lung.dds)
resultsNames(lung.dds)
lung.low_vs_CTR <- results(lung.dds, name = "condition_low_vs_CTR")
lung.high_vs_CTR <- results(lung.dds, name = "condition_high_vs_CTR")
summary(lung.low_vs_CTR)
summary(lung.high_vs_CTR)

# Generate volcano plot
EnhancedVolcano(lung.high_vs_CTR.df,
                lab = lung.high_vs_CTR.df$symbol,
                x = 'log2FoldChange',
                y = 'padj',
                xlab = bquote(~Log[2]~ 'fold change'),
                pCutoff = 0.05,
                FCcutoff = 1.0,
                pointSize = 3,
                labSize = 8,
                axisLabSize = 25,
                labCol = "black",
                labFace = "bold",
                col = c('black', 'pink', 'purple', 'lightskyblue'),
                colAlpha = 4/5,
                legendPosition = 'none',
                subtitle = NULL,
                caption = NULL,
                drawConnectors = FALSE,
                widthConnectors = 1.0,
                arrowheads = FALSE,
                #title = 'high dose TB vs. CTR - lung',
                title = NULL,
                xlim = c(-10, 15))

# Get normalized counts
norm_counts <- counts(lung.dds, normalized = TRUE)

# Map gene symbols again if needed
gene_symbols <- mapIds(org.Mm.eg.db,
                       keys = rownames(norm_counts),
                       column = "SYMBOL",
                       keytype = "ENSEMBL",
                       multiVals = "first")

# Add symbols to rownames for easier subsetting
rownames(norm_counts) <- gene_symbols

# Select genes of interest
genes_of_interest <- c("Cxcl10", "Ifng", "Nos2")

# Subset and reshape to long format
expr_long <- norm_counts[rownames(norm_counts) %in% genes_of_interest, ] %>%
  as.data.frame() %>%
  mutate(gene = rownames(.)) %>%
  pivot_longer(
    cols = -gene,
    names_to = "sample",
    values_to = "expression"
  )

# Extract condition info
sample_conditions <- colData(lung.dds) %>%
  as.data.frame() %>%
  dplyr::select(condition) %>%
  mutate(sample = rownames(.))

# Join with expression data
expr_long <- expr_long %>%
  left_join(sample_conditions, by = "sample")

# Generate box plot
ggplot(expr_long, aes(x = condition, y = expression)) +
  geom_boxplot(fill="lightskyblue", color="black") +
  geom_jitter(width = 0.2, alpha = 0.5) +
  facet_wrap(~ gene, scales = "free_y") +
  labs(y = "expression") +
  theme_bw() +
  theme(axis.title = element_text(size = 30),
        axis.text = element_text(size = 25),
        axis.text.x = element_text(hjust = 1, angle = 45),
        strip.text.x = element_text(size = 25, face = "bold"))

# Export session info to a text file
sink("Moreira-Teixeira_session_info.txt")
sessionInfo()
sink()