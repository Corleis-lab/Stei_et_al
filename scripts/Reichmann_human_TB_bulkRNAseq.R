library(tidyverse)
library(DESeq2)
library(EnhancedVolcano)
library(biomaRt)
library(org.Hs.eg.db)
library(AnnotationDbi)

# Load published data
Reichmann <- read.delim("./GSE174443_raw_counts_GRCh38.p13_NCBI.tsv")
Reichmann.metadata <- read.delim("./Reichmann_metadata.txt")

# Set up DESeq2 data set
rownames(Reichmann) <- Reichmann$GeneID
Reichmann$GeneID <- NULL
Reichmann.dds <- DESeqDataSetFromMatrix(countData = Reichmann,
                                        colData = Reichmann.metadata,
                                        design = ~condition)
Reichmann.dds$condition <- factor(Reichmann.dds$condition, levels = c("CTR", "Sarcoidosis", "TB"))
Reichmann.dds$condition <- relevel(Reichmann.dds$condition, ref = "CTR")

# Standard DESeq2 workflow
Reichmann.dds <- DESeq(Reichmann.dds)
Reichmann.res <- results(Reichmann.dds)
resultsNames(Reichmann.dds)
Reichmann.Sarcoidosis_vs_CTR <- results(Reichmann.dds, name = "condition_Sarcoidosis_vs_CTR")
Reichmann.TB_vs_CTR <- results(Reichmann.dds, name = "condition_TB_vs_CTR")
summary(Reichmann.Sarcoidosis_vs_CTR)
summary(Reichmann.TB_vs_CTR)

# Set up data frame for the volcano plot
Reichmann.TB_vs_CTR.df <- as.data.frame(Reichmann.TB_vs_CTR)
Reichmann.TB_vs_CTR.df$entrezgene_id <- rownames(Reichmann.TB_vs_CTR.df)
Reichmann.TB_vs_CTR.df$entrezgene_id <- as.integer(Reichmann.TB_vs_CTR.df$entrezgene_id)

# Translate NCBI gene IDs to gene symbols
ensdb_hm <- useEnsembl(biomart="genes", dataset="hsapiens_gene_ensembl", version=105)
entrezids <- Reichmann.TB_vs_CTR.df$entrezgene_id
symbols <- getBM(attributes = c("ensembl_gene_id", "entrezgene_id", "hgnc_symbol"),
                 filters = "entrezgene_id",
                 values = entrezids,
                 mart = ensdb_hm)
Reichmann.TB_vs_CTR.df <- left_join(Reichmann.TB_vs_CTR.df, symbols, by = "entrezgene_id")

# Generate volcano plot
EnhancedVolcano(Reichmann.TB_vs_CTR.df,
                lab = Reichmann.TB_vs_CTR.df$hgnc_symbol,
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
                col = c('black', 'pink', 'purple', 'firebrick'),
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
norm_counts <- counts(Reichmann.dds, normalized = TRUE)

# Map gene symbols again if needed
gene_symbols <- mapIds(org.Hs.eg.db,
                       keys = rownames(norm_counts),
                       column = "SYMBOL",
                       keytype = "ENTREZID",
                       multiVals = "first")

# Add symbols to rownames for easier subsetting
rownames(norm_counts) <- gene_symbols

# Select genes of interest
genes_of_interest <- c("CXCL10", "IFNG", "NOS2")

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
sample_conditions <- colData(Reichmann.dds) %>%
  as.data.frame() %>%
  dplyr::select(condition) %>%
  mutate(sample = rownames(.))

# Join with expression data
expr_long <- expr_long %>%
  left_join(sample_conditions, by = "sample")

# Check levels and keep conditions of interest
expr_long <- expr_long %>%
  filter(condition %in% c("CTR", "TB"))

# Generate box plot
ggplot(expr_long, aes(x = condition, y = expression)) +
  geom_boxplot(fill="firebrick", color="black") +
  geom_jitter(width = 0.2, alpha = 0.5) +
  facet_wrap(~ gene, scales = "free_y") +
  labs(y = "expression") +
  theme_bw() +
  theme(axis.title = element_text(size = 30),
        axis.text = element_text(size = 25),
        axis.text.x = element_text(hjust = 1, angle = 45),
        strip.text.x = element_text(size = 25, face = "bold"))

# Export session info to a text file
sink("Reichmann_session_info.txt")
sessionInfo()
sink()