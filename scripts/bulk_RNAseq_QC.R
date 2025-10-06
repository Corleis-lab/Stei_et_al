### Input reads per sample and read mapping summary ###

library(tidyverse)

# Run this script twice (once for mouse and once for human)
# Set directory containing STAR logs
log_dir <- "./logs/STAR/mouse" # change path for human data
log_files <- list.files(log_dir, pattern = "Log.final.out$", full.names = TRUE)

# Helper to parse a single STAR log
parse_star_log <- function(file) {
  lines <- readLines(file)
  
  # Split each line at the pipe symbol
  tab <- strsplit(lines, "\\|")
  
  # Filter only lines that split into 2 parts (metric + value)
  parsed <- tab[sapply(tab, length) == 2]
  
  # Trim whitespace
  metrics <- trimws(sapply(parsed, `[`, 1))
  values  <- trimws(sapply(parsed, `[`, 2))
  
  # Build named vector or tibble
  data <- setNames(values, metrics)
  
  # Return selected fields
  tibble(
    sample = gsub("Log.final.out$", "", basename(file)),
    input_reads = as.numeric(data["Number of input reads"]),
    uniquely_mapped_pct = as.numeric(sub("%", "", data["Uniquely mapped reads %"])),
    multimapped_pct = as.numeric(sub("%", "", data["% of reads mapped to multiple loci"])),
    unmapped_short_pct = as.numeric(sub("%", "", data["% of reads unmapped: too short"])),
    unmapped_other_pct = as.numeric(sub("%", "", data["% of reads unmapped: other"]))
  )
}

# Parse all logs
qc_metrics <- map_dfr(log_files, parse_star_log)

# Save QC table
write.table(qc_metrics, "STAR_log_summary.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

# Plot 1: Input reads per sample
p <- ggplot(qc_metrics, aes(x = fct_reorder(sample, input_reads), y = input_reads / 1e6)) +
  geom_col(fill = "lightskyblue") + # fill = "firebrick" for human
  theme_minimal(base_size = 12) +
  coord_flip() +
  labs(
    title = "Input Reads per Sample",
    x = "Sample",
    y = "Millions of Reads"
  ) +
  theme(
    axis.text.y = element_text(size = 8),
    plot.title = element_text(hjust = 0.5)
  )

ggsave("mouse_STAR_input_reads_barplot.pdf", p, width = 8, height = 0.25 * nrow(qc_metrics) + 2) # change file name for human

# Plot 2: Mapping breakdown
qc_long <- qc_metrics %>%
  pivot_longer(cols = ends_with("_pct"),
               names_to = "category",
               values_to = "percent") %>%
  mutate(category = recode(category,
                           "uniquely_mapped_pct" = "Unique",
                           "multimapped_pct" = "Multi",
                           "unmapped_short_pct" = "Unmapped (too short)",
                           "unmapped_other_pct" = "Unmapped (other)"))

p_map <- ggplot(qc_long, aes(x = fct_reorder(sample, -percent), y = percent, fill = category)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  scale_fill_brewer(palette = "Set2") +
  theme_minimal(base_size = 12) +
  labs(
    title = "Read Mapping Summary",
    x = "Sample",
    y = "Percent"
  ) +
  theme(
    axis.text.y = element_text(size = 8),
    legend.position = "bottom",
    plot.title = element_text(hjust = 0.5)
  )

ggsave("mouse_STAR_mapping_summary_barplot.pdf", p_map, width = 8, height = 0.25 * length(unique(qc_long$sample)) + 2) # change file name for human

# ===================================================================================================================================================

### Clean up column headers of count files ###
library(tidyverse)

# Base directory for count files
counts_dir <- "./counts"

# Species list
species_list <- c("mouse", "human")

# Helper function to clean column names
clean_colnames <- function(df) {
  colnames(df) <- colnames(df) %>%
    gsub("^.*/", "", .) %>%                    # remove path
    gsub("Aligned.sorted.bam$", "", .)        # remove BAM suffix
  return(df)
}

for (species in species_list) {
  message("Cleaning featureCounts output for: ", species)
  
  # TOTAL COUNTS
  total_file <- file.path(counts_dir, paste0(species, "_featureCounts.txt"))
  total_clean_file <- file.path(counts_dir, paste0(species, "_counts_clean.txt"))
  
  if (file.exists(total_file)) {
    total <- read.delim(total_file, comment.char = "#", check.names = FALSE)
    total_clean <- clean_colnames(total)
    write.table(total_clean, total_clean_file, sep = "\t", quote = FALSE, row.names = FALSE)
    message("leaned: ", total_clean_file)
  } else {
    warning("Missing file: ", total_file)
  }
  
  # rRNA COUNTS
  rrna_file <- file.path(counts_dir, paste0(species, "_rRNA_counts.txt"))
  rrna_clean_file <- file.path(counts_dir, paste0(species, "_rRNA_counts_clean.txt"))
  
  if (file.exists(rrna_file)) {
    rrna <- read.delim(rrna_file, comment.char = "#", check.names = FALSE)
    rrna_clean <- clean_colnames(rrna)
    write.table(rrna_clean, rrna_clean_file, sep = "\t", quote = FALSE, row.names = FALSE)
    message("Cleaned: ", rrna_clean_file)
  } else {
    warning("Missing file: ", rrna_file)
  }
}

# ===================================================================================================================================================

### Calculate rRNA reads ###
library(tidyverse)

# Species list
species_list <- c("mouse", "human")

# Species-specific colors
species_colors <- c("mouse" = "lightskyblue", "human" = "firebrick")

# Base directory
counts_dir <- "./counts"

# Loop through species
for (species in species_list) {
  message("Processing rRNA contamination for: ", species)
  
  # Filenames
  total_file <- file.path(counts_dir, paste0(species, "_counts_clean.txt"))
  rrna_file  <- file.path(counts_dir, paste0(species, "_rRNA_counts_clean.txt"))
  output_summary <- file.path(counts_dir, paste0(species, "_rRNA_percentages.tsv"))
  output_plot <- file.path(counts_dir, paste0(species, "_rRNA_percentages.pdf"))
  
  # Check if files exist
  if (!file.exists(total_file) || !file.exists(rrna_file)) {
    warning("Missing files for ", species, " - skipping.")
    next
  }
  
  # Load counts
  total <- read.delim(total_file, comment.char = "#", check.names = FALSE)
  rrna  <- read.delim(rrna_file, comment.char = "#", check.names = FALSE)
  
  # Extract count matrices (remove annotation columns)
  total_counts <- total[, 7:ncol(total)]
  rrna_counts  <- rrna[, 7:ncol(rrna)]
  
  # Sanity check: same sample names
  common_samples <- intersect(names(total_counts), names(rrna_counts))
  if (length(common_samples) == 0) {
    warning("No overlapping samples found in counts for ", species)
    next
  }
  
  total_sum <- colSums(total_counts[, common_samples])
  rrna_sum  <- colSums(rrna_counts[, common_samples])
  
  # Calculate percentage
  rRNA_pct <- (rrna_sum / total_sum) * 100
  rRNA_df <- tibble(Sample = names(rRNA_pct), rRNA_percent = round(rRNA_pct, 2))
  
  # Save summary
  write.table(rRNA_df, output_summary, sep = "\t", quote = FALSE, row.names = FALSE)
  
  # Plot with species-specific color
  plot <- ggplot(rRNA_df, aes(x = reorder(Sample, -rRNA_percent), y = rRNA_percent)) +
    geom_bar(stat = "identity", fill = species_colors[[species]]) +
    theme_minimal() +
    labs(
      title = paste("Estimated rRNA Contamination —", species),
      x = "Sample",
      y = "rRNA Reads (%)"
    ) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  ggsave(output_plot, plot, width = 8, height = 4)
  message("Saved summary: ", output_summary)
}

# ===================================================================================================================================================

### Genes per sample ###

library(DESeq2)
library(edgeR)
library(ggplot2)
library(reshape2)

# Parameters
files <- c("./counts/human_counts_clean.txt", "./counts/mouse_counts_clean.txt")
dataset_names <- gsub("(_)?counts_clean(er)?\\.txt$", "", basename(files))
output_dir <- "RNAseq_QC_Plots"
dir.create(output_dir, showWarnings = FALSE)
species_colors <- c("mouse" = "lightskyblue", "human" = "firebrick")

# Loop Through Files
for (i in seq_along(files)) {
  file <- files[i]
  dataset <- dataset_names[i]
  cat("Processing:", dataset, "\n")
  
  # Load & Clean Count Matrix
  counts_raw <- read.delim(file, header = TRUE, row.names = 1, check.names = FALSE)
  count_matrix <- counts_raw[, 6:ncol(counts_raw)]
  
  # Parse Sample Metadata
  samples <- colnames(count_matrix)
  sample_info <- data.frame(
    sample = sub("_.*", "", samples),
    condition = sub(".*_", "", samples)
  )
  rownames(sample_info) <- samples
  sample_info$sample <- factor(sample_info$sample)
  sample_info$condition <- factor(sample_info$condition)
  
  # Create DESeq2 Dataset
  dds <- DESeqDataSetFromMatrix(countData = count_matrix,
                                colData = sample_info,
                                design = ~ condition)
  dds <- DESeq(dds)
  res <- results(dds)
  
  # Create Output Folder
  dataset_dir <- file.path(output_dir, dataset)
  dir.create(dataset_dir, showWarnings = FALSE)
  
  # Plot 1: Detected Genes Per Sample
  dge <- DGEList(counts = count_matrix)
  cpm_vals <- cpm(dge)
  detected_genes <- colSums(cpm_vals > 1)
  
  # Identify species from dataset name
  species <- tolower(gsub("_.*", "", dataset))  # crude parse from filename
  
  # Bar color based on species
  bar_color <- species_colors[[species]]
  if (is.null(bar_color)) bar_color <- "gray70"  # fallback
  
  p1 <- ggplot(data.frame(Sample = names(detected_genes), Genes = detected_genes), 
               aes(x = Sample, y = Genes)) +
    geom_bar(stat = "identity", fill = bar_color) +
    theme_minimal() +
    labs(title = paste("Detected Genes per Sample:", dataset),
         y = "Number of Genes", x = "Sample") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  ggsave(file.path(dataset_dir, "detected_genes_per_sample.pdf"), p1, width = 8, height = 5)
}