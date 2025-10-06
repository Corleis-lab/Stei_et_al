library(DESeq2)
library(tibble)
library(stringr)
library(purrr)
library(GSEA)

# Function to Create DESeq2 Object and Write FPM
process_counts <- function(directory, species, condition_vector, output_fpm_file) {
  message("Processing ", species, "...")
  
  sample_files <- list.files(directory, pattern = ".txt", full.names = FALSE)
  sample_names <- str_remove(sample_files, "_counts.txt")
  
  if (length(condition_vector) != length(sample_files)) {
    stop("Condition vector length does not match number of sample files for ", species)
  }
  
  sample_table <- data.frame(
    sampleName = sample_names,
    fileName = sample_files,
    condition = factor(condition_vector)
  )
  
  dds <- DESeqDataSetFromHTSeqCount(
    sampleTable = sample_table,
    directory = directory,
    design = ~ condition
  )
  
  dds <- dds[rowSums(counts(dds)) >= 10, ]
  fpm_matrix <- round(fpm(dds, robust = FALSE), 3)
  
  write.table(fpm_matrix, file = output_fpm_file, quote = FALSE, sep = "\t")
  return(invisible(fpm_matrix))
}

# Function to Run GSEA Analysis
run_gsea <- function(input_gct, input_cls, input_chip, gs_db, output_dir, doc_string) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  
  GSEA(
    input.ds = input_gct,
    input.cls = input_cls,
    input.chip = input_chip,
    collapse.dataset = TRUE,
    collapse.mode = "sum",
    gs.db = gs_db,
    output.directory = output_dir,
    doc.string = doc_string,
    gsea.type = "GSEA"
  )
}

# Human: Count processing
human_conditions <- c("IFNy", "Mtb", "IFNyMtb", "CTR", "IFNy", "Mtb", "IFNyMtb", "CTR", "IFNy", "Mtb", "IFNyMtb", "CTR", "IFNy", "Mtb", "IFNyMtb", "CTR", "IFNy", "Mtb", "IFNyMtb")
human_fpm <- process_counts(
  directory = "./human_count_data (wo H1.1)",
  species = "human",
  condition_vector = human_conditions,
  output_fpm_file = "./human_fpm.txt"
)

# Mouse: Count processing
mouse_conditions <- c("CTR", "Mtb", "IFNyMtb", "IFNy", "Mtb", "IFNyMtb", "CTR", "IFNy", "Mtb", "IFNyMtb", "CTR", "IFNy", "Mtb", "IFNyMtb", "CTR", "IFNy", "Mtb", "IFNyMtb")
mouse_fpm <- process_counts(
  directory = "./mouse_count_data",
  species = "mouse",
  condition_vector = mouse_conditions,
  output_fpm_file = "./mouse_fpm.txt"
)

# GSEA Parameters
gsea_runs <- tribble(
  ~species, ~input_gct, ~input_cls, ~input_chip, ~gs_db, ~output_dir, ~doc_string,
  "mouse", "./exp_files/mouse_fpm.gct", "./cls_files/mouse.cls", "./chip_files/Mouse_Gene_Symbol_Remapping_MSigDB.v2024.1.Mm.chip", "./pathways/mouse/GOBP_ACUTE_INFLAMMATORY_RESPONSE.v2024.1.Mm.gmt", "./results/mouse_ACUTE_INFLAMMATORY_RESPONSE", "mouse",
  "mouse", "./exp_files/mouse_fpm.gct", "./cls_files/mouse.cls", "./chip_files/Mouse_Gene_Symbol_Remapping_MSigDB.v2024.1.Mm.chip", "./pathways/mouse/REACTOME_INTERFERON_GAMMA_SIGNALING.v2024.1.Mm.gmt", "./results/mouse_INTERFERON_GAMMA_SIGNALING", "mouse",
  "mouse", "./exp_files/mouse_orthologes.gct", "./cls_files/mouse.cls", "./chip_files/Human_Gene_Symbol_with_Remapping_MSigDB.v2024.1.Hs.chip", "./pathways/orthologes/ZAMORA_NOS2_TARGETS_UP.v2024.1.Hs.gmt", "./results/mouse_orthologes_NOS2_TARGETS", "mouse_orthologes",
  "mouse", "./exp_files/mouse_fpm.gct", "./cls_files/mouse.cls", "./chip_files/Mouse_Gene_Symbol_Remapping_MSigDB.v2024.1.Mm.chip", "./pathways/mouse/GOBP_INTERFERON_MEDIATED_SIGNALING_PATHWAY.v2024.1.Mm.gmt", "./results/mouse_INTERFERON_MEDIATED_SIGNALING_PATHWAY", "mouse",
  "human", "./exp_files/human_fpm.gct", "./cls_files/human.cls", "./chip_files/Human_Gene_Symbol_with_Remapping_MSigDB.v2024.1.Hs.chip", "./pathways/human/HALLMARK_INFLAMMATORY_RESPONSE.v2024.1.Hs.gmt", "./results/human_INFLAMMATORY_RESPONSE", "human",
  "human", "./exp_files/human_fpm.gct", "./cls_files/human.cls", "./chip_files/Human_Gene_Symbol_with_Remapping_MSigDB.v2024.1.Hs.chip", "./pathways/human/HALLMARK_INTERFERON_GAMMA_RESPONSE.v2024.1.Hs.gmt", "./results/human_INTERFERON_GAMMA_RESPONSE", "human",
  "human", "./exp_files/human_fpm.gct", "./cls_files/human.cls", "./chip_files/Human_Gene_Symbol_with_Remapping_MSigDB.v2024.1.Hs.chip", "./pathways/human/REACTOME_INTERFERON_GAMMA_SIGNALING.v2024.1.Hs.gmt", "./results/human_INTERFERON_GAMMA_SIGNALING", "human",
  "human", "./exp_files/human_orthologes.gct", "./cls_files/human_orthologes.cls", "./chip_files/Human_Gene_Symbol_with_Remapping_MSigDB.v2024.1.Hs.chip", "./pathways/orthologes/ZAMORA_NOS2_TARGETS_UP.v2024.1.Hs.gmt", "./results/human_orthologes_NOS2_TARGETS", "human_orthologes"
)

# Run GSEA
walk2(gsea_runs$doc_string, seq_len(nrow(gsea_runs)), function(doc, i) {
  row <- gsea_runs[i, ]
  message("Running GSEA for: ", doc)
  run_gsea(
    input_gct   = row$input_gct,
    input_cls   = row$input_cls,
    input_chip  = row$input_chip,
    gs_db       = row$gs_db,
    output_dir  = row$output_dir,
    doc_string  = doc
  )
})

# Export session info to a text file
sink("GSEA_session_info.txt")
sessionInfo()
sink()