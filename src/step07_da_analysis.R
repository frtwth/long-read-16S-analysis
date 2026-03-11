#!/usr/bin/env Rscript

##### Code by S. Jang
##### Last updated: 2026-03-09
##### Code definition: Differential abundance analysis

library(DESeq2)
library(dplyr)
library(readr)

###Step1.Define directories
dir_in <- "/Users/celina/Desktop/project/16s_metagenomics/06.DA_analysis/input_split"
dir_out <- "/Users/celina/Desktop/project/16s_metagenomics/06.DA_analysis/results"

dir.create(dir_out, recursive = TRUE, showWarnings = FALSE)

###Step2.List comparison folders
comp_dirs <- list.dirs(dir_in, recursive = FALSE, full.names = TRUE)

###Step3.Initialize summary lists
summary_list <- list()
skip_list <- list()

###Step4.Run DESeq2 for each comparison
for (dir_comp in comp_dirs) {

  comp_name <- basename(dir_comp)

  message("Running DESeq2 for: ", comp_name)

  ###Step4-1.Load data
  counts_path <- file.path(dir_comp, "counts.tsv")
  meta_path   <- file.path(dir_comp, "metadata.tsv")

  counts_df <- read_tsv(counts_path, show_col_types = FALSE)
  meta_df   <- read_tsv(meta_path, show_col_types = FALSE)

  ###Step4-2.Prepare count matrix
  counts_df <- as.data.frame(counts_df)
  rownames(counts_df) <- counts_df$genus
  counts_df$genus <- NULL

  counts_mat <- as.matrix(counts_df)

  ###Step4-3.Round estimated counts (EMU output)
  counts_mat <- round(counts_mat)

  ###Step4-4.Prepare metadata
  meta_df <- as.data.frame(meta_df)
  meta_df$condition <- factor(meta_df$condition, levels = c("ctrl", "case"))
  rownames(meta_df) <- meta_df$sample_id

  ###Step4-5.Match sample order
  meta_df <- meta_df[colnames(counts_mat), , drop = FALSE]
  stopifnot(all(colnames(counts_mat) == rownames(meta_df)))

  ###Step4-6.Create DESeq2 object
  dds <- DESeqDataSetFromMatrix(
    countData = counts_mat,
    colData   = meta_df,
    design    = ~ condition
  )

  ###Step4-7.Filter low abundance taxa
  keep <- rowSums(counts(dds)) > 10

  if (sum(keep) < 2) {
    message("Skipping ", comp_name, ": too few taxa after filtering.")

    skip_list[[length(skip_list) + 1]] <- data.frame(
      comparison = comp_name,
      reason = "Too few taxa after filtering"
    )

    next
  }

  dds <- dds[keep, ]

  ###Step4-8.Run DESeq2 with default size-factor normalization
  dds <- DESeq(dds)

  res <- results(dds)

  ###Step4-9.Convert result table safely
  res_df <- as.data.frame(res)
  res_df$genus <- rownames(res_df)
  res_df <- res_df %>%
    arrange(is.na(padj), padj)

  ###Step4-10.Get normalized counts
  norm_counts <- counts(dds, normalized = TRUE) %>%
    as.data.frame()

  norm_counts$genus <- rownames(norm_counts)

  ###Step4-11.Save size factors
  sf <- sizeFactors(dds)
  size_factors <- data.frame(
    sample = names(sf),
    size_factor = as.numeric(sf)
  )

  ###Step4-12.Output directory
  dir_comp_out <- file.path(dir_out, comp_name)
  dir.create(dir_comp_out, recursive = TRUE, showWarnings = FALSE)

  ###Step4-13.Save full results
  write_csv(
    res_df,
    file.path(dir_comp_out, "deseq2_results.csv")
  )

  ###Step4-14.Save significant taxa (padj < 0.05)
  sig_df_005 <- res_df %>%
    filter(!is.na(padj), padj < 0.05)

  write_csv(
    sig_df_005,
    file.path(dir_comp_out, "sig_padj0.05.csv")
  )

  ###Step4-15.Save significant taxa (padj < 0.01)
  sig_df_001 <- res_df %>%
    filter(!is.na(padj), padj < 0.01)

  write_csv(
    sig_df_001,
    file.path(dir_comp_out, "sig_padj0.01.csv")
  )

  ###Step4-16.Save normalized counts
  write_tsv(
    norm_counts,
    file.path(dir_comp_out, "normalized_counts.tsv")
  )

  ###Step4-17.Save size factors
  write_csv(
    size_factors,
    file.path(dir_comp_out, "size_factors.csv")
  )

  ###Step4-18.Store summary info
  summary_list[[comp_name]] <- data.frame(
    comparison = comp_name,
    normalization = "DESeq2_default",
    n_taxa_tested = nrow(res_df),
    sig_taxa_padj0.05 = sum(!is.na(res_df$padj) & res_df$padj < 0.05),
    sig_taxa_padj0.01 = sum(!is.na(res_df$padj) & res_df$padj < 0.01),
    up_padj0.05 = sum(!is.na(res_df$padj) & res_df$padj < 0.05 & res_df$log2FoldChange > 0),
    down_padj0.05 = sum(!is.na(res_df$padj) & res_df$padj < 0.05 & res_df$log2FoldChange < 0),
    min_size_factor = min(sf),
    median_size_factor = median(sf),
    max_size_factor = max(sf)
  )
}

###Step5.Make summary tables
if (length(summary_list) > 0) {
  summary_table <- bind_rows(summary_list) %>%
    arrange(desc(sig_taxa_padj0.05))
} else {
  summary_table <- data.frame(
    comparison = character(),
    normalization = character(),
    n_taxa_tested = numeric(),
    sig_taxa_padj0.05 = numeric(),
    sig_taxa_padj0.01 = numeric(),
    up_padj0.05 = numeric(),
    down_padj0.05 = numeric(),
    min_size_factor = numeric(),
    median_size_factor = numeric(),
    max_size_factor = numeric()
  )
}

write_csv(
  summary_table,
  file.path(dir_out, "DA_summary_table.csv")
)

###Step6.Save skipped comparisons
if (length(skip_list) > 0) {
  skip_table <- bind_rows(skip_list)
} else {
  skip_table <- data.frame(
    comparison = character(),
    reason = character()
  )
}

write_csv(
  skip_table,
  file.path(dir_out, "DA_skipped_comparisons.csv")
)
