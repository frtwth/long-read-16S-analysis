#!/usr/bin/env Rscript

##### Code by S. Jang
##### Last updated: 2026-03-10
##### Code definition: Prepare data for heatmap

library(dplyr)
library(readr)
library(stringr)
library(tidyr)

###Step1.Define directories
dir_in <- "/Users/celina/Desktop/project/16s_metagenomics/06.DA_analysis/results"
dir_out <- "/Users/celina/Desktop/project/16s_metagenomics/06.DA_analysis/heatmap_data"

dir.create(dir_out, recursive = TRUE, showWarnings = FALSE)

###Step2.List comparison result folders
comp_dirs <- list.dirs(dir_in, recursive = FALSE, full.names = TRUE)

###Step3.Initialize result lists
res_list <- list()
skip_list <- list()

###Step4.Read DESeq2 result table from each comparison folder
for (dir_comp in comp_dirs) {

  comp_name <- basename(dir_comp)

  message("Preparing heatmap data for: ", comp_name)

  ###Step4-1.Set result file path
  res_path <- file.path(dir_comp, "deseq2_results.csv")

  ###Step4-2.Skip if result file does not exist
  if (!file.exists(res_path)) {
    message("Skipping ", comp_name, ": deseq2_results.csv not found.")

    skip_list[[length(skip_list) + 1]] <- data.frame(
      comparison = comp_name,
      reason = "deseq2_results.csv not found"
    )

    next
  }

  ###Step4-3.Read result table
  res_df <- read_csv(res_path, show_col_types = FALSE)

  ###Step4-4.Check required columns
  required_cols <- c("genus", "log2FoldChange", "stat", "padj")
  missing_cols <- setdiff(required_cols, colnames(res_df))

  if (length(missing_cols) > 0) {
    message("Skipping ", comp_name, ": missing required columns.")

    skip_list[[length(skip_list) + 1]] <- data.frame(
      comparison = comp_name,
      reason = paste("Missing columns:", paste(missing_cols, collapse = ", "))
    )

    next
  }

  ###Step4-5.Parse comparison name
  comp_split <- str_split(comp_name, "_", simplify = TRUE)

  if (ncol(comp_split) != 4) {
    message("Skipping ", comp_name, ": invalid comparison name format.")

    skip_list[[length(skip_list) + 1]] <- data.frame(
      comparison = comp_name,
      reason = "Invalid comparison name format"
    )

    next
  }

  medium <- comp_split[1]
  time <- comp_split[2]
  agitation <- comp_split[3]
  compound <- comp_split[4]

  ###Step4-6.Add comparison information
  res_df <- res_df %>%
    mutate(
      comparison = comp_name,
      medium = medium,
      time = time,
      agitation = agitation,
      compound = compound
    ) %>%
    select(
      comparison,
      medium,
      time,
      agitation,
      compound,
      genus,
      log2FoldChange,
      stat,
      padj
    )

  res_list[[comp_name]] <- res_df
}

###Step5.Combine all result tables
if (length(res_list) > 0) {
  heatmap_all <- bind_rows(res_list)
} else {
  heatmap_all <- data.frame(
    comparison = character(),
    medium = character(),
    time = character(),
    agitation = character(),
    compound = character(),
    genus = character(),
    log2FoldChange = numeric(),
    stat = numeric(),
    padj = numeric()
  )
}

###Step6.Clean and order variables
heatmap_all <- heatmap_all %>%
  mutate(
    time = str_remove(time, "^t"),
    time = factor(time, levels = c("0", "24", "48")),
    medium = factor(medium, levels = c("GAM", "mBHI", "Schaedler")),
    agitation = factor(agitation, levels = c("Static", "Stirred")),
    compound = factor(compound, levels = c("Quercetin", "CGA"))
  )

###Step7.Save combined heatmap data
write_csv(
  heatmap_all,
  file.path(dir_out, "heatmap_data_all.csv")
)

###Step8.Select genera for Quercetin heatmap
selected_genera_quercetin <- heatmap_all %>%
  filter(
    compound == "Quercetin",
    !is.na(padj),
    padj < 0.05
  ) %>%
  distinct(genus) %>%
  arrange(genus)

###Step9.Make Quercetin heatmap dataset
heatmap_quercetin <- heatmap_all %>%
  filter(
    compound == "Quercetin",
    genus %in% selected_genera_quercetin$genus
  ) %>%
  mutate(
    sig_padj0.05 = !is.na(padj) & padj < 0.05,
    sig_padj0.01 = !is.na(padj) & padj < 0.01
  )

###Step10.Order Quercetin genera
genus_order_quercetin <- heatmap_quercetin %>%
  group_by(genus) %>%
  summarise(
    n_sig_padj0.05 = sum(sig_padj0.05, na.rm = TRUE),
    min_padj = if (all(is.na(padj))) NA_real_ else min(padj, na.rm = TRUE),
    max_abs_stat = if (all(is.na(stat))) NA_real_ else max(abs(stat), na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(n_sig_padj0.05), min_padj, desc(max_abs_stat), genus)

heatmap_quercetin <- heatmap_quercetin %>%
  mutate(
    genus = factor(genus, levels = rev(genus_order_quercetin$genus))
  )

###Step11.Save Quercetin outputs
write_csv(
  selected_genera_quercetin,
  file.path(dir_out, "selected_genera_quercetin.csv")
)

write_csv(
  heatmap_quercetin,
  file.path(dir_out, "heatmap_data_quercetin.csv")
)

write_csv(
  genus_order_quercetin,
  file.path(dir_out, "genus_order_quercetin.csv")
)

###Step12.Select genera for CGA heatmap
selected_genera_cga <- heatmap_all %>%
  filter(
    compound == "CGA",
    !is.na(padj),
    padj < 0.05
  ) %>%
  distinct(genus) %>%
  arrange(genus)

###Step13.Make CGA heatmap dataset
heatmap_cga <- heatmap_all %>%
  filter(
    compound == "CGA",
    genus %in% selected_genera_cga$genus
  ) %>%
  mutate(
    sig_padj0.05 = !is.na(padj) & padj < 0.05,
    sig_padj0.01 = !is.na(padj) & padj < 0.01
  )

###Step14.Order CGA genera
genus_order_cga <- heatmap_cga %>%
  group_by(genus) %>%
  summarise(
    n_sig_padj0.05 = sum(sig_padj0.05, na.rm = TRUE),
    min_padj = if (all(is.na(padj))) NA_real_ else min(padj, na.rm = TRUE),
    max_abs_stat = if (all(is.na(stat))) NA_real_ else max(abs(stat), na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(n_sig_padj0.05), min_padj, desc(max_abs_stat), genus)

heatmap_cga <- heatmap_cga %>%
  mutate(
    genus = factor(genus, levels = rev(genus_order_cga$genus))
  )

###Step15.Save CGA outputs
write_csv(
  selected_genera_cga,
  file.path(dir_out, "selected_genera_cga.csv")
)

write_csv(
  heatmap_cga,
  file.path(dir_out, "heatmap_data_cga.csv")
)

write_csv(
  genus_order_cga,
  file.path(dir_out, "genus_order_cga.csv")
)

###Step16.Save skipped comparisons
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
  file.path(dir_out, "heatmap_prepare_skipped.csv")
)

###Step17.Print summary
message("Heatmap data preparation completed.")
message("Total rows in combined heatmap data: ", nrow(heatmap_all))
message("Number of selected genera for Quercetin: ", nrow(selected_genera_quercetin))
message("Number of selected genera for CGA: ", nrow(selected_genera_cga))

