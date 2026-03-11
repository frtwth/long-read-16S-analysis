#!/usr/bin/env Rscript

##### Code by S. Jang
##### Last updated: 2026-03-09
##### Code definition: Table split for differential abundance analysis

library(data.table)
library(dplyr)
library(readr)

###Step1.Load data
genus_cts <- fread("/Users/jsy/Desktop/project/16s_metagenomics/04.split_taxonomy/emu-combined-genus-counts.tsv")
metadata_raw <- readRDS("/Users/jsy/Desktop/project/16s_metagenomics/metadata/metadata.df.rds")

dir_out <- "/Users/jsy/Desktop/project/16s_metagenomics/06.DA_analysis/input_split"
dir.create(dir_out, recursive = TRUE, showWarnings = FALSE)

###Step2.Preprocess genus counts table
genus_cts <- genus_cts %>%
  select(-family, -order, -class, -phylum, -superkingdom) %>%
  filter(!is.na(genus), genus != "")

colnames(genus_cts) <- gsub(
  "_trimmed_filtered_split$",
  "",
  colnames(genus_cts)
)

genus_cts[is.na(genus_cts)] <- 0

genus_cts <- as.data.frame(genus_cts)

###Step3.Check sample IDs
if (!all(metadata_raw$sample_id %in% colnames(genus_cts))) {
  missing_samples <- setdiff(metadata_raw$sample_id, colnames(genus_cts))
  stop(
    "Some sample IDs in metadata are missing from genus count table:\n",
    paste(missing_samples, collapse = ", ")
  )
}

###Step4.Define DA comparison pairs
comparison_list <- list(
  CGA = c("WaterControl", "CGA"),
  Quercetin = c("DMSOControl", "Quercetin")
)

###Step5.Get unique medium-time-stirred combinations
comb_df <- metadata_raw %>%
  distinct(medium, time, stirred) %>%
  arrange(medium, time, stirred)

###Step6.Initialize summary tables
manifest_list <- list()
skip_list <- list()

###Step7.Generate DA input folders
for (i in seq_len(nrow(comb_df))) {

  current_medium  <- comb_df$medium[i]
  current_time    <- comb_df$time[i]
  current_stirred <- comb_df$stirred[i]

  message(
    "Processing: ",
    current_medium, " / ",
    current_time, " / ",
    current_stirred
  )

  meta_base <- metadata_raw %>%
    filter(
      medium == current_medium,
      time == current_time,
      stirred == current_stirred
    )

  for (comp_name in names(comparison_list)) {

    ctrl_case <- comparison_list[[comp_name]]
    ctrl_name <- ctrl_case[1]
    case_name <- ctrl_case[2]

    meta_sub <- meta_base %>%
      filter(treatment %in% c(ctrl_name, case_name))

    ###Step7-1.Skip if either group is missing
    group_counts <- table(meta_sub$treatment)

    if (!(ctrl_name %in% names(group_counts)) || !(case_name %in% names(group_counts))) {
      skip_list[[length(skip_list) + 1]] <- data.frame(
        comparison_id = paste(current_medium, paste0("t", current_time), current_stirred, comp_name, sep = "_"),
        medium = current_medium,
        time = current_time,
        stirred = current_stirred,
        comparison = comp_name,
        reason = "One or both treatment groups missing"
      )
      next
    }

    ###Step7-2.Add DESeq2 condition column
    meta_sub <- meta_sub %>%
      mutate(
        condition = case_when(
          treatment == ctrl_name ~ "ctrl",
          treatment == case_name ~ "case"
        )
      )

    ###Step7-3.Reorder samples by condition
    meta_sub <- meta_sub %>%
      arrange(condition, replicate)

    sample_vec <- meta_sub$sample_id

    ###Step7-4.Subset count table
    counts_sub <- genus_cts %>%
      select(genus, all_of(sample_vec))

    ###Step7-5.Make output folder
    comparison_id <- paste(
      current_medium,
      paste0("t", current_time),
      current_stirred,
      comp_name,
      sep = "_"
    )

    dir_comp <- file.path(dir_out, comparison_id)
    dir.create(dir_comp, recursive = TRUE, showWarnings = FALSE)

    ###Step7-6.Save metadata
    meta_save <- meta_sub %>%
      select(sample_id, treatment, condition, medium, time, stirred, replicate)

    write_tsv(
      meta_save,
      file.path(dir_comp, "metadata.tsv")
    )

    ###Step7-7.Save counts
    write_tsv(
      counts_sub,
      file.path(dir_comp, "counts.tsv")
    )

    ###Step7-8.Save manifest info
    manifest_list[[length(manifest_list) + 1]] <- data.frame(
      comparison_id = comparison_id,
      medium = current_medium,
      time = current_time,
      stirred = current_stirred,
      comparison = comp_name,
      control = ctrl_name,
      treatment = case_name,
      n_ctrl = sum(meta_sub$condition == "ctrl"),
      n_case = sum(meta_sub$condition == "case"),
      n_total = nrow(meta_sub)
    )
  }
}

###Step8.Save summary files
dir_summary <- "/Users/jsy/Desktop/project/16s_metagenomics/06.DA_analysis/summary"
dir.create(dir_summary, recursive = TRUE, showWarnings = FALSE)

manifest_df <- bind_rows(manifest_list)
write_csv(
  manifest_df,
  file.path(dir_summary, "comparison_manifest.csv")
)

if (length(skip_list) > 0) {
  skip_df <- bind_rows(skip_list)
} else {
  skip_df <- data.frame(
    comparison_id = character(),
    medium = character(),
    time = character(),
    stirred = character(),
    comparison = character(),
    reason = character()
  )
}

write_csv(
  skip_df,
  file.path(dir_summary, "skipped_comparisons.csv")
)
