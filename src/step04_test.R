#!/usr/bin/env Rscript

library(readr)
library(dplyr)
library(tidyr)
library(stringr)
library(yaml)

args <- commandArgs(trailingOnly = TRUE)
cfg <- yaml::read_yaml(args[1])

dir_in  <- file.path(cfg$project$root, cfg$taxonomy_split$dir_in_rel)
dir_out <- file.path(cfg$project$root, cfg$taxonomy_split$dir_out_rel)

dir.create(dir_out, showWarnings = FALSE, recursive = TRUE)

tax_cols <- c("superkingdom","phylum","class","order","family","genus","species")
files <- list.files(dir_in, pattern = "_rel-abundance.tsv$", full.names = TRUE)

for (f in files) {
  df <- read_tsv(f, show_col_types = FALSE) %>%
    mutate(lineage = str_trim(lineage)) %>%
    filter(!is.na(lineage), lineage != "") %>%
    separate(lineage, into = tax_cols, sep = ";", fill = "right", extra = "drop", remove = FALSE) %>%
    mutate(across(all_of(tax_cols), ~ na_if(str_trim(.x), "")))

  out_file <- file.path(
    dir_out,
    paste0(sub("_rel-abundance.tsv$", "", basename(f)), "_split_rel-abundance.tsv")
  )

  write_tsv(df, out_file)
}

cat("\n===== warnings() =====\n")
print(warnings())
