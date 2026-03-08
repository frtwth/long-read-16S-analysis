#!/usr/bin/env Rscript

##### Code by S. Jang
##### Last updated: 2026-03-06
##### Code definition: Split EMU lineage column into taxonomy ranks

library(readr)
library(tidyr)
library(yaml)

args <- commandArgs(trailingOnly = TRUE)
config_file <- args[1]

config <- yaml::read_yaml(config_file)

root <- config$project$root
dir_in <- file.path(root, config$taxonomy_split$dir_in_rel)
dir_out <- file.path(root, config$taxonomy_split$dir_out_rel)

dir.create(dir_out, showWarnings = FALSE, recursive = TRUE)

files <- list.files(
  path = dir_in,
  pattern = "_rel-abundance.tsv$",
  recursive = TRUE,
  full.names = TRUE
)

if (length(files) == 0) {
  stop(paste("No *_rel-abundance.tsv files found in:", dir_in))
}

for (f in files) {
  df <- read_tsv(f, show_col_types = FALSE)

  df_split <- separate(
    data = df,
    col = lineage,
    into = c("superkingdom", "phylum", "class", "order", "family", "genus", "species"),
    sep = ";",
    fill = "right",
    remove = FALSE
  )
  
  df_split[df_split == ""] <- NA
  df_split[df_split == "NA"] <- NA

  sample_name <- basename(f)
  sample_name <- sub("_rel-abundance.tsv$", "", sample_name)

  out_file <- file.path(dir_out, paste0(sample_name, "_split_rel-abundance.tsv"))
  write_tsv(df_split, out_file)
}
