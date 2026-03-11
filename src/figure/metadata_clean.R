#!/usr/bin/Rscript

###MRI, jsy
###Final Edit: 2026-03-07
###Description: Clean metadata

library(data.table)
library(tidyverse)
#library(tibble)

###Load A Metadata
metadata.df <- fread("/Users/jsy/Desktop/project/16s_metagenomics/metadata/raw/metadata.csv", check.names=FALSE)
metadata.df[, sample_id := paste0(sample_id, "_Pl", Platte)]
metadata.df[, sample_id := gsub("repBC_", "repBC", sample_id)]
metadata.df$Platte <-NULL
#metadata.df <- column_to_rownames(as.data.frame(metadata), var = "sample_id")

###Clean metadata labels and set factor levels
metadata.df$medium <- factor(
  recode(metadata.df$medium,
    "BHI mod." = "mBHI"
  ),
  levels = c("mBHI","GAM","Schaedler")
)

metadata.df$treatment <- factor(
  recode(metadata.df$treatment,
    "Chlorogensaeure" = "CGA",
    "ddWasser (Kontrolle Chloro.)" = "WaterControl",
    "DMSO (Kontrolle Quer.)" = "DMSOControl",
    "nur Faezes im Medium" = "FecesOnly"
  ),
  levels = c("FecesOnly","WaterControl","CGA","DMSOControl","Quercetin")
)

metadata.df$stirred <- factor(
  recode(metadata.df$stirred,
    "ungeruehrt" = "Static",
    "geruehrt" = "Stirred"
  ),
  levels = c("Static","Stirred")
)

metadata.df$time <- factor(
  recode(metadata.df$time,
    "00h" = "0",
    "24h" = "24",
    "48h" = "48"
  ),
  levels = c("0","24","48")
)

metadata.df$replicate <- factor( rep(1:3, length.out = nrow(metadata.df)), levels 
= 1:3, labels = c("R1", "R2", "R3") )

metadata.df$condition <- with(
  metadata.df,
  paste(medium, treatment, stirred, time, sep = "_")
)

metadata.df <- metadata.df %>%
  dplyr::select(sample_id, medium, treatment, stirred, time, replicate, condition)

###Save RDS files
saveRDS(metadata.df, "/Users/jsy/Desktop/project/16s_metagenomics/metadata/metadata.df.rds")
