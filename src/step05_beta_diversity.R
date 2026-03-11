#!/usr/bin/env Rscript

##### Code by S. Jang
##### Last updated: 2026-03-08
##### Code definition: Beta diversity analysis

library(data.table)
library(dplyr)
library(vegan)

###Step1.Load data
genus_abundance <- fread("/Users/celina/Desktop/project/16s_metagenomics/04.split_taxonomy/emu-combined-genus.tsv")
metadata_raw <- readRDS("/Users/celina/Desktop/project/16s_metagenomics/metadata/metadata.df.rds")

###Step2.Preprocess genus abundance table
genus_abundance <- genus_abundance %>%
  select(-family, -order, -class, -phylum, -superkingdom) %>%
  filter(!is.na(genus), genus != "")

colnames(genus_abundance) <- gsub(
  "_trimmed_filtered_split$",
  "",
  colnames(genus_abundance)
)

genus_abundance[is.na(genus_abundance)] <- 0
genus_abundance <- as.data.frame(genus_abundance)

###Step3.Filter rare taxa
abundance_mat <- genus_abundance[, -1, drop = FALSE]
abundance_mat[abundance_mat < 0.0001] <- 0

genus_abundance[, -1] <- abundance_mat
genus_abundance <- genus_abundance[
  rowSums(genus_abundance[, -1, drop = FALSE]) > 0,
]

###Step4.Re-normalize within each sample
abundance_mat <- genus_abundance[, -1, drop = FALSE]
sample_sums <- colSums(abundance_mat)

stopifnot(all(sample_sums > 0))

genus_abundance[, -1] <- sweep(
  abundance_mat,
  2,
  sample_sums,
  "/"
)

###Step5.Prepare abundance matrix for beta diversity
abundance_df <- genus_abundance
rownames(abundance_df) <- abundance_df$genus
abundance_df$genus <- NULL

abundance_mat <- as.matrix(abundance_df)
storage.mode(abundance_mat) <- "numeric"

### Step 6. Prepare metadata
metadata <- as.data.frame(metadata_raw)
rownames(metadata) <- metadata$sample_id
metadata <- metadata[colnames(abundance_mat), , drop = FALSE]

metadata <- metadata %>%
  mutate(
    medium = factor(medium),
    time = factor(time),
    stirred = factor(stirred),
    treatment = factor(treatment)
  )

stopifnot(all(colnames(abundance_mat) == rownames(metadata)))

###Step7.Calculate Bray-Curtis distance
bray_dist <- vegdist(t(abundance_mat), method = "bray")

###Step8.PERMANOVA
permanova_res <- adonis2(
  bray_dist ~ medium + time + stirred + treatment,
  data = metadata,
  by = "margin"
)

print(permanova_res)

###Step9.Homogeneity of dispersion
bd_medium <- betadisper(bray_dist, metadata$medium)
bd_time <- betadisper(bray_dist, metadata$time)
bd_stirred <- betadisper(bray_dist, metadata$stirred)
bd_treatment <- betadisper(bray_dist, metadata$treatment)

dispersion_res <- list(
  medium = anova(bd_medium),
  time = anova(bd_time),
  stirred = anova(bd_stirred),
  treatment = anova(bd_treatment)
)

###Step10.PCoA coordinates
pcoa_res <- cmdscale(bray_dist, eig = TRUE, k = 2)

pcoa_df <- data.frame(
  sample_id = rownames(pcoa_res$points),
  PCoA1 = pcoa_res$points[, 1],
  PCoA2 = pcoa_res$points[, 2],
  stringsAsFactors = FALSE
) %>%
  left_join(
    metadata,
    by = "sample_id"
  )

pcoa_var <- round(100 * pcoa_res$eig / sum(pcoa_res$eig[pcoa_res$eig > 0]), 1)

###Step11.Save outputs
saveRDS(
  list(
    abundance_mat = abundance_mat,
    metadata = metadata,
    bray_dist = bray_dist,
    permanova_res = permanova_res,
    dispersion_res = dispersion_res,
    pcoa_df = pcoa_df,
    pcoa_var = pcoa_var
  ),
  file = "/Users/celina/Desktop/project/16s_metagenomics/05.beta_diversity/beta_diversity_genus.rds"
)
