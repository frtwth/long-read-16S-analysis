#!/usr/bin/env Rscript

##### Code by S. Jang
##### Last updated: 2026-03-08
##### Code definition: Beta diversity analysis

library(data.table)
library(dplyr)
library(vegan)
library(ggplot2)

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

###Step6.Prepare metadata
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

###Step7.Split analysis by medium
medium_levels <- levels(metadata$medium)

beta_by_medium <- vector("list", length(medium_levels))
names(beta_by_medium) <- medium_levels

for (current_medium in medium_levels) {

  message("Running beta-diversity analysis for medium: ", current_medium)

  ###Step7-1.Subset metadata
  metadata_sub <- metadata %>%
    filter(medium == current_medium)

  sample_ids <- rownames(metadata_sub)

  ###Step7-2.Subset abundance matrix
  abundance_mat_sub <- abundance_mat[, sample_ids, drop = FALSE]

  stopifnot(all(colnames(abundance_mat_sub) == rownames(metadata_sub)))

  ###Step7-3.Calculate Bray-Curtis distance
  bray_dist_sub <- vegdist(t(abundance_mat_sub), method = "bray")

  ###Step7-4.PERMANOVA within medium
  permanova_res_sub <- adonis2(
    bray_dist_sub ~ time + stirred + treatment,
    data = metadata_sub,
    by = "margin"
  )

  print(permanova_res_sub)

  ###Step7-5.Homogeneity of dispersion
  bd_time_sub <- betadisper(bray_dist_sub, metadata_sub$time)
  bd_stirred_sub <- betadisper(bray_dist_sub, metadata_sub$stirred)
  bd_treatment_sub <- betadisper(bray_dist_sub, metadata_sub$treatment)

  dispersion_res_sub <- list(
    time = anova(bd_time_sub),
    stirred = anova(bd_stirred_sub),
    treatment = anova(bd_treatment_sub)
  )

  ###Step7-6.PCoA coordinates
  pcoa_res_sub <- cmdscale(bray_dist_sub, eig = TRUE, k = 2)

  pcoa_df_sub <- data.frame(
    sample_id = rownames(pcoa_res_sub$points),
    PCoA1 = pcoa_res_sub$points[, 1],
    PCoA2 = pcoa_res_sub$points[, 2],
    stringsAsFactors = FALSE
  ) %>%
    left_join(
      metadata_sub,
      by = "sample_id"
    )

  pcoa_var_sub <- round(
    100 * pcoa_res_sub$eig / sum(pcoa_res_sub$eig[pcoa_res_sub$eig > 0]),
    1
  )

  ###Step7-7.Store results
  beta_by_medium[[current_medium]] <- list(
    medium = current_medium,
    sample_ids = sample_ids,
    metadata = metadata_sub,
    abundance_mat = abundance_mat_sub,
    bray_dist = bray_dist_sub,
    permanova_res = permanova_res_sub,
    dispersion_res = dispersion_res_sub,
    pcoa_df = pcoa_df_sub,
    pcoa_var = pcoa_var_sub
  )
}

###Step8.Save outputs
saveRDS(
  beta_by_medium,
  file = "/Users/celina/Desktop/project/16s_metagenomics/05.beta_diversity/beta_diversity_by_medium_genus.rds"
)
