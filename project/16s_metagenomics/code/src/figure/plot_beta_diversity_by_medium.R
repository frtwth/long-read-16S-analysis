#!/usr/bin/Rscript

#### Code by S. Jang
#### Last updated: 2026-03-08
#### Code definition: Plot Beta diversity results

library(dplyr)
library(ggplot2)
library(tibble)

source("/Users/celina/Desktop/project/16s_metagenomics/scripts/plot_theme.R")

### Step 1. Set paths
rds_file <- "/Users/celina/Desktop/project/16s_metagenomics/05.beta_diversity/beta_diversity_by_medium_genus.rds"
fig_dir  <- "/Users/celina/Desktop/project/16s_metagenomics/results/figures/beta_diversity_by_medium"
tab_dir  <- "/Users/celina/Desktop/project/16s_metagenomics/results/tables/beta_diversity_by_medium"

dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(tab_dir, showWarnings = FALSE, recursive = TRUE)

### Step 2. Load results
beta_res_by_medium <- readRDS(rds_file)
medium_names <- names(beta_res_by_medium)

### Step 3. Define treatment colors
treatment_cols <- c(
  "FecesOnly"    = "#4C78A8",
  "WaterControl" = "#9C755F",
  "CGA"          = "#F58518",
  "DMSOControl"  = "#B8B8B8",
  "Quercetin"    = "#54A24B"
)

### Step 4. Initialize output lists
pcoa_plot_list  <- list()
perma_table_list <- list()
perma_plot_list <- list()

### Step 5. Generate plots for each medium
for (current_medium in medium_names) {

  message("Plotting beta-diversity results for medium: ", current_medium)

  current_res   <- beta_res_by_medium[[current_medium]]
  pcoa_df       <- current_res$pcoa_df
  pcoa_var      <- current_res$pcoa_var
  permanova_res <- current_res$permanova_res

  ### Step 5-1. Convert variables to factors
  pcoa_df <- pcoa_df %>%
    mutate(
      medium = factor(medium),
      time = factor(time, levels = c("0", "24", "48")),
      stirred = factor(stirred),
      treatment = factor(treatment)
    )

  ### Step 5-2. Adjust treatment colors if needed
  treatment_levels <- levels(pcoa_df$treatment)

  if (!all(treatment_levels %in% names(treatment_cols))) {
    missing_levels <- setdiff(treatment_levels, names(treatment_cols))
    extra_cols <- setNames(rep("#666666", length(missing_levels)), missing_levels)
    treatment_cols_use <- c(treatment_cols, extra_cols)
  } else {
    treatment_cols_use <- treatment_cols
  }

  ### Step 5-3. Prepare PERMANOVA table
  perma_df <- as.data.frame(permanova_res) %>%
    rownames_to_column("Factor")

  p_col <- grep("Pr", colnames(perma_df), value = TRUE)

  perma_table <- perma_df %>%
    mutate(
      SumOfSqs = sprintf("%.3f", SumOfSqs),
      R2 = sprintf("%.4f", R2),
      F = ifelse(is.na(F), "-", sprintf("%.2f", F)),
      p = ifelse(
        is.na(.data[[p_col]]),
        "-",
        ifelse(.data[[p_col]] <= 0.001, "<0.001", sprintf("%.3f", .data[[p_col]]))
      )
    ) %>%
    select(Factor, Df, SumOfSqs, R2, F, p)

  perma_table_list[[current_medium]] <- perma_table
  print(perma_table)

  write.csv(
    perma_table,
    file.path(tab_dir, paste0("permanova_table_", current_medium, "_genus.csv")),
    row.names = FALSE
  )

  ### Step 5-4. PCoA plot
  p_pcoa_medium <- ggplot(
    pcoa_df,
    aes(PCoA1, PCoA2, color = treatment, shape = time)
  ) +
    geom_point(size = 2.8, alpha = 0.9) +
    scale_color_manual(
      values = treatment_cols_use,
      name = "Treatment",
      breaks = treatment_levels
    ) +
    scale_shape_manual(
      values = time_shapes,
      name = "Time [hrs]"
    ) +
    labs(
      title = current_medium,
      x = paste0("PCoA1 (", pcoa_var[1], "%)"),
      y = paste0("PCoA2 (", pcoa_var[2], "%)")
    ) +
    theme_16S()

  pcoa_plot_list[[current_medium]] <- p_pcoa_medium
  print(p_pcoa_medium)

  ggsave(
    filename = file.path(fig_dir, paste0("pcoa_", current_medium, "_genus.pdf")),
    plot = p_pcoa_medium,
    width = 6.5,
    height = 5.0
  )

  ### Step 5-5. PERMANOVA bar plot
  perma_bar_df <- perma_df %>%
    filter(Factor %in% c("time", "treatment", "stirred")) %>%
    mutate(
      Factor = factor(Factor, levels = c("time", "treatment", "stirred")),
      p_sig = case_when(
        .data[[p_col]] <= 0.001 ~ "***",
        .data[[p_col]] <= 0.01  ~ "**",
        .data[[p_col]] <= 0.05  ~ "*",
        .data[[p_col]] <= 0.1   ~ ".",
        TRUE ~ ""
      )
    )

  factor_cols_sub <- factor_cols[c("time", "treatment", "stirred")]
  y_max <- max(perma_bar_df$R2, na.rm = TRUE)

  p_permanova_medium <- ggplot(
    perma_bar_df,
    aes(x = Factor, y = R2, fill = Factor)
  ) +
    geom_col(width = 0.72, color = "black", linewidth = 0.4) +
    geom_text(
      aes(label = sprintf("%.3f", R2)),
      vjust = -0.35,
      size = 3.6
    ) +
    geom_text(
      aes(label = p_sig),
      vjust = -1.2,
      size = 5.5
    ) +
    scale_fill_manual(
      values = factor_cols_sub,
      breaks = c("time", "treatment", "stirred"),
      drop = FALSE,
      name = "Factor"
    ) +
    coord_cartesian(ylim = c(0, y_max * 1.20)) +
    labs(
      title = current_medium,
      x = NULL,
      y = expression(R^2)
    ) +
    theme_16S() +
    theme(legend.position = "none")

  perma_plot_list[[current_medium]] <- p_permanova_medium
  print(p_permanova_medium)

  ggsave(
    filename = file.path(fig_dir, paste0("permanova_", current_medium, "_genus.pdf")),
    plot = p_permanova_medium,
    width = 5.5,
    height = 4.5
  )
}
