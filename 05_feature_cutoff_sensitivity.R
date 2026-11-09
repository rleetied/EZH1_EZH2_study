suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
})

root <- "RKO_TAZ_v_VAL_project/Final_figure_analysis/"
outdir <- file.path(root, "outputs", "FIG2_FEATURE_CUTOFF_SENSITIVITY_1K_10K")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

paths <- c(
  RKO = "RKO_TAZ_v_VAL_project/RNA-seq/RKO/RKO.DAC_TAZ_VAL_treatments_D3.ann.stats.0.25.txt",
  SW620 = "RKO_TAZ_v_VAL_project/RNA-seq/SW620/SW620.DAC_TAZ_VAL_treatments_D3.ann.stats.0.25.txt"
)
treatments <- c("Veh", "DAC", "TAZ", "VAL", "TD", "VD")
sample_cols <- unlist(lapply(treatments, function(z) paste0(z, 1:3)))

read_one <- function(path) {
  fread(path)[gene_biotype == "protein_coding" & !is.na(ENSG)]
}
rko <- read_one(paths[["RKO"]])
sw620 <- read_one(paths[["SW620"]])
common <- intersect(rko$ENSG, sw620$ENSG)
rko <- rko[match(common, ENSG)]
sw620 <- sw620[match(common, ENSG)]
stopifnot(all(rko$ENSG == sw620$ENSG))

make_expr <- function(x, line) {
  m <- log2(as.matrix(x[, ..sample_cols]) + 0.25)
  rownames(m) <- x$ENSG
  colnames(m) <- paste(line, sample_cols, sep = "_")
  m
}
er <- make_expr(rko, "RKO")
es <- make_expr(sw620, "SW620")

# Use one common ranking based on the mean of the within-line variances. This
# prevents the higher-variance cell line from determining the feature space.
variance_rko <- apply(er, 1, var)
variance_sw620 <- apply(es, 1, var)
mean_within_line_variance <- (variance_rko + variance_sw620) / 2
ranked <- order(mean_within_line_variance, decreasing = TRUE)

ranking <- data.table(
  rank = seq_along(ranked),
  ENSG = rko$ENSG[ranked],
  gene_symbol = rko$ext_gene[ranked],
  variance_RKO = variance_rko[ranked],
  variance_SW620 = variance_sw620[ranked],
  mean_within_line_variance = mean_within_line_variance[ranked]
)
fwrite(ranking, file.path(outdir, "COMMON_PROTEIN_CODING_GENE_VARIANCE_RANKING.tsv"), sep = "\t")

calc_distances <- function(m, idx, line, pair) {
  mm <- m[idx, , drop = FALSE]
  veh <- rowMeans(mm[, grepl("_Veh[123]$", colnames(mm)), drop = FALSE])
  rbindlist(lapply(pair, function(tr) {
    jj <- grepl(paste0("_", tr, "[123]$"), colnames(mm))
    data.table(
      cell_line = line,
      treatment = tr,
      sample = colnames(mm)[jj],
      RMS_distance_from_VEH = apply(
        mm[, jj, drop = FALSE], 2,
        function(z) sqrt(mean((z - veh)^2))
      )
    )
  }))
}

analyze_pair <- function(idx, cutoff, pair, comparison) {
  d <- rbindlist(list(
    calc_distances(er, idx, "RKO", pair),
    calc_distances(es, idx, "SW620", pair)
  ))
  d[, cell_line := factor(cell_line, levels = c("RKO", "SW620"))]
  d[, treatment := factor(treatment, levels = pair)]

  fit <- lm(RMS_distance_from_VEH ~ cell_line * treatment, data = d)
  co <- coef(summary(fit))
  V <- vcov(fit)
  bb <- coef(fit)
  interaction_name <- grep(":", rownames(co), value = TRUE)
  interaction_estimate <- co[interaction_name, "Estimate"]
  interaction_se <- co[interaction_name, "Std. Error"]
  critical <- qt(0.975, df.residual(fit))

  contrast_p <- function(line) {
    cc <- if (line == "RKO") c(0, 0, 1, 0) else c(0, 0, 1, 1)
    estimate <- sum(cc * bb)
    se <- sqrt(drop(t(cc) %*% V %*% cc))
    2 * pt(-abs(estimate / se), df = df.residual(fit))
  }

  line_stats <- d[, .(
    first_mean = mean(RMS_distance_from_VEH[treatment == pair[1]]),
    second_mean = mean(RMS_distance_from_VEH[treatment == pair[2]])
  ), by = cell_line]
  line_stats[, `:=`(
    fold = second_mean / first_mean,
    within_line_p = vapply(as.character(cell_line), contrast_p, numeric(1))
  )]

  summary <- data.table(
    cutoff = cutoff,
    n_genes = length(idx),
    comparison = comparison,
    first_treatment = pair[1],
    second_treatment = pair[2],
    interaction_definition = "(second-first) SW620 minus (second-first) RKO",
    interaction_estimate = interaction_estimate,
    interaction_CI_low = interaction_estimate - critical * interaction_se,
    interaction_CI_high = interaction_estimate + critical * interaction_se,
    interaction_p = co[interaction_name, "Pr(>|t|)"],
    RKO_fold = line_stats[cell_line == "RKO", fold],
    RKO_within_line_p = line_stats[cell_line == "RKO", within_line_p],
    SW620_fold = line_stats[cell_line == "SW620", fold],
    SW620_within_line_p = line_stats[cell_line == "SW620", within_line_p]
  )
  list(summary = summary, samples = cbind(cutoff = cutoff, comparison = comparison, d))
}

cutoffs <- c(1000L, 2000L, 5000L, 10000L)
cutoff_labels <- paste0("Top ", format(cutoffs, big.mark = ","))
pairs <- list(c("TAZ", "VAL"), c("TD", "VD"))
comparison_names <- c("VAL vs TAZ", "VD vs TD")

results <- list()
k <- 1L
for (i in seq_along(cutoffs)) {
  idx <- ranked[seq_len(cutoffs[i])]
  for (j in seq_along(pairs)) {
    results[[k]] <- analyze_pair(idx, cutoff_labels[i], pairs[[j]], comparison_names[j])
    k <- k + 1L
  }
}

stats <- rbindlist(lapply(results, `[[`, "summary"))
samples <- rbindlist(lapply(results, `[[`, "samples"), fill = TRUE)
stats[, cutoff := factor(cutoff, levels = cutoff_labels)]
stats[, comparison := factor(comparison, levels = comparison_names)]
fwrite(stats, file.path(outdir, "FEATURE_CUTOFF_1K_2K_5K_10K_STATISTICS.tsv"), sep = "\t")
fwrite(samples, file.path(outdir, "FEATURE_CUTOFF_1K_2K_5K_10K_SAMPLE_DISTANCES.tsv"), sep = "\t")

# Panel A: formal cell-line x treatment interaction with 95% CI.
p_interaction <- ggplot(stats, aes(cutoff, interaction_estimate, group = comparison)) +
  geom_hline(yintercept = 0, color = "gray40", linetype = "dashed", linewidth = 0.55) +
  geom_errorbar(aes(ymin = interaction_CI_low, ymax = interaction_CI_high),
    width = 0.12, linewidth = 0.85, color = "black") +
  geom_line(color = "black", linewidth = 0.65) +
  geom_point(aes(fill = comparison), shape = 21, color = "black",
    size = 3.8, stroke = 0.75) +
  geom_text(aes(label = sprintf("P = %.3g", interaction_p)),
    nudge_y = 0.075, size = 3.1, fontface = "bold") +
  facet_wrap(~comparison, nrow = 1, scales = "fixed") +
  scale_fill_manual(values = c("VAL vs TAZ" = "#008f8f", "VD vs TD" = "#ff84ff")) +
  labs(
    title = "Cell-line interaction across variable-gene cutoffs",
    subtitle = "Estimate = treatment difference in SW620 minus treatment difference in RKO; bars show 95% CI",
    x = NULL,
    y = "Cell-line x treatment interaction estimate"
  ) +
  guides(fill = "none") +
  theme_classic(base_size = 11) +
  theme(
    panel.border = element_rect(fill = NA, color = "black", linewidth = 0.75),
    axis.line = element_blank(),
    strip.background = element_rect(fill = "white", color = "black", linewidth = 0.75),
    strip.text = element_text(face = "bold"),
    plot.title = element_text(face = "bold"),
    axis.text.x = element_text(angle = 25, hjust = 1)
  )

# Panel B: intuitive within-line fold effect for the same models.
folds <- rbindlist(list(
  stats[, .(cutoff, comparison, cell_line = "RKO", fold = RKO_fold,
    within_line_p = RKO_within_line_p)],
  stats[, .(cutoff, comparison, cell_line = "SW620", fold = SW620_fold,
    within_line_p = SW620_within_line_p)]
))
folds[, cell_line := factor(cell_line, levels = c("RKO", "SW620"))]

p_fold <- ggplot(folds, aes(cutoff, fold, group = cell_line, color = cell_line)) +
  geom_hline(yintercept = 1, color = "gray40", linetype = "dashed", linewidth = 0.55) +
  geom_line(linewidth = 0.85) +
  geom_point(size = 3.2) +
  facet_wrap(~comparison, nrow = 1) +
  scale_color_manual(values = c(RKO = "black", SW620 = "gray60")) +
  scale_y_continuous(limits = c(1.0, 1.34), breaks = seq(1.0, 1.3, 0.1)) +
  labs(
    title = "Relative transcriptional displacement is stable across feature spaces",
    subtitle = "Fold = dual-inhibitor RMS distance / selective-inhibitor RMS distance within each cell line",
    x = "Common protein-coding genes ranked by mean within-line variance",
    y = "Fold difference in RMS distance",
    color = NULL
  ) +
  theme_classic(base_size = 11) +
  theme(
    panel.border = element_rect(fill = NA, color = "black", linewidth = 0.75),
    axis.line = element_blank(),
    strip.background = element_rect(fill = "white", color = "black", linewidth = 0.75),
    strip.text = element_text(face = "bold"),
    plot.title = element_text(face = "bold"),
    legend.position = "top",
    axis.text.x = element_text(angle = 25, hjust = 1)
  )

combined <- p_interaction / p_fold +
  plot_layout(heights = c(1.05, 1)) +
  plot_annotation(
    title = "Figure 2 sensitivity to the variable-gene cutoff",
    subtitle = "Primary analysis: top 1,000 genes; sensitivity analyses: top 2,000, 5,000 and 10,000 genes",
    tag_levels = "A",
    theme = theme(
      plot.title = element_text(face = "bold", size = 15),
      plot.subtitle = element_text(size = 10.5),
      plot.tag = element_text(face = "bold", size = 13)
    )
  )

ggsave(file.path(outdir, "FIG2_FEATURE_CUTOFF_SENSITIVITY_1K_2K_5K_10K.pdf"),
  combined, width = 10.2, height = 8.2, useDingbats = FALSE)
ggsave(file.path(outdir, "FIG2_FEATURE_CUTOFF_SENSITIVITY_1K_2K_5K_10K.png"),
  combined, width = 10.2, height = 8.2, dpi = 450, bg = "white")

print(stats)
print(combined)
