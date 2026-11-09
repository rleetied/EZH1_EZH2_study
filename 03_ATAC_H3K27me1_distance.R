#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
})

root <- getwd()
outdir <- file.path(root, "outputs/VAL_BIASED_ATAC_GAIN_near_TAZ_H3K27me1")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

input_file <- file.path(root,
  "H3K27me1_ATAC_TAZ_VAL_integration/results/distance_windows",
  "ATAC_increased_H3K27me1_distance_windows_percent.tsv")

context_levels <- c(
  "TAZ-specific H3K27me1 gain",
  "TAZ-maintained / VAL-only H3K27me1 loss"
)
windows_keep <- c(seq(0, 1000, 100), seq(2000, 10000, 1000))
cols <- c(TAZ = "#ec357a", VAL = "#008f8f")

d <- fread(input_file)[
  treatment %in% c("TAZ", "VAL") &
    H3K27me1_class %in% c("H3_TAZ_gain", "H3_VAL_loss") &
    window_bp %in% windows_keep
]
d[, relationship := fifelse(
  H3K27me1_class == "H3_TAZ_gain",
  context_levels[1], context_levels[2]
)]
d[, relationship := factor(relationship, levels = context_levels)]
d[, treatment := factor(treatment, levels = c("TAZ", "VAL"))]

# Endpoint statistics: OR quantifies relative enrichment; count ratio quantifies magnitude.
endpoint <- d[window_bp == 10000L]
stats <- rbindlist(lapply(context_levels, function(rel) {
  z <- endpoint[relationship == rel]
  val <- z[treatment == "VAL"]
  taz <- z[treatment == "TAZ"]
  ft <- fisher.test(matrix(c(
    val$n_ATAC_within, val$n_ATAC_total - val$n_ATAC_within,
    taz$n_ATAC_within, taz$n_ATAC_total - taz$n_ATAC_within
  ), nrow = 2, byrow = TRUE))
  data.table(
    relationship = rel,
    odds_ratio = unname(ft$estimate),
    CI_lower = unname(ft$conf.int[1]),
    CI_upper = unname(ft$conf.int[2]),
    fisher_p_value = ft$p.value,
    val_count = val$n_ATAC_within,
    taz_count = taz$n_ATAC_within,
    raw_count_ratio = val$n_ATAC_within / taz$n_ATAC_within,
    total_gain_ratio = val$n_ATAC_total / taz$n_ATAC_total
  )
}))
stats[, relationship := factor(relationship, levels = context_levels)]
stats[, percent_label := sprintf(
  "10-kb OR %.2f (95%% CI %.2f-%.2f)\nFisher p = %s",
  odds_ratio, CI_lower, CI_upper,
  fifelse(fisher_p_value < 0.001,
          format(fisher_p_value, digits = 2, scientific = TRUE),
          sprintf("%.3f", fisher_p_value))
)]
stats[, count_label := sprintf(
  "10 kb: VAL %s vs TAZ %s\n%.2f-fold more VAL peaks",
  scales::comma(val_count), scales::comma(taz_count), raw_count_ratio
)]

percent_ann <- merge(
  stats[, .(relationship, label = percent_label)],
  d[, .(y = max(pct_ATAC_within) * 0.96), by = relationship],
  by = "relationship"
)
count_ann <- merge(
  stats[, .(relationship, label = count_label)],
  d[, .(y = max(n_ATAC_within) * 0.96), by = relationship],
  by = "relationship"
)

base_theme <- theme_classic(base_size = 11.5) +
  theme(
    strip.background = element_rect(fill = "white", color = "black", linewidth = 0.7),
    strip.text = element_text(face = "bold"),
    axis.line = element_line(color = "black", linewidth = 0.7),
    axis.ticks = element_line(color = "black", linewidth = 0.65),
    axis.text = element_text(color = "black"),
    plot.title = element_text(face = "bold"),
    legend.position = "top",
    legend.title = element_blank()
  )

distance_scale <- scale_x_continuous(
  breaks = c(0, 1000, 5000, 10000),
  labels = c("Direct\noverlap", "1 kb", "5 kb", "10 kb")
)

p_pct <- ggplot(d, aes(window_bp, pct_ATAC_within,
                       color = treatment, group = treatment)) +
  geom_line(linewidth = 1.0) +
  geom_point(size = 2.0) +
  geom_label(data = percent_ann,
             aes(x = 5200, y = y, label = label),
             inherit.aes = FALSE, size = 3.2, fontface = "bold",
             label.size = 0, fill = scales::alpha("white", 0.82)) +
  facet_wrap(~relationship, nrow = 1, scales = "free_y") +
  scale_color_manual(values = cols) + distance_scale +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.15))) +
  labs(
    x = NULL,
    y = "Increased ATAC peaks (%)",
    title = "Relative enrichment",
    subtitle = "Fraction of each treatment's total accessibility gains near the H3K27me1 state"
  ) + base_theme

p_count <- ggplot(d, aes(window_bp, n_ATAC_within,
                         color = treatment, group = treatment)) +
  geom_line(linewidth = 1.0) +
  geom_point(size = 2.0) +
  geom_label(data = count_ann,
             aes(x = 5200, y = y, label = label),
             inherit.aes = FALSE, size = 3.2, fontface = "bold",
             label.size = 0, fill = scales::alpha("white", 0.82)) +
  facet_wrap(~relationship, nrow = 1, scales = "free_y") +
  scale_color_manual(values = cols) + distance_scale +
  scale_y_continuous(labels = scales::label_comma(),
                     expand = expansion(mult = c(0.02, 0.15))) +
  labs(
    x = "Maximum edge-to-edge distance",
    y = "Number of increased ATAC peaks",
    title = "Absolute magnitude",
    subtitle = "Cumulative raw counts preserve VAL's substantially broader accessibility response"
  ) + base_theme + guides(color = "none")

combined <- p_pct / p_count +
  plot_annotation(
    title = "VAL accessibility gains are both broader and selectively enriched near the adaptive H3K27me1 state",
    subtitle = paste0(
      "Percentages and odds ratios test relative association; raw counts show response magnitude. ",
      "Both views use the same cumulative distance windows."
    ),
    tag_levels = "A",
    theme = theme(
      plot.title = element_text(face = "bold", size = 15),
      plot.subtitle = element_text(size = 10.5),
      plot.tag = element_text(face = "bold", size = 13)
    )
  )

stem <- "ATAC_H3K27me1_DISTANCE_PERCENT_and_RAW_MAGNITUDE"
ggsave(file.path(outdir, paste0(stem, ".png")), combined,
       width = 12.5, height = 9.2, dpi = 600, bg = "white")
ggsave(file.path(outdir, paste0(stem, ".pdf")), combined,
       width = 12.5, height = 9.2, device = pdf)

fwrite(stats, file.path(outdir,
  "ATAC_H3K27me1_10kb_OR_and_magnitude_summary.tsv"), sep = "\t")
cat("Saved percentage + magnitude figure to:", outdir, "\n")
print(stats)
