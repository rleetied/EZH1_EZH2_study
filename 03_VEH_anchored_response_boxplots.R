suppressPackageStartupMessages({library(data.table); library(ggplot2); library(patchwork)})

root <- "RKO_TAZ_v_VAL_project/Final_figure_analysis/"
indir <- file.path(root, "outputs", "VEH_ANCHORED_SINGLE_COMBO_COMBINED")
outdir <- file.path(root, "outputs", "VEH_ANCHORED_RELATIVE_ADVANTAGE_BOXPLOTS")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

d <- fread(file.path(indir, "SAMPLE_DISTANCES.tsv"))
st <- fread(file.path(indir, "WITHIN_CELL_LINE_COMPARISON_STATS.tsv"))
d[, cell_line := factor(cell_line, levels = c("RKO", "SW620"))]
d[, treatment := factor(treatment, levels = c("TAZ", "VAL", "TD", "VD"))]
cols <- c(TAZ = "#ec357a", VAL = "#008f8f", TD = "#921c52", VD = "#ff84ff")

make_panel <- function(pair, comparison_name, title, interaction_p, bracket_y, label_y) {
  z <- copy(d[treatment %chin% pair])
  # Normalize only for display. Inferential statistics and labels below come
  # from the original untransformed RMS-distance interaction models.
  ref <- z[treatment == pair[1], .(
    reference_mean = mean(RMS_distance_from_VEH)), by = cell_line]
  z <- merge(z, ref, by = "cell_line")
  z[, relative_distance := RMS_distance_from_VEH / reference_mean]
  z[, base_x := as.numeric(cell_line)]
  z[, xpos := base_x + fifelse(as.character(treatment) == pair[1], -0.18, 0.18)]

  ann <- copy(st[comparison == comparison_name])
  br_y <- bracket_y
  lab_y <- label_y
  ann[, base_x := as.numeric(factor(cell_line, levels = c("RKO", "SW620")))]
  ann[, `:=`(
    x1_plot = base_x - 0.18,
    x2_plot = base_x + 0.18,
    bracket_y_plot = br_y,
    tick_y_plot = br_y - 0.020,
    label_y_plot = lab_y,
    panel_label = sprintf("%.2fx; P = %.3g", fold, p_value)
  )]

  ggplot(z, aes(xpos, relative_distance, color = treatment, fill = treatment,
      group = interaction(cell_line, treatment))) +
    geom_hline(yintercept = 1, color = "gray55", linetype = "dashed", linewidth = 0.55) +
    geom_boxplot(outlier.shape = NA, width = 0.28, alpha = 0.12, linewidth = 0.75) +
    geom_point(position = position_jitter(width = 0.035, height = 0),
      size = 3.1, alpha = 0.80) +
    geom_segment(data = ann,
      aes(x = x1_plot, xend = x2_plot, y = bracket_y_plot, yend = bracket_y_plot),
      inherit.aes = FALSE, color = "black", linewidth = 0.65) +
    geom_segment(data = ann,
      aes(x = x1_plot, xend = x1_plot, y = tick_y_plot, yend = bracket_y_plot),
      inherit.aes = FALSE, color = "black", linewidth = 0.65) +
    geom_segment(data = ann,
      aes(x = x2_plot, xend = x2_plot, y = tick_y_plot, yend = bracket_y_plot),
      inherit.aes = FALSE, color = "black", linewidth = 0.65) +
    geom_text(data = ann,
      aes(x = base_x, y = label_y_plot, label = panel_label),
      inherit.aes = FALSE, fontface = "bold", size = 3.6) +
    scale_x_continuous(breaks = c(1, 2), labels = c("RKO", "SW620"),
      limits = c(0.55, 2.45)) +
    scale_y_continuous(breaks = c(0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4)) +
    scale_color_manual(values = cols[pair], breaks = pair) +
    scale_fill_manual(values = cols[pair], breaks = pair) +
    coord_cartesian(ylim = c(0.8, 1.4), clip = "off") +
    labs(
      title = title,
      subtitle = sprintf("Cell-line x treatment interaction: P = %.3g", interaction_p),
      x = NULL,
      y = sprintf("RMS distance relative to %s mean", pair[1]),
      color = NULL, fill = NULL
    ) +
    guides(fill = "none", color = guide_legend(override.aes = list(size = 4, alpha = 1))) +
    theme_classic(base_size = 11) +
    theme(
      aspect.ratio = 1,
      panel.border = element_rect(fill = NA, color = "black", linewidth = 0.8),
      axis.line = element_blank(),
      legend.position = "top",
      legend.justification = "left",
      plot.title = element_text(face = "bold", size = 13),
      plot.subtitle = element_text(size = 9.5),
      axis.text.x = element_text(face = "bold", color = "black")
    )
}

p_single <- make_panel(
  c("TAZ", "VAL"), "VAL/TAZ", "Single agents",
  st[comparison == "VAL/TAZ", interaction_p][1],
  bracket_y = 1.325, label_y = 1.375
)
p_combo <- make_panel(
  c("TD", "VD"), "VD/TD", "Combinations with DAC",
  st[comparison == "VD/TD", interaction_p][1],
  bracket_y = 1.325, label_y = 1.375
)

combined <- (p_single | p_combo) +
  plot_annotation(
    title = "Relative dual-inhibitor advantage across cell-line contexts",
    subtitle = paste(
      "Displayed values are normalized to the line-specific selective-inhibitor mean;",
      "P values are from models fit to untransformed RMS distances"
    ),
    tag_levels = "A",
    theme = theme(
      plot.title = element_text(face = "bold", size = 16),
      plot.subtitle = element_text(size = 10.5),
      plot.tag = element_text(face = "bold", size = 13)
    )
  )

ggsave(file.path(outdir, "RELATIVE_VAL_TAZ_AND_VD_TD_ADVANTAGE_BOXPLOTS.pdf"), combined,
  width = 9.5, height = 5.8, useDingbats = FALSE)
ggsave(file.path(outdir, "RELATIVE_VAL_TAZ_AND_VD_TD_ADVANTAGE_BOXPLOTS.png"), combined,
  width = 9.5, height = 5.8, dpi = 450, bg = "white")
fwrite(d, file.path(outdir, "SOURCE_UNTRANSFORMED_SAMPLE_DISTANCES.tsv"), sep = "\t")

print(combined)
