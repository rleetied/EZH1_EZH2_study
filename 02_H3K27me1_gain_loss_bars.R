suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

input_dir <- "RKO_TAZ_v_VAL_project/20260731_final_called_regions/H3K27me1/treatment_vs_VEH"
outdir <- file.path(
  "RKO_TAZ_v_VAL_project/Final_figure_analysis/outputs",
  "H3K27me1_TAZ_VAL_gain_loss_bars"
)
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

files <- data.table(
  Treatment = factor(c("TAZ", "VAL", "TAZ", "VAL"), levels = c("TAZ", "VAL")),
  Direction = factor(c("Increase", "Increase", "Decrease", "Decrease"),
                     levels = c("Decrease", "Increase")),
  file = c(
    "TAZ_called.GRanges_H3K27me1_gain.var_1.logFC_0.5.TAZ.bed",
    "VAL_called.GRanges_H3K27me1_gain.var_1.logFC_0.5.VAL.bed",
    "VEH_called.GRanges_H3K27me1_lost.var_1.logFC_0.5.TAZ.bed",
    "VEH_called.GRanges_H3K27me1_lost.var_1.logFC_0.5.VAL.bed"
  )
)

summarize_bed <- function(filename) {
  x <- fread(file.path(input_dir, filename), header = FALSE, select = 1:3,
             col.names = c("chrom", "start", "end"))
  stopifnot(all(x$end >= x$start))
  data.table(
    Peak_number = nrow(x),
    Coverage_bp = sum(x$end - x$start),
    Coverage_Mb = sum(x$end - x$start) / 1e6
  )
}

summary_dt <- cbind(files, rbindlist(lapply(files$file, summarize_bed)))
fwrite(summary_dt, file.path(outdir, "H3K27me1_gain_loss_peak_number_and_coverage.tsv"),
       sep = "\t")

# Explicit label positions keep very small gain segments (particularly VAL)
# readable without changing the scale or geometry of the stacked bars.
summary_dt[, total_peaks := sum(Peak_number), by = Treatment]
summary_dt[, increase_peaks := Peak_number[Direction == "Increase"], by = Treatment]
summary_dt[, peak_label_y := fifelse(
  Direction == "Increase",
  fifelse(Peak_number / total_peaks < 0.08,
          Peak_number + total_peaks * 0.035,
          Peak_number / 2),
  increase_peaks + Peak_number / 2
)]
summary_dt[, total_coverage := sum(Coverage_Mb), by = Treatment]
summary_dt[, increase_coverage := Coverage_Mb[Direction == "Increase"], by = Treatment]
summary_dt[, coverage_label_y := fifelse(
  Direction == "Increase",
  fifelse(Coverage_Mb / total_coverage < 0.08,
          Coverage_Mb + total_coverage * 0.035,
          Coverage_Mb / 2),
  increase_coverage + Coverage_Mb / 2
)]
summary_dt[, label_color := fifelse(
  Direction == "Increase" & Peak_number / total_peaks < 0.08,
  "#9E1711", "black"
)]

direction_colors <- c(
  "Increase" = "#000000",
  "Decrease" = "#B3B3B3"
)

base_theme <- theme_classic(base_size = 13) +
  theme(
    axis.title = element_text(face = "bold"),
    axis.text = element_text(color = "black"),
    axis.line = element_line(color = "black", linewidth = 0.7),
    axis.ticks = element_line(color = "black", linewidth = 0.65),
    legend.position = "top",
    legend.title = element_blank(),
    plot.title = element_text(face = "bold", size = 15),
    plot.margin = margin(12, 16, 10, 10)
  )

p_number <- ggplot(summary_dt,
                   aes(x = Treatment, y = Peak_number, fill = Direction)) +
  geom_col(width = 0.68, color = "black", linewidth = 0.55) +
  geom_text(aes(y = peak_label_y, label = scales::comma(Peak_number),
                color = label_color),
            fontface = "bold", size = 4) +
  scale_color_identity() +
  scale_fill_manual(values = direction_colors, breaks = c("Increase", "Decrease")) +
  scale_y_continuous(labels = scales::label_comma(),
                     expand = expansion(mult = c(0, 0.08))) +
  labs(
    title = "H3K27me1 differential regions",
    x = NULL,
    y = "Number of regions"
  ) +
  base_theme

p_coverage <- ggplot(summary_dt,
                     aes(x = Treatment, y = Coverage_Mb, fill = Direction)) +
  geom_col(width = 0.68, color = "black", linewidth = 0.55) +
  geom_text(aes(y = coverage_label_y, label = sprintf("%.2f", Coverage_Mb),
                color = label_color),
            fontface = "bold", size = 4) +
  scale_color_identity() +
  scale_fill_manual(values = direction_colors, breaks = c("Increase", "Decrease")) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.08))) +
  labs(
    title = "H3K27me1 genomic coverage",
    x = NULL,
    y = "Coverage (Mb)"
  ) +
  base_theme

ggsave(file.path(outdir, "H3K27me1_gain_loss_peak_number_stacked_bar.png"),
       p_number, width = 6.2, height = 5.4, dpi = 600, bg = "white")
ggsave(file.path(outdir, "H3K27me1_gain_loss_peak_number_stacked_bar.pdf"),
       p_number, width = 6.2, height = 5.4, device = "pdf",
       useDingbats = FALSE, bg = "white")
ggsave(file.path(outdir, "H3K27me1_gain_loss_coverage_Mb_stacked_bar.png"),
       p_coverage, width = 6.2, height = 5.4, dpi = 600, bg = "white")
ggsave(file.path(outdir, "H3K27me1_gain_loss_coverage_Mb_stacked_bar.pdf"),
       p_coverage, width = 6.2, height = 5.4, device = "pdf",
       useDingbats = FALSE, bg = "white")

# Diverging versions: decreases extend left of zero and increases extend right.
diverging_dt <- copy(summary_dt)
diverging_dt[, Treatment := factor(as.character(Treatment),
                                   levels = c("VAL", "TAZ"))]
diverging_dt[, Peak_number_signed := fifelse(
  Direction == "Decrease", -Peak_number, Peak_number)]
diverging_dt[, Coverage_Mb_signed := fifelse(
  Direction == "Decrease", -Coverage_Mb, Coverage_Mb)]

diverging_theme <- theme_classic(base_size = 13) +
  theme(
    axis.title = element_text(face = "bold"),
    axis.text = element_text(color = "black"),
    axis.line.y = element_blank(),
    axis.line.x = element_line(color = "black", linewidth = 0.7),
    axis.ticks = element_line(color = "black", linewidth = 0.65),
    legend.position = "bottom",
    legend.title = element_blank(),
    plot.title = element_text(face = "bold", size = 15),
    plot.margin = margin(12, 18, 10, 12)
  )

peak_limit <- max(abs(diverging_dt$Peak_number_signed)) * 1.12
p_number_diverging <- ggplot(
  diverging_dt,
  aes(y = Treatment, x = Peak_number_signed, fill = Direction)
) +
  geom_col(width = 0.62, color = "black", linewidth = 0.55) +
  geom_vline(xintercept = 0, color = "black", linewidth = 0.8) +
  geom_text(
    data = diverging_dt[Peak_number / total_peaks >= 0.08],
    aes(label = scales::comma(Peak_number), color = Direction),
    position = position_stack(vjust = 0.5),
    fontface = "bold", size = 4, show.legend = FALSE
  ) +
  scale_fill_manual(values = direction_colors, breaks = c("Decrease", "Increase")) +
  scale_color_manual(values = c("Decrease" = "black", "Increase" = "white")) +
  scale_x_continuous(
    limits = c(-peak_limit, peak_limit),
    labels = function(x) scales::comma(abs(x)),
    expand = expansion(mult = c(0.01, 0.01)),
    position = "top"
  ) +
  labs(
    title = "H3K27me1 differential regions",
    x = "Number of regions",
    y = NULL
  ) +
  diverging_theme

coverage_limit <- max(abs(diverging_dt$Coverage_Mb_signed)) * 1.12
p_coverage_diverging <- ggplot(
  diverging_dt,
  aes(y = Treatment, x = Coverage_Mb_signed, fill = Direction)
) +
  geom_col(width = 0.62, color = "black", linewidth = 0.55) +
  geom_vline(xintercept = 0, color = "black", linewidth = 0.8) +
  geom_text(
    data = diverging_dt[Coverage_Mb / total_coverage >= 0.08],
    aes(label = sprintf("%.2f", Coverage_Mb), color = Direction),
    position = position_stack(vjust = 0.5),
    fontface = "bold", size = 4, show.legend = FALSE
  ) +
  scale_fill_manual(values = direction_colors, breaks = c("Decrease", "Increase")) +
  scale_color_manual(values = c("Decrease" = "black", "Increase" = "white")) +
  scale_x_continuous(
    limits = c(-coverage_limit, coverage_limit),
    labels = function(x) abs(x),
    expand = expansion(mult = c(0.01, 0.01)),
    position = "top"
  ) +
  labs(
    title = "H3K27me1 genomic coverage",
    x = "Coverage (Mb)",
    y = NULL
  ) +
  diverging_theme

ggsave(file.path(outdir, "H3K27me1_gain_loss_peak_number_DIVERGING.png"),
       p_number_diverging, width = 8.2, height = 4.8, dpi = 600, bg = "white")
ggsave(file.path(outdir, "H3K27me1_gain_loss_peak_number_DIVERGING.pdf"),
       p_number_diverging, width = 8.2, height = 4.8, device = "pdf",
       useDingbats = FALSE, bg = "white")
ggsave(file.path(outdir, "H3K27me1_gain_loss_coverage_Mb_DIVERGING.png"),
       p_coverage_diverging, width = 8.2, height = 4.8, dpi = 600, bg = "white")
ggsave(file.path(outdir, "H3K27me1_gain_loss_coverage_Mb_DIVERGING.pdf"),
       p_coverage_diverging, width = 8.2, height = 4.8, device = "pdf",
       useDingbats = FALSE, bg = "white")
ggsave(file.path(outdir,
                 "H3K27me1_gain_loss_peak_number_DIVERGING_EDITABLE.pdf"),
       p_number_diverging, width = 8.2, height = 4.8, device = "pdf",
       useDingbats = FALSE, bg = "white")
ggsave(file.path(outdir,
                 "H3K27me1_gain_loss_coverage_Mb_DIVERGING_EDITABLE.pdf"),
       p_coverage_diverging, width = 8.2, height = 4.8, device = "pdf",
       useDingbats = FALSE, bg = "white")

cat("Saved outputs to:", outdir, "\n")
print(summary_dt[, .(Treatment, Direction, Peak_number, Coverage_Mb)])
