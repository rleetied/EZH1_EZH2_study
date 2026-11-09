#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
})

root <- "RKO_TAZ_v_VAL_project/Final_figure_analysis/"
out <- file.path(root, "outputs", "TD_VD_ATAC_GAIN_VENN_SHARED_20PCT")
dir.create(out, recursive = TRUE, showWarnings = FALSE)

membership <- read.delim(
  file.path(root, "outputs", "ATAC_treatment_overlaps", "increased",
            "all_increased_peaks_membership.tsv"),
  stringsAsFactors = FALSE, check.names = FALSE
)
stats <- read.delim(
  file.path(root, "outputs", "ATAC_VAL_DAC_specificity_20pct", "tables",
            "VAL_DAC_primary_peaks_full_classification.tsv"),
  stringsAsFactors = FALSE, check.names = FALSE
)

td <- membership$TAZ_DAC == 1
vd <- membership$VAL_DAC == 1
td_only <- td & !vd
shared <- td & vd
vd_only <- vd & !td

shared_ids <- membership$peak_id[shared]
s <- stats[match(shared_ids, stats$peak_id), ]
stopifnot(nrow(s) == sum(shared), !anyNA(s$peak_id))

# Literal 20%-greater rule requested for the shared increased-accessibility set.
# Both effects are positive because membership requires each treatment's
# significant increased-accessibility call.
s$magnitude_class <- ifelse(
  s$VAL_DAC_vs_Veh_log2FC >= 1.20 * s$TAZ_DAC_vs_Veh_log2FC,
  "VD at least 20% greater",
  "VD less than 20% greater"
)

breakdown <- as.data.frame(table(s$magnitude_class), stringsAsFactors = FALSE)
names(breakdown) <- c("category", "peaks")
breakdown$category <- factor(
  breakdown$category,
  levels = c("VD at least 20% greater", "VD less than 20% greater")
)
breakdown <- breakdown[order(breakdown$category), ]
breakdown$percent_shared <- 100 * breakdown$peaks / sum(breakdown$peaks)
breakdown$label <- sprintf("%s\n%.1f%%", scales::comma(breakdown$peaks),
                           breakdown$percent_shared)

counts <- data.frame(
  region = c("TD only", "Shared", "VD only"),
  peaks = c(sum(td_only), sum(shared), sum(vd_only))
)

# Area-proportional two-circle Venn. Circle areas encode total set sizes and
# geometric intersection encodes the observed shared count.
r1 <- sqrt(sum(td) / pi)
r2 <- sqrt(sum(vd) / pi)
overlap_area <- function(d, a, b) {
  if (d <= abs(b - a)) return(pi * min(a, b)^2)
  if (d >= a + b) return(0)
  p1 <- a^2 * acos((d^2 + a^2 - b^2) / (2 * d * a))
  p2 <- b^2 * acos((d^2 + b^2 - a^2) / (2 * d * b))
  p3 <- 0.5 * sqrt((-d + a + b) * (d + a - b) *
                   (d - a + b) * (d + a + b))
  p1 + p2 - p3
}
d <- uniroot(function(z) overlap_area(z, r1, r2) - sum(shared),
             c(abs(r2 - r1) + 1e-6, r1 + r2 - 1e-6))$root
scale_factor <- max(r1, r2)
r1 <- r1 / scale_factor
r2 <- r2 / scale_factor
d <- d / scale_factor
c1 <- -d / 2
c2 <- d / 2
circle <- function(cx, r, id) {
  a <- seq(0, 2 * pi, length.out = 500)
  data.frame(x = cx + r * cos(a), y = r * sin(a), id = id)
}
p1 <- circle(c1, r1, "TD")
p2 <- circle(c2, r2, "VD")

p_venn <- ggplot() +
  geom_polygon(data = p1, aes(x, y, group = id),
               fill = "#931100", color = "#5D0B00", alpha = 0.55,
               linewidth = 1.05) +
  geom_polygon(data = p2, aes(x, y, group = id),
               fill = "#ff84ff", color = "#A52AA5", alpha = 0.55,
               linewidth = 1.05) +
  annotate("text", x = c1 - r1 * 0.72, y = 0,
           label = scales::comma(sum(td_only)), fontface = "bold", size = 6) +
  annotate("text", x = (c1 + c2) / 2, y = 0,
           label = scales::comma(sum(shared)), fontface = "bold", size = 6) +
  annotate("text", x = c2 + r2 * 0.58, y = 0,
           label = scales::comma(sum(vd_only)), fontface = "bold", size = 6) +
  annotate("text", x = c1 - r1 * 0.35, y = r1 + 0.20,
           label = paste0("TD\n", scales::comma(sum(td)), " peaks"),
           fontface = "bold", size = 4.5, color = "#5D0B00") +
  annotate("text", x = c2 + r2 * 0.32, y = r2 + 0.28,
           label = paste0("VD\n", scales::comma(sum(vd)), " peaks"),
           fontface = "bold", size = 4.5, color = "#A52AA5") +
  coord_equal(
    xlim = c(min(c1-r1, c2-r2)-0.24, max(c1+r1, c2+r2)+0.24),
    ylim = c(-max(r1,r2)-0.25, max(r1,r2)+0.38), clip = "off"
  ) +
  labs(title = "A  Increased accessibility overlap", x = NULL, y = NULL) +
  theme_void(base_size = 12) +
  theme(plot.title = element_text(face = "bold", size = 15),
        plot.margin = margin(8, 8, 8, 8))

p_bar <- ggplot(breakdown, aes(x = 1, y = peaks, fill = category)) +
  geom_col(width = 0.62, color = "black", linewidth = 0.65) +
  geom_text(aes(label = label), position = position_stack(vjust = 0.5),
            fontface = "bold", size = 4.5,
            color = c("white", "black")) +
  scale_fill_manual(
    values = c("VD at least 20% greater" = "#A52AA5",
               "VD less than 20% greater" = "#E5E5E5"),
    labels = c("VD ≥20% greater than TD", "VD <20% greater than TD")
  ) +
  scale_y_continuous(labels = scales::comma,
                     expand = expansion(mult = c(0, 0.06))) +
  scale_x_continuous(breaks = NULL) +
  labs(title = "B  Magnitude within shared peaks",
       subtitle = paste0(scales::comma(sum(shared)),
                         " peaks increased under both TD and VD"),
       x = NULL, y = "Shared increased-accessibility peaks", fill = NULL) +
  theme_classic(base_size = 12) +
  theme(plot.title = element_text(face = "bold", size = 15),
        plot.subtitle = element_text(color = "gray30", size = 10.5),
        axis.title.y = element_text(face = "bold"),
        axis.line.x = element_blank(), axis.ticks.x = element_blank(),
        legend.position = "bottom",
        legend.text = element_text(size = 10),
        plot.margin = margin(8, 8, 8, 8))

combined <- p_venn + p_bar +
  plot_layout(widths = c(1.35, 0.85)) +
  plot_annotation(
    title = "VD amplifies accessibility at a subset of loci shared with TD",
    subtitle = paste0(
      "Significant increases: adjusted P < 0.05 and log2FC ≥ 1 versus Vehicle; ",
      "the bar applies a literal VD log2FC ≥ 1.2 × TD log2FC rule"
    ),
    theme = theme(
      plot.title = element_text(face = "bold", size = 18),
      plot.subtitle = element_text(size = 11, color = "gray25"),
      plot.background = element_rect(fill = "white", color = NA)
    )
  )

ggsave(file.path(out, "TD_VD_ATAC_GAIN_PROPORTIONAL_VENN_SHARED_20PCT.png"),
       combined, width = 13, height = 6.8, dpi = 600, bg = "white")
ggsave(file.path(out, "TD_VD_ATAC_GAIN_PROPORTIONAL_VENN_SHARED_20PCT.pdf"),
       combined, width = 13, height = 6.8, device = pdf)

write.table(counts,
            file.path(out, "TD_VD_ATAC_gain_venn_counts.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(breakdown,
            file.path(out, "TD_VD_shared_ATAC_gain_20pct_breakdown.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(
  cbind(s[, c("peak_id", "Chr", "bed_start", "bed_end",
              "TAZ_DAC_vs_Veh_log2FC", "TAZ_DAC_vs_Veh_padj",
              "VAL_DAC_vs_Veh_log2FC", "VAL_DAC_vs_Veh_padj")],
        magnitude_class = s$magnitude_class),
  file.path(out, "TD_VD_shared_ATAC_gain_peak_level_20pct_classification.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

cat("Saved to:", out, "\n")
print(counts)
print(breakdown)
