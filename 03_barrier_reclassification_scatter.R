#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

root <- "RKO_TAZ_v_VAL_project/Final_figure_analysis/"
parent <- file.path(root, "outputs", "FINAL_UP_DOWN_rGREAT_RNA_PRIORITY_LINKS",
  "Combo_expression_scatter", "barrier_reclassification")
expression_treatment <- Sys.getenv("EXPRESSION_TREATMENT", "VD")
stopifnot(expression_treatment %in% c("TD", "VD"))
accessibility_treatment <- Sys.getenv("ACCESSIBILITY_TREATMENT", "ORIGINAL")
stopifnot(accessibility_treatment %in% c("ORIGINAL", "TD", "VD"))
outdir <- file.path(parent, "combined_K27me1_induced_bivalency",
  paste0(expression_treatment,
    "_EXPRESSION_", accessibility_treatment,
    "_ACCESSIBILITY_FOR_BOTH_CLASSES"))
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

d <- fread(file.path(parent, "barrier_reclassified_complete_data.tsv.gz"))
d[barrier_class %in% c(
  "Adaptive H3K27me1 + induced bivalency",
  "H3K27me1 + induced bivalency"
), barrier_class := "H3K27me1 + induced bivalency"]

# Optionally replace the original class-defining ATAC value with a single
# treatment's ATAC response at the same representative peak. The row grouping
# remains VD-favored versus Comparable, but both axes can now reflect TD.
if (accessibility_treatment != "ORIGINAL") {
  stats <- fread(file.path(root, "outputs", "ATAC_VAL_DAC_specificity_20pct",
    "tables", "VAL_DAC_primary_peaks_full_classification.tsv"))
  atac_col <- if (accessibility_treatment == "TD") {
    "TAZ_DAC_vs_Veh_log2FC"
  } else {
    "VAL_DAC_vs_Veh_log2FC"
  }
  atac <- unique(stats[, .(peak_id,
    Treatment_ATAC_log2FC = get(atac_col))], by = "peak_id")
  d <- merge(d, atac, by = "peak_id", all.x = TRUE, sort = FALSE)
  d[, ATAC_log2FC := Treatment_ATAC_log2FC]
}

# Use VD-vs-VEH RNA expression for both accessibility-response strata. This
# keeps the y-axis biology constant while the rows differ only in whether the
# linked ATAC response was VD-favored or comparable between TD and VD.
rna_file <- paste0(
  "RKO_TAZ_v_VAL_project/",
  "PR001414_ROTS_RNA/trimmed_fastq/count_matrix/",
  "RKO_TAZ_VAL_RNA.0.25_norm_counts.stats.ann.txt"
)
rna <- as.data.table(read.table(rna_file, header = TRUE, check.names = FALSE,
  stringsAsFactors = FALSE, quote = "\""))
ambiguous <- rna[, .N, by = ext_gene][N > 1L, ext_gene]
rna_expression <- rna[!ext_gene %in% ambiguous,
  .(gene = ext_gene,
    Expression_RNA_log2FC = get(paste0("LogFC.", expression_treatment)),
    Expression_RNA_FDR = get(paste0("FDR.", expression_treatment)),
    Expression_RNA_logCPM = get(paste0("logCPM.", expression_treatment)))]
d <- merge(d, rna_expression, by = "gene", all.x = TRUE, sort = FALSE)
d[, `:=`(
  RNA_log2FC = Expression_RNA_log2FC,
  RNA_FDR = Expression_RNA_FDR,
  RNA_logCPM = Expression_RNA_logCPM,
  Treatment = expression_treatment,
  displayed = Expression_RNA_log2FC <= -0.5 | Expression_RNA_log2FC >= 1,
  RNA_FDR_significant = !is.na(Expression_RNA_FDR) &
    Expression_RNA_FDR < 0.05 &
    (Expression_RNA_log2FC <= -0.5 | Expression_RNA_log2FC >= 1)
)]
d[accessibility_panel == "VD", accessibility_panel := "VD-favored"]
d[accessibility_panel == "Combo", accessibility_panel := "Comparable"]

class_levels <- c(
  "Adaptive H3K27me1 barrier",
  "H3K27me1 + induced bivalency",
  "Induced bivalency",
  "H3K27me1 loss"
)
d[, barrier_class := factor(barrier_class, levels = class_levels)]
d[, accessibility_panel := factor(accessibility_panel,
  levels = c("VD-favored", "Comparable"))]
d[, regulatory_class := factor(regulatory_class, levels = c("TSS", "CRE"))]
plot_d <- d[displayed == TRUE]

# Category-specific transparency: soften the abundant adaptive-barrier signal
# and strengthen the combined-barrier signal at final figure size.
plot_d[, base_alpha := fcase(
  barrier_class == "Adaptive H3K27me1 barrier", 0.10,
  barrier_class == "H3K27me1 + induced bivalency", 0.34,
  default = 0.22
)]
plot_d[, significant_alpha := fcase(
  barrier_class == "Adaptive H3K27me1 barrier", 0.21,
  barrier_class == "H3K27me1 + induced bivalency", 0.58,
  default = 0.38
)]

colors <- c(
  "Adaptive H3K27me1 barrier" = "#B5177E",
  "H3K27me1 + induced bivalency" = "#B8860B",
  "Induced bivalency" = "#0072B2",
  "H3K27me1 loss" = "#5C5C5C"
)

plot_d[, quadrant := fcase(
  ATAC_log2FC <= -1 & RNA_log2FC >= 1, "ATAC loss / RNA gain",
  ATAC_log2FC <= -1 & RNA_log2FC <= -0.5, "ATAC loss / RNA loss",
  ATAC_log2FC >= 1 & RNA_log2FC >= 1, "ATAC gain / RNA gain",
  ATAC_log2FC >= 1 & RNA_log2FC <= -0.5, "ATAC gain / RNA loss"
)]
qc <- plot_d[, .(n = uniqueN(gene)),
  by = .(accessibility_panel, regulatory_class, quadrant)]
aq <- CJ(
  accessibility_panel = levels(d$accessibility_panel),
  regulatory_class = levels(d$regulatory_class),
  quadrant = c("ATAC loss / RNA gain", "ATAC loss / RNA loss",
               "ATAC gain / RNA gain", "ATAC gain / RNA loss"), unique = TRUE)
qc <- merge(aq, qc,
  by = c("accessibility_panel", "regulatory_class", "quadrant"), all.x = TRUE)
qc[is.na(n), n := 0L]
qc[, `:=`(
  label_x = fifelse(grepl("ATAC loss", quadrant), -3.25, 5.05),
  label_y = fifelse(grepl("RNA gain", quadrant), 11.35, -3.45),
  hjust_value = fifelse(grepl("ATAC loss", quadrant), 0, 1),
  vjust_value = fifelse(grepl("RNA gain", quadrant), 1, 0),
  label = paste0("n = ", format(n, big.mark = ","))
)]
qc[, accessibility_panel := factor(accessibility_panel,
  levels = levels(d$accessibility_panel))]
qc[, regulatory_class := factor(regulatory_class,
  levels = levels(d$regulatory_class))]

summary <- d[, .(
  total_genes = uniqueN(gene), displayed_genes = sum(displayed),
  displayed_RNA_FDR_lt_0.05 = sum(RNA_FDR_significant)
), by = .(barrier_class, accessibility_panel, regulatory_class)]
fwrite(d, file.path(outdir,
  "combined_barrier_reclassified_complete_data.tsv.gz"), sep = "\t")
fwrite(summary, file.path(outdir,
  "combined_barrier_reclassification_counts.tsv"), sep = "\t")

p <- ggplot(plot_d, aes(ATAC_log2FC, RNA_log2FC)) +
  geom_vline(xintercept = 0, color = "gray18", linewidth = 0.58) +
  geom_hline(yintercept = 0, color = "gray18", linewidth = 0.58) +
  geom_vline(xintercept = c(-1, 1), color = "gray48",
    linetype = "dashed", linewidth = 0.48) +
  geom_hline(yintercept = c(-0.5, 1), color = "gray48",
    linetype = "dashed", linewidth = 0.48) +
  geom_point(aes(color = barrier_class, alpha = base_alpha), size = 2.20,
    shape = 16) +
  geom_point(data = plot_d[RNA_FDR_significant == TRUE],
    aes(color = barrier_class, alpha = significant_alpha), size = 3.10,
    shape = 16) +
  geom_label(data = qc,
    aes(x = label_x, y = label_y, label = label,
      hjust = hjust_value, vjust = vjust_value),
    inherit.aes = FALSE, size = 3.15, fontface = "bold",
    linewidth = 0, fill = scales::alpha("white", 0.72),
    label.padding = grid::unit(0.10, "lines"), color = "gray15") +
  facet_grid(rows = vars(accessibility_panel),
    cols = vars(regulatory_class), drop = FALSE, switch = "y") +
  scale_color_manual(values = colors, drop = FALSE) +
  scale_alpha_identity() +
  scale_x_continuous(breaks = c(-3, -1, 0, 1, 3, 5),
    expand = expansion(mult = c(0.02, 0.04))) +
  scale_y_continuous(breaks = c(-4, -2, -1, 0, 1, 4, 8, 12),
    expand = expansion(mult = c(0.02, 0.03)), position = "right") +
  coord_cartesian(xlim = c(-3.5, 5.35), ylim = c(-4, 12), clip = "on") +
  labs(title = paste0("Chromatin-barrier classification | ",
      expression_treatment, " accessibility and expression"),
    x = paste0(accessibility_treatment,
      " Accessibility (log2 FC vs VEH)"),
    y = paste0(expression_treatment,
      " Gene Expression (log2 FC vs VEH)"), color = NULL) +
  guides(color = guide_legend(
    override.aes = list(size = 3.0, alpha = 0.72), ncol = 2)) +
  theme_classic(base_size = 12) +
  theme(plot.title = element_text(hjust = 0.5, size = 14),
    strip.background = element_rect(fill = "white", color = "black", linewidth = 0.8),
    strip.text = element_text(face = "bold", size = 11),
    strip.placement = "outside",
    panel.border = element_rect(fill = NA, color = "black", linewidth = 0.85),
    panel.spacing = grid::unit(0.08, "lines"), axis.line = element_blank(),
    axis.ticks = element_line(color = "black", linewidth = 0.65),
    axis.text = element_text(color = "black", size = 9.5),
    axis.title = element_text(color = "black", size = 12),
    legend.position = "top", legend.text = element_text(size = 9.4),
    legend.margin = margin(t = -4, b = 0), plot.margin = margin(8, 10, 8, 8))

png <- file.path(outdir, "FINAL_barrier_reclassified_COMBINED_bivalency_FACET_LEFT_AXIS_RIGHT.png")
pdf <- file.path(outdir, "FINAL_barrier_reclassified_COMBINED_bivalency_FACET_LEFT_AXIS_RIGHT.pdf")
ggsave(png, p, width = 7.2, height = 6.5, dpi = 600, bg = "white")
ggsave(pdf, p, width = 7.2, height = 6.5, device = "pdf",
  useDingbats = FALSE, bg = "white")

p0 <- p
p0$layers <- Filter(function(layer) !inherits(layer$geom, "GeomLabel"), p0$layers)
p0 <- p0 + labs(title = NULL, x = NULL, y = NULL, color = NULL) +
  theme(legend.position = "none", strip.text = element_blank(),
    axis.text = element_blank(), axis.title = element_blank(),
    plot.title = element_blank())
textless <- file.path(outdir,
  "FINAL_barrier_reclassified_COMBINED_bivalency_FACET_LEFT_AXIS_RIGHT_TEXTLESS.png")
ggsave(textless, p0, width = 6.2, height = 6.2, dpi = 600, bg = "white")
textless_pdf <- file.path(outdir,
  "FINAL_barrier_reclassified_COMBINED_bivalency_ILLUSTRATOR_TEXTLESS.pdf")
ggsave(textless_pdf, p0, width = 6.2, height = 6.2, device = "pdf",
  useDingbats = FALSE, bg = "white")

cat("Saved combined-bivalency version to:\n", outdir, "\n")
print(summary)
