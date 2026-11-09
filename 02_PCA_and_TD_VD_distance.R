suppressPackageStartupMessages({library(data.table); library(ggplot2); library(patchwork)})

root <- "RKO_TAZ_v_VAL_project/Final_figure_analysis/"
outdir <- file.path(root, "outputs", "RKO_SW620_PCA_TD_VD_DISTANCE")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
paths <- c(
  RKO = "RKO_TAZ_v_VAL_project/RNA-seq/RKO/RKO.DAC_TAZ_VAL_treatments_D3.ann.stats.0.25.txt",
  SW620 = "RKO_TAZ_v_VAL_project/RNA-seq/SW620/SW620.DAC_TAZ_VAL_treatments_D3.ann.stats.0.25.txt"
)
tx <- c("Veh", "DAC", "TAZ", "VAL", "TD", "VD")
sample_cols <- unlist(lapply(tx, function(z) paste0(z, 1:3)))
tx_cols <- c(Veh = "#777777", DAC = "#00c7f5", TAZ = "#ec357a",
  VAL = "#008f8f", TD = "#921c52", VD = "#ff84ff")

read_one <- function(path) {
  x <- fread(path)
  x[gene_biotype == "protein_coding" & !is.na(ENSG)]
}
rko <- read_one(paths[["RKO"]]); sw <- read_one(paths[["SW620"]])
common <- intersect(rko$ENSG, sw$ENSG)
rko <- rko[match(common, ENSG)]; sw <- sw[match(common, ENSG)]
stopifnot(all(rko$ENSG == sw$ENSG))

make_expr <- function(x, line) {
  m <- log2(as.matrix(x[, ..sample_cols]) + 0.25)
  rownames(m) <- x$ENSG; colnames(m) <- paste(line, sample_cols, sep = "_")
  m
}
er <- make_expr(rko, "RKO"); es <- make_expr(sw, "SW620")

# Select one common feature space using mean within-line variance.
v <- (apply(er, 1, var) + apply(es, 1, var)) / 2
sel <- names(sort(v, decreasing = TRUE))[1:1000]

meta <- data.table(sample = c(colnames(er), colnames(es)))
meta[, c("cell_line", "sample_short") := tstrsplit(sample, "_", fixed = TRUE)]
meta[, treatment := sub("[123]$", "", sample_short)]
meta[, replicate := as.integer(sub(".*([123])$", "\\1", sample_short))]
meta[, treatment := factor(treatment, levels = tx)]

analyze_line <- function(m, line) {
  m <- m[sel, , drop = FALSE]
  pc <- prcomp(t(m), center = TRUE, scale. = FALSE)
  ve <- 100 * pc$sdev^2 / sum(pc$sdev^2)
  scores <- data.table(sample = rownames(pc$x), PC1 = pc$x[, 1], PC2 = pc$x[, 2])
  scores <- merge(scores, meta[cell_line == line], by = "sample")
  cent <- scores[, .(PC1 = mean(PC1), PC2 = mean(PC2)), by = treatment]
  tdvd <- cent[treatment %chin% c("TD", "VD")]
  pc2d <- sqrt(diff(tdvd$PC1)^2 + diff(tdvd$PC2)^2)
  td <- rowMeans(m[, grepl("_TD[123]$", colnames(m)), drop = FALSE])
  vd <- rowMeans(m[, grepl("_VD[123]$", colnames(m)), drop = FALSE])
  full_d <- sqrt(sum((vd - td)^2))
  rms_d <- sqrt(mean((vd - td)^2))
  # Mean replicate-to-centroid distance in the same 1,000-gene space.
  within <- c(
    apply(m[, grepl("_TD[123]$", colnames(m)), drop = FALSE], 2,
      function(z) sqrt(sum((z - td)^2))),
    apply(m[, grepl("_VD[123]$", colnames(m)), drop = FALSE], 2,
      function(z) sqrt(sum((z - vd)^2)))
  )
  list(scores = scores, centroids = cent, ve = ve,
    stats = data.table(cell_line = line, PC1_PC2_distance = pc2d,
      full_1000_gene_euclidean = full_d, RMS_gene_distance = rms_d,
      pooled_within_dispersion = mean(within),
      separation_ratio = full_d / mean(within)))
}
a <- analyze_line(er, "RKO"); b <- analyze_line(es, "SW620")
stats <- rbindlist(list(a$stats, b$stats))
fold <- stats[cell_line == "RKO", RMS_gene_distance] / stats[cell_line == "SW620", RMS_gene_distance]

# Use identical axis limits so visual separation is directly comparable.
all_scores <- rbindlist(list(a$scores, b$scores), fill = TRUE)
pad_range <- function(x, fraction = 0.06) {
  r <- range(x, finite = TRUE)
  r + c(-1, 1) * diff(r) * fraction
}
common_xlim <- pad_range(all_scores$PC1)
common_ylim <- pad_range(all_scores$PC2)

make_pca <- function(z, line) {
  scores <- z$scores; cent <- z$centroids
  tdvd <- cent[treatment %chin% c("TD", "VD")]
  lab <- sprintf("TD-VD centroid distance\nPC1/2 = %.2f; full = %.2f",
    z$stats$PC1_PC2_distance, z$stats$full_1000_gene_euclidean)
  ggplot(scores, aes(PC1, PC2, color = treatment)) +
    geom_point(data = scores[!treatment %chin% c("TD", "VD")], size = 2.5, alpha = 0.28) +
    geom_point(data = scores[treatment %chin% c("TD", "VD")], size = 4.2, alpha = 0.92) +
    geom_segment(data = data.table(x = tdvd[treatment == "TD", PC1],
        y = tdvd[treatment == "TD", PC2], xend = tdvd[treatment == "VD", PC1],
        yend = tdvd[treatment == "VD", PC2]),
      aes(x = x, y = y, xend = xend, yend = yend), inherit.aes = FALSE,
      color = "black", linewidth = 1.1, arrow = arrow(length = unit(0.13, "inches"))) +
    geom_point(data = tdvd, aes(PC1, PC2, fill = treatment), shape = 23,
      color = "black", size = 5, stroke = 1.0) +
    geom_label(aes(x = -Inf, y = Inf, label = lab), inherit.aes = FALSE,
      hjust = -0.04, vjust = 1.08, size = 3.1, fontface = "bold",
      fill = "white", linewidth = 0) +
    scale_color_manual(values = tx_cols) + scale_fill_manual(values = tx_cols) +
    coord_cartesian(xlim = common_xlim, ylim = common_ylim, clip = "off") +
    labs(title = line, x = sprintf("PC1 (%.1f%%)", z$ve[1]),
      y = sprintf("PC2 (%.1f%%)", z$ve[2]), color = NULL) +
    theme_classic(base_size = 11) +
    theme(panel.border = element_rect(fill = NA, color = "black"), axis.line = element_blank(),
      plot.title = element_text(face = "bold"), legend.position = "none")
}
p_rko <- make_pca(a, "RKO")
p_sw <- make_pca(b, "SW620")

p_dist <- ggplot(stats, aes(cell_line, RMS_gene_distance, fill = cell_line)) +
  geom_col(width = 0.58, color = "black", linewidth = 0.7) +
  geom_text(aes(label = sprintf("%.3f", RMS_gene_distance)), vjust = -0.55,
    size = 4.2, fontface = "bold") +
  scale_fill_manual(values = c(RKO = "black", SW620 = "gray65")) +
  coord_cartesian(ylim = c(0, max(stats$RMS_gene_distance) * 1.22), clip = "off") +
  labs(title = sprintf("VD-TD transcriptome separation is %.2f-fold greater in RKO", fold),
    subtitle = "Euclidean centroid distance normalized per gene; identical 1,000-gene feature space",
    x = NULL, y = "RMS VD-TD centroid distance") +
  guides(fill = "none") + theme_classic(base_size = 11) +
  theme(panel.border = element_rect(fill = NA, color = "black"), axis.line = element_blank(),
    plot.title = element_text(face = "bold"))

combined <- (p_rko | p_sw) / p_dist +
  plot_layout(heights = c(1, 0.72)) +
  plot_annotation(title = "RNA-seq treatment structure and VD-TD separation",
    subtitle = "PCA uses the same 1,000 common, most-variable protein-coding genes in both cell lines",
    tag_levels = "A", theme = theme(plot.title = element_text(face = "bold", size = 16)))

ggsave(file.path(outdir, "RKO_SW620_PCA_TD_VD_DISTANCE.pdf"), combined,
  width = 10, height = 8.5, useDingbats = FALSE)
ggsave(file.path(outdir, "RKO_SW620_PCA_TD_VD_DISTANCE.png"), combined,
  width = 10, height = 8.5, dpi = 450, bg = "white")
fwrite(stats, file.path(outdir, "TD_VD_EUCLIDEAN_DISTANCE_STATISTICS.tsv"), sep = "\t")
fwrite(rbindlist(list(a$scores, b$scores), fill = TRUE),
  file.path(outdir, "PCA_SAMPLE_COORDINATES.tsv"), sep = "\t")
print(stats)
cat(sprintf("RKO/SW620 RMS distance ratio: %.3f\n", fold))
