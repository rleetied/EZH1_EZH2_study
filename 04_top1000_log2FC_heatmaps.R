suppressPackageStartupMessages({library(data.table); library(pheatmap)})

root <- "RKO_TAZ_v_VAL_project/Final_figure_analysis/"
outdir <- file.path(root, "outputs", "FEATURE_CUTOFF_SENSITIVITY_AND_HEATMAPS")
paths <- c(
  RKO = "RKO_TAZ_v_VAL_project/RNA-seq/RKO/RKO.DAC_TAZ_VAL_treatments_D3.ann.stats.0.25.txt",
  SW620 = "RKO_TAZ_v_VAL_project/RNA-seq/SW620/SW620.DAC_TAZ_VAL_treatments_D3.ann.stats.0.25.txt"
)
tx <- c("Veh", "DAC", "TAZ", "VAL", "TD", "VD")
sample_cols <- unlist(lapply(tx, function(z) paste0(z, 1:3)))

read_one <- function(path) fread(path)[gene_biotype == "protein_coding" & !is.na(ENSG)]
rko <- read_one(paths[["RKO"]])
sw <- read_one(paths[["SW620"]])
common <- intersect(rko$ENSG, sw$ENSG)
rko <- rko[match(common, ENSG)]
sw <- sw[match(common, ENSG)]

make_log_expr <- function(x) {
  m <- log2(as.matrix(x[, ..sample_cols]) + 0.25)
  rownames(m) <- x$ENSG
  m
}
er <- make_log_expr(rko)
es <- make_log_expr(sw)

# Recreate the joint within-line variance ranking used throughout the analysis.
mean_variance <- (apply(er, 1, var) + apply(es, 1, var)) / 2
top_ensg <- rko$ENSG[order(mean_variance, decreasing = TRUE)[1:1000]]
er <- er[top_ensg, , drop = FALSE]
es <- es[top_ensg, , drop = FALSE]

# Row-wise z-scores across all 18 samples in each cell line. Because subtracting
# a gene's VEH mean is a constant shift, z-scoring the log2 expression values is
# equivalent to z-scoring its sample-resolved log2FC values, while retaining the
# three VEH replicates as the baseline pattern.
r_fc <- t(scale(t(er)))
s_fc <- t(scale(t(es)))
r_fc[!is.finite(r_fc)] <- 0
s_fc[!is.finite(s_fc)] <- 0

# Derive one row order from both cell lines and use it for both heatmaps.
order_matrix <- cbind(r_fc, s_fc)
order_matrix[!is.finite(order_matrix)] <- 0
row_order <- hclust(dist(order_matrix), method = "ward.D2")$order
r_fc <- r_fc[row_order, , drop = FALSE]
s_fc <- s_fc[row_order, , drop = FALSE]

cap <- 2.5
clip <- function(m) {
  out <- m
  out[out < -cap] <- -cap
  out[out > cap] <- cap
  out
}
colors <- colorRampPalette(c("#2166ac", "#f7f7f7", "#b2182b"))(100)
breaks <- seq(-cap, cap, length.out = 101)

treatment <- sub("[123]$", "", sample_cols)
annotation_col <- data.frame(Treatment = factor(treatment, levels = tx), row.names = sample_cols)
annotation_colors <- list(Treatment = c(
  Veh = "#777777", DAC = "#00c7f5", TAZ = "#ec357a",
  VAL = "#008f8f", TD = "#921c52", VD = "#ff84ff"
))

heat_args <- list(
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  show_rownames = FALSE,
  border_color = NA,
  color = colors,
  breaks = breaks,
  annotation_col = annotation_col,
  annotation_colors = annotation_colors,
  gaps_col = c(3, 6, 9, 12, 15),
  fontsize_col = 8,
  angle_col = 45,
  legend_breaks = c(-2.5, 0, 2.5),
  legend_labels = c("<= -2.5", "0", ">= 2.5")
)

save_heatmap <- function(mat, line) {
  title <- sprintf("%s: row-z-scored expression across all 18 samples", line)
  pdf(file.path(outdir, sprintf("TOP1000_ALL_18_SAMPLES_ROW_ZSCORE_HEATMAP_%s.pdf", line)),
    width = 9.0, height = 9.5, useDingbats = FALSE)
  do.call(pheatmap, c(list(mat = clip(mat), main = title), heat_args))
  dev.off()
  png(file.path(outdir, sprintf("TOP1000_ALL_18_SAMPLES_ROW_ZSCORE_HEATMAP_%s.png", line)),
    width = 4050, height = 4275, res = 450, bg = "white")
  do.call(pheatmap, c(list(mat = clip(mat), main = title), heat_args))
  dev.off()
}

save_heatmap(r_fc, "RKO")
save_heatmap(s_fc, "SW620")

row_key <- data.table(
  heatmap_row = seq_len(1000),
  ENSG = rownames(r_fc),
  gene_symbol = rko$ext_gene[match(rownames(r_fc), rko$ENSG)]
)
fwrite(row_key, file.path(outdir, "TOP1000_ALL_SAMPLE_HEATMAP_ROW_ORDER.tsv"), sep = "\t")
fwrite(cbind(row_key, as.data.table(r_fc)),
  file.path(outdir, "TOP1000_ALL_18_SAMPLES_ROW_ZSCORE_MATRIX_RKO.tsv"), sep = "\t")
fwrite(cbind(row_key, as.data.table(s_fc)),
  file.path(outdir, "TOP1000_ALL_18_SAMPLES_ROW_ZSCORE_MATRIX_SW620.tsv"), sep = "\t")
