suppressPackageStartupMessages({library(data.table); library(ggplot2); library(patchwork)})

root <- "RKO_TAZ_v_VAL_project/Final_figure_analysis/"
outdir <- file.path(root, "outputs", "EZH1_EZH2_TRANSCRIPTS_RKO_SW620")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
paths <- c(
  RKO = "RKO_TAZ_v_VAL_project/RNA-seq/RKO/RKO.DAC_TAZ_VAL_treatments_D3.ann.stats.0.25.txt",
  SW620 = "RKO_TAZ_v_VAL_project/RNA-seq/SW620/SW620.DAC_TAZ_VAL_treatments_D3.ann.stats.0.25.txt"
)
tx <- c("Veh", "DAC", "TAZ", "VAL", "TD", "VD")
sample_cols <- unlist(lapply(tx, function(z) paste0(z, 1:3)))
tx_cols <- c(Veh = "#777777", DAC = "#00c7f5", TAZ = "#ec357a",
  VAL = "#008f8f", TD = "#921c52", VD = "#ff84ff")

read_one <- function(path, line) {
  x <- fread(path)
  x <- x[ext_gene %chin% c("EZH1", "EZH2")]
  stopifnot(nrow(x) == 2)
  x[, cell_line := line]
  x
}
wide <- rbindlist(Map(read_one, paths, names(paths)), fill = TRUE)

expr <- melt(wide,
  id.vars = c("cell_line", "ENSG", "ext_gene"),
  measure.vars = sample_cols,
  variable.name = "sample", value.name = "normalized_expression")
expr[, treatment := sub("[123]$", "", sample)]
expr[, replicate := as.integer(sub(".*([123])$", "\\1", sample))]
expr[, treatment := factor(treatment, levels = tx)]
expr[, log2_expression := log2(normalized_expression + 0.25)]

means <- expr[, .(mean_log2 = mean(log2_expression),
  sd_log2 = sd(log2_expression)), by = .(cell_line, ext_gene, treatment)]

p_abs <- ggplot(expr, aes(treatment, log2_expression, color = treatment, group = 1)) +
  geom_line(data = means, aes(y = mean_log2), color = "gray35", linewidth = 0.65) +
  geom_point(position = position_jitter(width = 0.08, height = 0), size = 2.7, alpha = 0.75) +
  geom_point(data = means, aes(y = mean_log2), shape = 21, fill = "white",
    size = 3.2, stroke = 0.9) +
  facet_grid(ext_gene ~ cell_line) +
  scale_color_manual(values = tx_cols) +
  labs(title = "EZH1 and EZH2 transcript abundance across treatments",
    subtitle = "Points are biological replicates; open circles indicate treatment means",
    x = NULL, y = "log2(normalized expression + 0.25)") +
  guides(color = "none") + theme_classic(base_size = 11) +
  theme(panel.border = element_rect(fill = NA, color = "black"), axis.line = element_blank(),
    strip.background = element_rect(fill = "white", color = "black"),
    strip.text = element_text(face = "bold"), axis.text.x = element_text(angle = 35, hjust = 1),
    plot.title = element_text(face = "bold"))

fc <- rbindlist(lapply(seq_len(nrow(wide)), function(i) {
  rbindlist(lapply(tx[-1], function(z) data.table(
    cell_line = wide$cell_line[i], ext_gene = wide$ext_gene[i], treatment = z,
    log2FC = wide[[paste0("LogFC.", z)]][i], FDR = wide[[paste0("FDR.", z)]][i]
  )))
}))
fc[, treatment := factor(treatment, levels = tx[-1])]
fc[, significant := FDR < 0.05]

p_fc <- ggplot(fc, aes(treatment, log2FC, color = treatment)) +
  geom_hline(yintercept = 0, color = "gray35", linewidth = 0.5) +
  geom_segment(aes(xend = treatment, y = 0, yend = log2FC), color = "gray70", linewidth = 0.75) +
  geom_point(aes(shape = significant), size = 3.2, stroke = 0.9) +
  facet_grid(ext_gene ~ cell_line) +
  scale_color_manual(values = tx_cols) +
  scale_shape_manual(values = c(`TRUE` = 16, `FALSE` = 1),
    labels = c(`TRUE` = "FDR < 0.05", `FALSE` = "Not significant")) +
  labs(title = "Treatment-associated transcript changes",
    subtitle = "Each treatment is compared with VEH within its cell line",
    x = NULL, y = "log2 fold change versus VEH", shape = NULL) +
  guides(color = "none") + theme_classic(base_size = 11) +
  theme(panel.border = element_rect(fill = NA, color = "black"), axis.line = element_blank(),
    strip.background = element_rect(fill = "white", color = "black"),
    strip.text = element_text(face = "bold"), axis.text.x = element_text(angle = 35, hjust = 1),
    legend.position = "top", plot.title = element_text(face = "bold"))

combined <- p_abs / p_fc + plot_annotation(tag_levels = "A")
ggsave(file.path(outdir, "EZH1_EZH2_TRANSCRIPTS_RKO_SW620.pdf"), combined,
  width = 9.5, height = 10, useDingbats = FALSE)
ggsave(file.path(outdir, "EZH1_EZH2_TRANSCRIPTS_RKO_SW620.png"), combined,
  width = 9.5, height = 10, dpi = 450, bg = "white")

fwrite(expr, file.path(outdir, "EZH1_EZH2_REPLICATE_EXPRESSION.tsv"), sep = "\t")
fwrite(means, file.path(outdir, "EZH1_EZH2_TREATMENT_MEANS.tsv"), sep = "\t")
fwrite(fc, file.path(outdir, "EZH1_EZH2_TREATMENT_VS_VEH_STATS.tsv"), sep = "\t")
print(means)
print(fc)
