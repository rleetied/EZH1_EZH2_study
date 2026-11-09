suppressPackageStartupMessages({library(data.table); library(ggplot2)})

root <- "RKO_TAZ_v_VAL_project/Final_figure_analysis/"
outdir <- file.path(root, "outputs", "FIG2_EZH1_TRANSCRIPTIONAL_DRAFT")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
cols <- c(RKO = "black", SW620 = "gray60")

read_one <- function(line) {
  x <- fread(file.path(root, "outputs", paste0(line, "_VD_vs_TD_MSIGDB_ENRICHMENT"),
    "PRERANKED_GSEA_ALL_REQUESTED_COLLECTIONS.tsv"))
  x[, cell_line := line]
  x
}
x <- rbindlist(lapply(c("RKO", "SW620"), read_one), fill = TRUE)

sel <- data.table(
  pathway = c(
    "BENPORATH_ES_WITH_H3K27ME3", "BENPORATH_PRC2_TARGETS",
    "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION", "HALLMARK_COAGULATION",
    "HALLMARK_INFLAMMATORY_RESPONSE", "HALLMARK_KRAS_SIGNALING_UP",
    "NABA_CORE_MATRISOME", "REACTOME_EXTRACELLULAR_MATRIX_ORGANIZATION",
    "REACTOME_COLLAGEN_CHAIN_TRIMERIZATION",
    "REACTOME_NON_INTEGRIN_MEMBRANE_ECM_INTERACTIONS",
    "KEGG_CALCIUM_SIGNALING_PATHWAY"
  ),
  label = c(
    "H3K27me3-associated ES genes", "PRC2 targets",
    "EMT", "Coagulation", "Inflammatory response", "KRAS signaling up",
    "Core matrisome", "Extracellular matrix organization",
    "Collagen-chain trimerization", "Non-integrin membrane–ECM interactions",
    "Calcium signaling"
  ),
  group = c(rep("PRC2", 2), rep("Hallmark", 4), rep("ECM", 4), "KEGG")
)
d <- merge(x, sel, by = "pathway")
d[, significant := !is.na(padj) & padj < 0.05]
d[, group := factor(group, levels = c("PRC2", "Hallmark", "ECM", "KEGG"))]
# Preserve intentional biological ordering within each facet.
d[, label := factor(label, levels = rev(sel$label))]

p <- ggplot(d, aes(NES, label, color = cell_line, group = label)) +
  geom_vline(xintercept = 0, color = "gray30", linewidth = 0.5) +
  geom_line(color = "gray75", linewidth = 0.75) +
  geom_point(aes(shape = significant), size = 3.5, stroke = 0.9) +
  scale_color_manual(values = cols) +
  scale_shape_manual(values = c(`TRUE` = 16, `FALSE` = 1),
    labels = c(`TRUE` = "FDR < 0.05", `FALSE` = "Not significant")) +
  coord_cartesian(xlim = c(-0.75, 4.6)) +
  labs(title = "VD-versus-TD pathway separation in RKO and SW620",
    subtitle = "Positive NES favors VD; pathways are ordered by biological theme",
    x = "Normalized enrichment score", y = NULL, color = NULL, shape = NULL) +
  theme_classic(base_size = 11) +
  theme(panel.border = element_rect(fill = NA, color = "black", linewidth = 0.65),
    axis.line = element_blank(), legend.position = "top",
    plot.title = element_text(face = "bold"), axis.text.y = element_text(color = "black"))

ggsave(file.path(outdir, "FIG2_SELECTED_PATHWAY_DUMBBELL.pdf"), p,
  width = 8.2, height = 8.4, useDingbats = FALSE)
ggsave(file.path(outdir, "FIG2_SELECTED_PATHWAY_DUMBBELL.png"), p,
  width = 8.2, height = 8.4, dpi = 450, bg = "white")
fwrite(d[, .(group, pathway, label, cell_line, NES, padj, significant)],
  file.path(outdir, "FIG2_SELECTED_PATHWAY_DUMBBELL_DATA.tsv"), sep = "\t")
