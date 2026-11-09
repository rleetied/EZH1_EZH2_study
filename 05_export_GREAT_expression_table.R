#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(data.table))

root <- "RKO_TAZ_v_VAL_project/Final_figure_analysis/"
indir <- file.path(root, "outputs", "FINAL_UP_DOWN_rGREAT_RNA_PRIORITY_LINKS")
outdir <- file.path(indir, "FINAL_GREAT_GENE_EXPRESSION_TABLE")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

d <- fread(file.path(indir,
  "ONE_representative_up_or_down_peak_per_gene_PRIORITY_RESOLVED.tsv"))

d[, accessibility_class := fcase(
  accessibility_class == "VD-favored", "VD-favored",
  accessibility_class == "Comparable TD/VD", "Comparable",
  default = accessibility_class
)]

d[, chromatin_barrier := fcase(
  H3K27me1_state == "Loss", "H3K27me1 loss",
  H3K27me1_state == "Gain" & accessibility_class == "VD-favored",
    "Adaptive H3K27me1 barrier",
  H3K27me1_state == "Gain" & accessibility_class == "Comparable",
    "H3K27me1 + induced bivalency",
  H3K27me1_state == "Maintain" & accessibility_class == "VD-favored",
    "H3K27me1 + induced bivalency",
  H3K27me1_state == "Maintain" & accessibility_class == "Comparable",
    "Induced bivalency",
  default = NA_character_
)]

d[, `:=`(
  TD_expression_displayed = !is.na(LogFC.TD) &
    (LogFC.TD <= -0.5 | LogFC.TD >= 1),
  VD_expression_displayed = !is.na(LogFC.VD) &
    (LogFC.VD <= -0.5 | LogFC.VD >= 1),
  TD_expression_significant = !is.na(FDR.TD) & FDR.TD < 0.05 &
    (LogFC.TD <= -0.5 | LogFC.TD >= 1),
  VD_expression_significant = !is.na(FDR.VD) & FDR.VD < 0.05 &
    (LogFC.VD <= -0.5 | LogFC.VD >= 1)
)]

setnames(d, c("LogFC.TD", "FDR.TD", "logCPM.TD",
              "LogFC.VD", "FDR.VD", "logCPM.VD"),
  c("TD_RNA_log2FC", "TD_RNA_FDR", "TD_RNA_logCPM",
    "VD_RNA_log2FC", "VD_RNA_FDR", "VD_RNA_logCPM"))

keep <- c(
  "gene", "ENSG", "ens_gene", "gene_biotype",
  "peak_id", "chr", "start", "end", "distance_to_TSS",
  "absolute_distance_bp", "regulatory_class", "region_set",
  "accessibility_direction", "accessibility_class", "H3K27me1_state",
  "chromatin_barrier",
  "TAZ_vs_Veh_log2FC", "TAZ_vs_Veh_padj",
  "VAL_vs_Veh_log2FC", "VAL_vs_Veh_padj",
  "DAC_vs_Veh_log2FC", "DAC_vs_Veh_padj",
  "TAZ_DAC_vs_Veh_log2FC", "TAZ_DAC_vs_Veh_padj",
  "VAL_DAC_vs_Veh_log2FC", "VAL_DAC_vs_Veh_padj",
  "TD_RNA_log2FC", "TD_RNA_FDR", "TD_RNA_logCPM",
  "VD_RNA_log2FC", "VD_RNA_FDR", "VD_RNA_logCPM",
  "TD_expression_displayed", "TD_expression_significant",
  "VD_expression_displayed", "VD_expression_significant",
  "RNA_match_status"
)
d <- d[, ..keep]
setorder(d, regulatory_class, accessibility_class, chromatin_barrier, gene)

fwrite(d, file.path(outdir,
  "FINAL_GREAT_gene_peak_ATAC_RNA_chromatin_barrier_table.tsv"), sep = "\t")
fwrite(d[RNA_match_status == "unique_RNA_match"], file.path(outdir,
  "FINAL_GREAT_gene_peak_ATAC_RNA_unique_RNA_matches.tsv"), sep = "\t")
fwrite(d[TD_expression_displayed | VD_expression_displayed], file.path(outdir,
  "FINAL_GREAT_gene_peak_expression_displayed_in_either_plot.tsv"), sep = "\t")

# Concise table supporting the final VD-centered scatter plot.
vd <- d[RNA_match_status == "unique_RNA_match" & VD_expression_displayed == TRUE,
  .(
    gene, ENSG, gene_biotype,
    peak_id, chr, start, end,
    regulatory_class, distance_to_TSS, absolute_distance_bp,
    accessibility_direction, accessibility_class,
    H3K27me1_state, chromatin_barrier,
    VD_ATAC_log2FC = VAL_DAC_vs_Veh_log2FC,
    VD_ATAC_padj = VAL_DAC_vs_Veh_padj,
    VD_RNA_log2FC, VD_RNA_FDR, VD_RNA_logCPM,
    VD_expression_significant
  )]
setorder(vd, regulatory_class, accessibility_class, chromatin_barrier, gene)
stopifnot(!anyDuplicated(vd$gene))
fwrite(vd, file.path(outdir,
  "FINAL_VD_GREAT_gene_expression_link_table.tsv"), sep = "\t")

cat("All representative GREAT-linked genes:", nrow(d), "\n")
cat("Unique RNA matches:", sum(d$RNA_match_status == "unique_RNA_match"), "\n")
cat("Displayed under TD or VD expression thresholds:",
  sum(d$TD_expression_displayed | d$VD_expression_displayed, na.rm = TRUE), "\n")
cat("Genes displayed in final VD plot:", nrow(vd), "\n")
cat("Saved to:", outdir, "\n")
