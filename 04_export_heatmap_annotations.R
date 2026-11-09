#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(data.table))

root <- "RKO_TAZ_v_VAL_project/Final_figure_analysis/"
indir <- file.path(
  root, "outputs", "FINAL_UP_DOWN_rGREAT_RNA_PRIORITY_LINKS",
  "Combo_expression_scatter", "barrier_reclassification",
  "combined_K27me1_induced_bivalency"
)

input <- file.path(indir, "combined_barrier_reclassified_complete_data.tsv.gz")
output <- file.path(indir, "HEATMAP_gene_chromatin_barrier_annotations.tsv")

d <- fread(input)
stopifnot(!anyDuplicated(d$gene))

heatmap_annotations <- d[, .(
  Gene = gene,
  Chromatin_barrier = as.character(barrier_class),
  Regulatory_class = as.character(regulatory_class),
  Accessibility_class = fifelse(accessibility_panel == "Combo", "Comparable", "VD")
)]

barrier_order <- c(
  "Adaptive H3K27me1 barrier",
  "H3K27me1 + induced bivalency",
  "Induced bivalency",
  "H3K27me1 loss"
)
heatmap_annotations[, barrier_rank := match(Chromatin_barrier, barrier_order)]
heatmap_annotations[, regulatory_rank := match(Regulatory_class, c("TSS", "CRE"))]
heatmap_annotations[, accessibility_rank := match(Accessibility_class, c("VD", "Comparable"))]
setorder(heatmap_annotations, barrier_rank, regulatory_rank, accessibility_rank, Gene)
heatmap_annotations[, c("barrier_rank", "regulatory_rank", "accessibility_rank") := NULL]

fwrite(heatmap_annotations, output, sep = "\t")

cat("Saved:", output, "\n")
cat("Rows:", nrow(heatmap_annotations), "\n")
cat("Unique genes:", uniqueN(heatmap_annotations$Gene), "\n")
cat("Duplicated genes:", sum(duplicated(heatmap_annotations$Gene)), "\n")
