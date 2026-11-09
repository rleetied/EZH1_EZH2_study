#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(DESeq2)
  library(ggplot2)
  library(pheatmap)
})

args <- commandArgs(trailingOnly = TRUE)

if (interactive() && length(args) < 1L) {
  message("Choose ATAC_consensus_featureCounts.txt in the file dialog.")
  count_file <- normalizePath(file.choose(), mustWork = TRUE)

  # The expected layout is canonical/counts/count_file.txt.
  # Place results in canonical/ATAC_QC.
  output_dir <- file.path(dirname(dirname(count_file)), "ATAC_QC")
} else {
  if (length(args) < 1L) {
    stop(
      "Usage: Rscript 01_ATAC_QC_PCA.R ",
      "/path/to/ATAC_consensus_featureCounts.txt [output_directory]"
    )
  }

  count_file <- normalizePath(args[[1]], mustWork = TRUE)
  output_dir <- if (length(args) >= 2L) args[[2]] else "ATAC_QC"
}

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

message("Reading featureCounts output: ", count_file)

fc <- read.delim(
  count_file,
  comment.char = "#",
  check.names = FALSE,
  stringsAsFactors = FALSE
)

required_annotation_columns <- c(
  "Geneid", "Chr", "Start", "End", "Strand", "Length"
)

if (!all(required_annotation_columns %in% colnames(fc))) {
  stop("The input does not have the expected featureCounts annotation columns.")
}

annotation <- fc[, required_annotation_columns, drop = FALSE]
count_columns <- setdiff(colnames(fc), required_annotation_columns)
counts <- as.matrix(fc[, count_columns, drop = FALSE])
storage.mode(counts) <- "integer"
rownames(counts) <- annotation$Geneid

sample_names <- basename(colnames(counts))
sample_names <- sub(
  "_PE\\.sorted\\.dedup\\.noMT\\.bt2\\.bam$",
  "",
  sample_names
)
colnames(counts) <- sample_names

expected_samples <- c(
  "A1RKOVeh", "A2RKOTAZ", "A3RKOVAL",
  "A4RKOD100", "A5RKOTD", "A6RKOVD",
  "B10RKOD100", "B11RKOTD", "B12RKOVD",
  "B7RKOVeh", "B8RKOTAZ", "B9RKOVAL",
  "C13RKOVeh", "C14RKOTAZ", "C15RKOVAL",
  "C16RKOD100", "C17RKOTD", "C18RKOVD"
)

if (!setequal(sample_names, expected_samples)) {
  stop(
    "Sample-name mismatch.\nObserved: ",
    paste(sample_names, collapse = ", "),
    "\nExpected: ",
    paste(expected_samples, collapse = ", ")
  )
}

get_condition <- function(sample) {
  if (grepl("Veh$", sample)) {
    return("Veh")
  }
  if (grepl("TAZ$", sample)) {
    return("TAZ")
  }
  if (grepl("VAL$", sample)) {
    return("VAL")
  }
  if (grepl("D100$", sample)) {
    return("DAC")
  }
  if (grepl("TD$", sample)) {
    return("TAZ_DAC")
  }
  if (grepl("VD$", sample)) {
    return("VAL_DAC")
  }
  stop("Could not assign condition for sample: ", sample)
}

sample_data <- data.frame(
  sample = sample_names,
  block = factor(substr(sample_names, 1L, 1L), levels = c("A", "B", "C")),
  condition = factor(
    vapply(sample_names, get_condition, character(1)),
    levels = c("Veh", "TAZ", "VAL", "DAC", "TAZ_DAC", "VAL_DAC")
  ),
  row.names = sample_names,
  stringsAsFactors = FALSE
)

if (anyNA(sample_data$block) || anyNA(sample_data$condition)) {
  stop("At least one sample has missing block or condition metadata.")
}

write.table(
  sample_data,
  file = file.path(output_dir, "sample_metadata.tsv"),
  sep = "\t",
  quote = FALSE,
  col.names = NA
)

write.table(
  annotation,
  file = file.path(output_dir, "consensus_peak_annotation.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

raw_library_sizes <- colSums(counts)

# Retain peaks with at least 10 fragments in at least 3 samples.
keep <- rowSums(counts >= 10L) >= 3L
filtered_counts <- counts[keep, , drop = FALSE]

filter_summary <- data.frame(
  metric = c(
    "consensus_peaks_input",
    "peaks_retained",
    "peaks_removed",
    "retained_percent"
  ),
  value = c(
    nrow(counts),
    nrow(filtered_counts),
    sum(!keep),
    100 * mean(keep)
  )
)

write.table(
  filter_summary,
  file = file.path(output_dir, "peak_filter_summary.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

message(
  "Retained ", nrow(filtered_counts), " of ", nrow(counts),
  " consensus peaks (", sprintf("%.2f", 100 * mean(keep)), "%)."
)

dds <- DESeqDataSetFromMatrix(
  countData = filtered_counts,
  colData = sample_data,
  design = ~ block + condition
)

dds <- DESeq(dds)
vsd <- vst(dds, blind = FALSE)

qc_table <- data.frame(
  sample = colnames(dds),
  block = colData(dds)$block,
  condition = colData(dds)$condition,
  raw_consensus_fragments = raw_library_sizes[colnames(dds)],
  filtered_consensus_fragments = colSums(counts(dds)),
  DESeq2_size_factor = sizeFactors(dds),
  row.names = NULL
)

write.table(
  qc_table,
  file = file.path(output_dir, "sample_QC_and_size_factors.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

saveRDS(dds, file.path(output_dir, "ATAC_dds_fitted.rds"))
saveRDS(vsd, file.path(output_dir, "ATAC_vst.rds"))

pca_data <- plotPCA(
  vsd,
  intgroup = c("condition", "block"),
  returnData = TRUE
)

percent_variance <- round(100 * attr(pca_data, "percentVar"), 1)
pca_data$sample <- rownames(pca_data)

write.table(
  pca_data,
  file = file.path(output_dir, "PCA_coordinates.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

pca_plot <- ggplot(
  pca_data,
  aes(
    x = PC1,
    y = PC2,
    color = condition,
    shape = block,
    label = sample
  )
) +
  geom_point(size = 4) +
  geom_text(vjust = -0.9, size = 3, show.legend = FALSE) +
  xlab(paste0("PC1: ", percent_variance[[1]], "% variance")) +
  ylab(paste0("PC2: ", percent_variance[[2]], "% variance")) +
  coord_equal() +
  theme_classic(base_size = 12) +
  theme(legend.position = "right")

ggsave(
  file.path(output_dir, "ATAC_PCA.pdf"),
  plot = pca_plot,
  width = 9,
  height = 7
)

ggsave(
  file.path(output_dir, "ATAC_PCA.png"),
  plot = pca_plot,
  width = 9,
  height = 7,
  dpi = 300
)

correlation_matrix <- cor(assay(vsd), method = "pearson")

write.table(
  correlation_matrix,
  file = file.path(output_dir, "VST_sample_correlations.tsv"),
  sep = "\t",
  quote = FALSE,
  col.names = NA
)

heatmap_annotation <- data.frame(
  condition = sample_data[colnames(correlation_matrix), "condition"],
  block = sample_data[colnames(correlation_matrix), "block"],
  row.names = colnames(correlation_matrix)
)

pdf(
  file.path(output_dir, "VST_sample_correlation_heatmap.pdf"),
  width = 10,
  height = 9
)
pheatmap(
  correlation_matrix,
  annotation_col = heatmap_annotation,
  annotation_row = heatmap_annotation,
  border_color = NA,
  main = "ATAC-seq VST sample correlations"
)
dev.off()

png(
  file.path(output_dir, "VST_sample_correlation_heatmap.png"),
  width = 3000,
  height = 2700,
  res = 300
)
pheatmap(
  correlation_matrix,
  annotation_col = heatmap_annotation,
  annotation_row = heatmap_annotation,
  border_color = NA,
  main = "ATAC-seq VST sample correlations"
)
dev.off()

library_plot_data <- qc_table

library_plot <- ggplot(
  library_plot_data,
  aes(
    x = reorder(sample, raw_consensus_fragments),
    y = raw_consensus_fragments,
    fill = condition
  )
) +
  geom_col() +
  coord_flip() +
  scale_y_continuous(labels = scales::label_comma()) +
  xlab(NULL) +
  ylab("Fragments assigned to consensus peaks") +
  theme_classic(base_size = 12)

ggsave(
  file.path(output_dir, "consensus_fragment_library_sizes.pdf"),
  plot = library_plot,
  width = 9,
  height = 7
)

ggsave(
  file.path(output_dir, "consensus_fragment_library_sizes.png"),
  plot = library_plot,
  width = 9,
  height = 7,
  dpi = 300
)

capture.output(
  sessionInfo(),
  file = file.path(output_dir, "R_sessionInfo.txt")
)

message("Finished successfully.")
message("Review these first:")
message("  ", file.path(output_dir, "ATAC_PCA.png"))
message("  ", file.path(output_dir, "VST_sample_correlation_heatmap.png"))
message("  ", file.path(output_dir, "sample_QC_and_size_factors.tsv"))
