#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(DESeq2)
  library(ggplot2)
  library(pheatmap)
})

# Statistical thresholds used for classification and summary tables.
alpha <- 0.05
lfc_threshold <- 1
top_n_heatmap <- 50L

args <- commandArgs(trailingOnly = TRUE)

if (interactive() && length(args) < 1L) {
  message("Choose ATAC_dds_fitted.rds in the file dialog.")
  dds_file <- normalizePath(file.choose(), mustWork = TRUE)
  qc_dir <- dirname(dds_file)
  output_dir <- file.path(dirname(qc_dir), "ATAC_differential")
} else {
  if (length(args) < 1L) {
    stop(
      "Usage: Rscript 02_ATAC_differential_accessibility.R ",
      "/path/to/ATAC_dds_fitted.rds [output_directory]"
    )
  }

  dds_file <- normalizePath(args[[1]], mustWork = TRUE)
  qc_dir <- dirname(dds_file)
  output_dir <- if (length(args) >= 2L) {
    args[[2]]
  } else {
    file.path(dirname(qc_dir), "ATAC_differential")
  }
}

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(output_dir, "tables"), showWarnings = FALSE)
dir.create(file.path(output_dir, "plots"), showWarnings = FALSE)
dir.create(file.path(output_dir, "beds"), showWarnings = FALSE)
dir.create(file.path(output_dir, "heatmaps"), showWarnings = FALSE)

write_gzip_table <- function(x, file) {
  connection <- gzfile(file, open = "wt")
  on.exit(close(connection), add = TRUE)

  write.table(
    x,
    file = connection,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
}

message("Loading fitted DESeq2 object: ", dds_file)
dds <- readRDS(dds_file)

if (!inherits(dds, "DESeqDataSet")) {
  stop("Selected RDS file is not a DESeqDataSet.")
}

required_conditions <- c(
  "Veh", "TAZ", "VAL", "DAC", "TAZ_DAC", "VAL_DAC"
)

observed_conditions <- levels(colData(dds)$condition)

if (!all(required_conditions %in% observed_conditions)) {
  stop(
    "Required condition levels are missing. Observed levels: ",
    paste(observed_conditions, collapse = ", ")
  )
}

annotation_file <- file.path(qc_dir, "consensus_peak_annotation.tsv")

if (!file.exists(annotation_file)) {
  stop("Could not find peak annotation beside the DDS: ", annotation_file)
}

annotation <- read.delim(
  annotation_file,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

annotation <- annotation[
  match(rownames(dds), annotation$Geneid),
  ,
  drop = FALSE
]

if (anyNA(annotation$Geneid) || !identical(annotation$Geneid, rownames(dds))) {
  stop("Peak annotation could not be matched cleanly to the DESeq2 rows.")
}

annotation$bed_start <- annotation$Start - 1L
annotation$bed_end <- annotation$End

message("Model: ", paste(as.character(design(dds)), collapse = " "))
message("Coefficients:")
message(paste("  ", resultsNames(dds), collapse = "\n"))

# Save normalized counts once for downstream browsing or custom plotting.
normalized_counts <- as.data.frame(counts(dds, normalized = TRUE))
normalized_counts <- cbind(
  peak_id = rownames(normalized_counts),
  normalized_counts
)

write_gzip_table(
  normalized_counts,
  file = file.path(
    output_dir,
    "tables",
    "DESeq2_normalized_counts.tsv.gz"
  )
)

# Use the already-created VST object when available.
vst_file <- file.path(qc_dir, "ATAC_vst.rds")
vsd <- if (file.exists(vst_file)) {
  readRDS(vst_file)
} else {
  vst(dds, blind = FALSE)
}

sample_annotation <- as.data.frame(
  colData(dds)[, c("condition", "block"), drop = FALSE]
)

condition_colors <- c(
  Veh = "#F8766D",
  TAZ = "#C49A00",
  VAL = "#00BA38",
  DAC = "#00BFC4",
  TAZ_DAC = "#619CFF",
  VAL_DAC = "#F564E3"
)

contrast_specs <- list(
  TAZ_vs_Veh = list(
    kind = "simple",
    contrast = c("condition", "TAZ", "Veh")
  ),
  VAL_vs_Veh = list(
    kind = "simple",
    contrast = c("condition", "VAL", "Veh")
  ),
  DAC_vs_Veh = list(
    kind = "simple",
    contrast = c("condition", "DAC", "Veh")
  ),
  TAZ_DAC_vs_Veh = list(
    kind = "simple",
    contrast = c("condition", "TAZ_DAC", "Veh")
  ),
  VAL_DAC_vs_Veh = list(
    kind = "simple",
    contrast = c("condition", "VAL_DAC", "Veh")
  ),
  TAZ_DAC_vs_TAZ = list(
    kind = "simple",
    contrast = c("condition", "TAZ_DAC", "TAZ")
  ),
  VAL_DAC_vs_VAL = list(
    kind = "simple",
    contrast = c("condition", "VAL_DAC", "VAL")
  ),
  TAZ_DAC_vs_DAC = list(
    kind = "simple",
    contrast = c("condition", "TAZ_DAC", "DAC")
  ),
  VAL_DAC_vs_DAC = list(
    kind = "simple",
    contrast = c("condition", "VAL_DAC", "DAC")
  ),
  TAZ_DAC_nonadditivity = list(
    kind = "linear",
    contrast = list(
      "condition_TAZ_DAC_vs_Veh",
      c("condition_TAZ_vs_Veh", "condition_DAC_vs_Veh")
    )
  ),
  VAL_DAC_nonadditivity = list(
    kind = "linear",
    contrast = list(
      "condition_VAL_DAC_vs_Veh",
      c("condition_VAL_vs_Veh", "condition_DAC_vs_Veh")
    )
  )
)

make_result <- function(spec) {
  if (identical(spec$kind, "simple")) {
    return(
      results(
        dds,
        contrast = spec$contrast,
        alpha = alpha,
        independentFiltering = TRUE
      )
    )
  }

  results(
    dds,
    contrast = spec$contrast,
    alpha = alpha,
    independentFiltering = TRUE
  )
}

classify_results <- function(result_table) {
  result_table$classification <- "Not_significant"

  result_table$classification[
    !is.na(result_table$padj) &
      result_table$padj < alpha &
      result_table$log2FoldChange >= lfc_threshold
  ] <- "Increased"

  result_table$classification[
    !is.na(result_table$padj) &
      result_table$padj < alpha &
      result_table$log2FoldChange <= -lfc_threshold
  ] <- "Decreased"

  result_table$classification <- factor(
    result_table$classification,
    levels = c("Not_significant", "Increased", "Decreased")
  )

  result_table
}

make_ma_plot <- function(result_table, contrast_name) {
  plot_data <- result_table
  plot_data$baseMean_for_plot <- pmax(plot_data$baseMean, 0.1)

  ggplot(
    plot_data,
    aes(
      x = baseMean_for_plot,
      y = log2FoldChange,
      color = classification
    )
  ) +
    geom_point(alpha = 0.35, size = 0.7) +
    geom_hline(
      yintercept = c(-lfc_threshold, lfc_threshold),
      linetype = "dashed",
      color = "grey40"
    ) +
    scale_x_log10(labels = scales::label_comma()) +
    scale_color_manual(
      values = c(
        Not_significant = "grey70",
        Increased = "#C51B7D",
        Decreased = "#2C7BB6"
      )
    ) +
    xlab("Mean normalized fragment count") +
    ylab("log2 fold change") +
    ggtitle(contrast_name) +
    theme_classic(base_size = 12)
}

make_volcano_plot <- function(result_table, contrast_name) {
  plot_data <- result_table
  plot_data$minus_log10_padj <- -log10(
    pmax(plot_data$padj, .Machine$double.xmin)
  )
  plot_data$minus_log10_padj[is.na(plot_data$padj)] <- NA_real_

  ggplot(
    plot_data,
    aes(
      x = log2FoldChange,
      y = minus_log10_padj,
      color = classification
    )
  ) +
    geom_point(alpha = 0.4, size = 0.7, na.rm = TRUE) +
    geom_vline(
      xintercept = c(-lfc_threshold, lfc_threshold),
      linetype = "dashed",
      color = "grey40"
    ) +
    geom_hline(
      yintercept = -log10(alpha),
      linetype = "dashed",
      color = "grey40"
    ) +
    scale_color_manual(
      values = c(
        Not_significant = "grey70",
        Increased = "#C51B7D",
        Decreased = "#2C7BB6"
      )
    ) +
    xlab("log2 fold change") +
    ylab(expression(-log[10]("adjusted p-value"))) +
    ggtitle(contrast_name) +
    theme_classic(base_size = 12)
}

make_top_peak_heatmap <- function(result_table, contrast_name) {
  eligible <- which(
    !is.na(result_table$padj) &
      result_table$padj < alpha
  )

  if (length(eligible) < 2L) {
    message("Skipping heatmap for ", contrast_name, ": fewer than 2 peaks.")
    return(invisible(NULL))
  }

  ordered <- eligible[
    order(result_table$padj[eligible], na.last = NA)
  ]
  selected <- head(ordered, top_n_heatmap)
  peak_ids <- result_table$peak_id[selected]

  heatmap_matrix <- assay(vsd)[peak_ids, , drop = FALSE]
  heatmap_matrix <- t(scale(t(heatmap_matrix)))
  heatmap_matrix[!is.finite(heatmap_matrix)] <- 0

  pdf(
    file.path(
      output_dir,
      "heatmaps",
      paste0(contrast_name, "_top", length(selected), "_heatmap.pdf")
    ),
    width = 10,
    height = 10
  )
  pheatmap(
    heatmap_matrix,
    annotation_col = sample_annotation,
    annotation_colors = list(condition = condition_colors),
    show_rownames = FALSE,
    border_color = NA,
    clustering_distance_cols = "correlation",
    main = paste0(contrast_name, ": top ", length(selected), " peaks")
  )
  dev.off()
}

result_list <- list()
summary_rows <- list()

for (contrast_name in names(contrast_specs)) {
  message("Processing contrast: ", contrast_name)

  res <- make_result(contrast_specs[[contrast_name]])
  result_table <- as.data.frame(res)
  result_table$peak_id <- rownames(result_table)

  result_table <- cbind(
    annotation[
      match(result_table$peak_id, annotation$Geneid),
      c(
        "Geneid", "Chr", "Start", "End", "Strand",
        "Length", "bed_start", "bed_end"
      ),
      drop = FALSE
    ],
    result_table[
      ,
      c(
        "baseMean", "log2FoldChange", "lfcSE",
        "stat", "pvalue", "padj", "peak_id"
      ),
      drop = FALSE
    ]
  )

  result_table <- classify_results(result_table)
  result_table <- result_table[
    order(result_table$padj, -abs(result_table$log2FoldChange), na.last = TRUE),
    ,
    drop = FALSE
  ]

  result_list[[contrast_name]] <- result_table

  write_gzip_table(
    result_table,
    file = file.path(
      output_dir,
      "tables",
      paste0(contrast_name, "_all_peaks.tsv.gz")
    )
  )

  significant <- result_table[
    result_table$classification != "Not_significant",
    ,
    drop = FALSE
  ]

  write.table(
    significant,
    file = file.path(
      output_dir,
      "tables",
      paste0(contrast_name, "_significant_peaks.tsv")
    ),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )

  increased_bed <- significant[
    significant$classification == "Increased",
    c("Chr", "bed_start", "bed_end", "Geneid", "log2FoldChange", "padj"),
    drop = FALSE
  ]

  decreased_bed <- significant[
    significant$classification == "Decreased",
    c("Chr", "bed_start", "bed_end", "Geneid", "log2FoldChange", "padj"),
    drop = FALSE
  ]

  write.table(
    increased_bed,
    file = file.path(
      output_dir,
      "beds",
      paste0(contrast_name, "_increased.bed")
    ),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    col.names = FALSE
  )

  write.table(
    decreased_bed,
    file = file.path(
      output_dir,
      "beds",
      paste0(contrast_name, "_decreased.bed")
    ),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    col.names = FALSE
  )

  increased_count <- sum(
    result_table$classification == "Increased",
    na.rm = TRUE
  )
  decreased_count <- sum(
    result_table$classification == "Decreased",
    na.rm = TRUE
  )

  summary_rows[[contrast_name]] <- data.frame(
    contrast = contrast_name,
    total_peaks = nrow(result_table),
    peaks_with_padj = sum(!is.na(result_table$padj)),
    increased = increased_count,
    decreased = decreased_count,
    significant_total = increased_count + decreased_count,
    alpha = alpha,
    absolute_log2FC_threshold = lfc_threshold,
    stringsAsFactors = FALSE
  )

  ma_plot <- make_ma_plot(result_table, contrast_name)
  volcano_plot <- make_volcano_plot(result_table, contrast_name)

  ggsave(
    file.path(output_dir, "plots", paste0(contrast_name, "_MA.pdf")),
    ma_plot,
    width = 7,
    height = 6
  )

  ggsave(
    file.path(output_dir, "plots", paste0(contrast_name, "_MA.png")),
    ma_plot,
    width = 7,
    height = 6,
    dpi = 300
  )

  ggsave(
    file.path(output_dir, "plots", paste0(contrast_name, "_volcano.pdf")),
    volcano_plot,
    width = 7,
    height = 6
  )

  ggsave(
    file.path(output_dir, "plots", paste0(contrast_name, "_volcano.png")),
    volcano_plot,
    width = 7,
    height = 6,
    dpi = 300
  )

  make_top_peak_heatmap(result_table, contrast_name)
}

contrast_summary <- do.call(rbind, summary_rows)
rownames(contrast_summary) <- NULL

write.table(
  contrast_summary,
  file = file.path(output_dir, "differential_peak_summary.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

