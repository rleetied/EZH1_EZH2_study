#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(GenomicRanges)
  library(ggplot2)
})

root <- "RKO_TAZ_v_VAL_project/Final_figure_analysis/"
master <- "RKO_TAZ_v_VAL_project/20260731_final_called_regions"
relax_td_p01 <- identical(Sys.getenv("RELAX_TD_P01"), "1")
combine_vd_favored <- identical(Sys.getenv("COMBINE_VD_FAVORED"), "1")
out_name <- if (combine_vd_favored) {
  "TD_VD_COMPARABLE_VS_VD_FAVORED_H3K27me1_DISTANCE"
} else if (relax_td_p01) {
  "TD_VD_MAGNITUDE_CLASSES_TD_PADJ01_H3K27me1_DISTANCE"
} else {
  "TD_VD_MAGNITUDE_CLASSES_H3K27me1_DISTANCE"
}
out <- file.path(root, "outputs", out_name)
dir.create(out, recursive = TRUE, showWarnings = FALSE)

shared_file <- file.path(root, "outputs", "TD_VD_ATAC_GAIN_VENN_SHARED_20PCT",
                         "TD_VD_shared_ATAC_gain_peak_level_20pct_classification.tsv")
membership_file <- file.path(root, "outputs", "ATAC_treatment_overlaps", "increased",
                             "all_increased_peaks_membership.tsv")
full_stats_file <- file.path(root, "outputs", "ATAC_VAL_DAC_specificity_20pct", "tables",
                             "VAL_DAC_primary_peaks_full_classification.tsv")
me1_paths <- c(
  TAZ_gain = file.path(master, "H3K27me1",
    "TAZ_called.GRanges_H3K27me1_gain.var_1.logFC_0.5.TAZ.bed"),
  VAL_only_loss = file.path(master, "H3K27me1",
    "TAZ_called.GRanges_H3K27me1_lost.var_1.logFC_0.5.VAL.bed"),
  Shared_TAZ_VAL_loss = file.path(master, "H3K27me1",
    "VEH_called.GRanges_H3K27me1_lost.var_1.logFC_0.5.TAZ.bed")
)
stopifnot(file.exists(shared_file), file.exists(membership_file),
          file.exists(full_stats_file), all(file.exists(me1_paths)))

read_bed3 <- function(path) {
  x <- fread(path, header = FALSE, select = 1:3)
  setnames(x, c("chr", "start", "end"))
  x[, `:=`(chr = as.character(chr), start = as.integer(start), end = as.integer(end))]
  x
}

class_levels <- c("Comparable shared", "Shared; VD >=20% stronger", "VD only")
if (relax_td_p01) {
  # Start from all primary VD increases (padj < 0.05, log2FC >= 1). A TD
  # response is detectable when it is positive and has padj < 0.10; those
  # loci are then split by the literal VD >= 1.2 * TD magnitude rule.
  z <- fread(full_stats_file)[primary_direction == "Increased"]
  z[, td_detected := !is.na(TAZ_DAC_vs_Veh_padj) &
       TAZ_DAC_vs_Veh_padj < 0.10 & TAZ_DAC_vs_Veh_log2FC > 0]
  z[, class := fifelse(
    !td_detected, "VD only",
    fifelse(VAL_DAC_vs_Veh_log2FC >= 1.20 * TAZ_DAC_vs_Veh_log2FC,
            "Shared; VD >=20% stronger", "Comparable shared")
  )]
  d <- z[, .(peak_id, chr = Chr, start = bed_start, end = bed_end, class)]
  expected <- c(6101L, 7093L, 2480L)
} else {
  shared <- fread(shared_file)
  shared[, class := fifelse(
    magnitude_class == "VD at least 20% greater",
    "Shared; VD >=20% stronger",
    "Comparable shared"
  )]
  shared <- shared[, .(peak_id, chr = Chr, start = bed_start, end = bed_end, class)]

  membership <- fread(membership_file)
  vd_only <- membership[TAZ_DAC == 0L & VAL_DAC == 1L,
                        .(peak_id, chr, start, end, class = "VD only")]
  d <- rbindlist(list(shared, vd_only), use.names = TRUE)
  expected <- c(5654L, 4239L, 5781L)
}
observed <- table(factor(d$class, levels = class_levels))
stopifnot(all(as.integer(observed) == expected), !anyDuplicated(d$peak_id))

if (combine_vd_favored) {
  d[class != "Comparable shared", class := "VD-favored"]
  class_levels <- c("Comparable shared", "VD-favored")
  expected <- c(5654L, 10020L)
  observed <- table(factor(d$class, levels = class_levels))
  stopifnot(all(as.integer(observed) == expected))
}

bed_dir <- file.path(out, "BEDs")
dir.create(bed_dir, recursive = TRUE, showWarnings = FALSE)
for (cl in class_levels) {
  fwrite(d[class == cl, .(chr, start, end, peak_id)],
         file.path(bed_dir, paste0(gsub("[^A-Za-z0-9]+", "_", cl), ".bed")),
         sep = "\t", col.names = FALSE)
}

me1 <- lapply(me1_paths, read_bed3)
distance_vector <- function(query, subject) {
  q <- GRanges(query$chr, IRanges(query$start + 1L, query$end))
  s <- GRanges(subject$chr, IRanges(subject$start, subject$end))
  h <- distanceToNearest(q, s, ignore.strand = TRUE)
  ans <- rep(NA_integer_, length(q))
  ans[queryHits(h)] <- mcols(h)$distance
  ans
}

windows <- c(seq(0L, 1000L, 100L), seq(2000L, 10000L, 1000L))
curve <- rbindlist(lapply(class_levels, function(cl) {
  q <- d[class == cl]
  rbindlist(lapply(names(me1), function(st) {
    dist <- distance_vector(q, me1[[st]])
    rbindlist(lapply(windows, function(w) data.table(
      class = cl, state = st, window_bp = w,
      near = sum(dist <= w, na.rm = TRUE), total = nrow(q),
      percent = 100 * mean(dist <= w, na.rm = TRUE)
    )))
  }))
}))

state_levels <- c("TAZ_gain", "VAL_only_loss", "Shared_TAZ_VAL_loss")
state_labels <- c(
  TAZ_gain = "TAZ-specific H3K27me1 gain",
  VAL_only_loss = "TAZ-maintained / VAL-only H3K27me1 loss",
  Shared_TAZ_VAL_loss = "Shared TAZ/VAL H3K27me1 loss"
)
curve[, `:=`(
  class = factor(class, levels = class_levels),
  state = factor(state, levels = state_levels),
  H3K27me1_state = factor(unname(state_labels[as.character(state)]),
                          levels = unname(state_labels[state_levels]))
)]
fwrite(curve, file.path(out, "TD_VD_magnitude_classes_H3K27me1_cumulative_distance.tsv"), sep = "\t")

endpoint <- curve[window_bp == 10000L]
endpoint[, `:=`(
  raw_label = scales::comma(near),
  pct_label = sprintf("%.1f%%", percent)
)]

comparisons <- if (combine_vd_favored) {
  list(c("VD-favored", "Comparable shared"))
} else {
  list(
    c("Shared; VD >=20% stronger", "Comparable shared"),
    c("VD only", "Comparable shared"),
    c("VD only", "Shared; VD >=20% stronger")
  )
}
stats <- rbindlist(lapply(state_levels, function(st) rbindlist(lapply(comparisons, function(cp) {
  a <- endpoint[state == st & class == cp[1]]
  b <- endpoint[state == st & class == cp[2]]
  ft <- fisher.test(matrix(c(a$near, a$total-a$near, b$near, b$total-b$near),
                           nrow = 2, byrow = TRUE))
  data.table(
    state = st, group_1 = cp[1], group_2 = cp[2],
    group_1_near = a$near, group_1_total = a$total,
    group_2_near = b$near, group_2_total = b$total,
    OR_group1_vs_group2 = unname(ft$estimate),
    CI_lower = unname(ft$conf.int[1]), CI_upper = unname(ft$conf.int[2]),
    log2OR = log2(unname(ft$estimate)), p = ft$p.value
  )
}))))
fwrite(stats, file.path(out, "TD_VD_magnitude_classes_H3K27me1_10kb_Fisher_tests.tsv"), sep = "\t")

cols <- c(
  "Comparable shared" = "#931100",
  "Shared; VD >=20% stronger" = "#B23AB8",
  "VD only" = "#ff84ff",
  "VD-favored" = "#D682FF"
)

base_theme <- theme_bw(base_size = 11) + theme(
  panel.border = element_rect(color = "black", fill = NA, linewidth = 0.9),
  panel.grid.major = element_line(color = "gray92", linewidth = 0.3),
  panel.grid.minor = element_blank(),
  strip.background = element_rect(fill = "white", color = "black", linewidth = 0.9),
  strip.text = element_text(face = "bold"),
  axis.text = element_text(color = "black"),
  legend.position = "top"
)

make_plot <- function(yvar, ylabel, suffix, label_col) {
  class_n <- setNames(as.integer(observed), class_levels)
  legend_labels <- setNames(
    paste0(class_levels, " (n = ", scales::comma(class_n), ")"),
    class_levels
  )
  subtitle_text <- if (relax_td_p01) {
    "TD response defined by adjusted P < 0.10 and positive log2FC; detectable TD responses partitioned by VD >= 1.2 x TD"
  } else if (combine_vd_favored) {
    "VD-favored combines VD-only peaks with shared peaks whose VD accessibility is at least 20% stronger than TD"
  } else {
    "Strict overlap partitioned by whether VD accessibility is at least 20% stronger than TD"
  }
  p <- ggplot(curve, aes(window_bp / 1000, .data[[yvar]], color = class, group = class)) +
    geom_line(linewidth = 1.05) + geom_point(size = 1.8) +
    geom_text(data = endpoint, aes(label = .data[[label_col]]),
              hjust = -0.10, size = 3.0, fontface = "bold", show.legend = FALSE) +
    facet_wrap(~H3K27me1_state, nrow = 1, scales = "fixed") +
    scale_color_manual(values = cols, labels = legend_labels) +
    scale_x_continuous(breaks = c(0, 1, 3, 5, 7, 10), limits = c(0, 11.2),
                       expand = expansion(mult = c(0, 0))) +
    scale_y_continuous(labels = if (yvar == "percent") function(x) paste0(x, "%") else scales::label_comma(),
                       expand = expansion(mult = c(0, 0.09))) +
    labs(
      title = "TD/VD increased-accessibility classes near H3K27me1 dynamics",
      subtitle = subtitle_text,
      x = expression("Accessibility distance to "*Delta*"H3K27me1 (kb)"),
      y = ylabel, color = NULL
    ) + base_theme +
    theme(plot.title = element_text(face = "bold", size = 14),
          plot.subtitle = element_text(color = "gray30"))
  ggsave(file.path(out, paste0("TD_VD_MAGNITUDE_CLASSES_H3K27me1_DISTANCE_", suffix, ".png")),
         p, width = 14, height = 5.8, dpi = 600, bg = "white")
  ggsave(file.path(out, paste0("TD_VD_MAGNITUDE_CLASSES_H3K27me1_DISTANCE_", suffix, ".pdf")),
         p, width = 14, height = 5.8, device = pdf)
}

make_plot("near", "# Accessibility increases", "RAW_COUNTS", "raw_label")
make_plot("percent", "% Accessibility increases", "PERCENT", "pct_label")
