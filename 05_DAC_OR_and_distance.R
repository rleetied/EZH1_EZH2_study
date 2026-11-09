#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(GenomicRanges)
  library(ggplot2)
  library(patchwork)
})

root <- "RKO_TAZ_v_VAL_project/Final_figure_analysis/"
master <- "RKO_TAZ_v_VAL_project/20260731_final_called_regions"
out <- file.path(root, "outputs", "DAC_TD_VD_COMPARABLE_VS_VD_FAVORED_H3K27me1_OR_DISTANCE")
dir.create(out, recursive = TRUE, showWarnings = FALSE)

read_bed <- function(path) {
  x <- fread(path, header = FALSE)
  x <- x[, 1:min(4L, ncol(x)), with = FALSE]
  setnames(x, c("chr", "start", "end", if (ncol(x) == 4L) "peak_id"))
  if (ncol(x) == 3L) x[, peak_id := paste(chr, start, end, sep = ":")]
  x[, `:=`(chr = as.character(chr), start = as.integer(start), end = as.integer(end))]
  x
}

bed_dir <- file.path(root, "outputs", "TD_VD_COMPARABLE_VS_VD_FAVORED_H3K27me1_DISTANCE", "BEDs")
groups <- list(
  DAC = fread(file.path(root, "outputs", "ATAC_treatment_overlaps", "increased",
                        "all_increased_peaks_membership.tsv"))[DAC == 1L,
      .(chr, start, end, peak_id)],
  `Comparable shared` = read_bed(file.path(bed_dir, "Comparable_shared.bed")),
  `VD-favored` = read_bed(file.path(bed_dir, "VD_favored.bed"))
)
expected <- c(DAC = 3711L, `Comparable shared` = 5654L, `VD-favored` = 10020L)
stopifnot(all(vapply(groups, nrow, integer(1)) == expected))

# Export only the two response-class BEDs requested for downstream work.
bed_out <- file.path(out, "BEDs")
dir.create(bed_out, recursive = TRUE, showWarnings = FALSE)
fwrite(groups[["Comparable shared"]][, .(chr, start, end, peak_id)],
       file.path(bed_out, "Comparable_TD_VD_accessibility_gain.bed"),
       sep = "\t", col.names = FALSE)
fwrite(groups[["VD-favored"]][, .(chr, start, end, peak_id)],
       file.path(bed_out, "VD_favored_accessibility_gain.bed"),
       sep = "\t", col.names = FALSE)

me1_paths <- c(
  TAZ_gain = file.path(master, "H3K27me1",
    "TAZ_called.GRanges_H3K27me1_gain.var_1.logFC_0.5.TAZ.bed"),
  VAL_only_loss = file.path(master, "H3K27me1",
    "TAZ_called.GRanges_H3K27me1_lost.var_1.logFC_0.5.VAL.bed"),
  Shared_TAZ_VAL_loss = file.path(master, "H3K27me1",
    "VEH_called.GRanges_H3K27me1_lost.var_1.logFC_0.5.TAZ.bed")
)
stopifnot(all(file.exists(me1_paths)))
me1 <- lapply(me1_paths, read_bed)

distance_vector <- function(query, subject) {
  q <- GRanges(query$chr, IRanges(query$start + 1L, query$end))
  s <- GRanges(subject$chr, IRanges(subject$start + 1L, subject$end))
  h <- distanceToNearest(q, s, ignore.strand = TRUE)
  ans <- rep(NA_integer_, length(q))
  ans[queryHits(h)] <- mcols(h)$distance
  ans
}

# Export the two requested accessibility classes stratified by their nearby
# H3K27me1 state. State proximity is evaluated independently at <=10 kb.
state_bed_out <- file.path(out, "BEDs", "by_H3K27me1_state")
dir.create(state_bed_out, recursive = TRUE, showWarnings = FALSE)
class_bed_groups <- groups[c("Comparable shared", "VD-favored")]
state_manifest <- rbindlist(lapply(names(class_bed_groups), function(cl) {
  q <- class_bed_groups[[cl]]
  gain_dist <- distance_vector(q, me1[["TAZ_gain"]])
  maintain_dist <- distance_vector(q, me1[["VAL_only_loss"]])
  gain_keep <- !is.na(gain_dist) & gain_dist <= 10000L
  # Mutually exclusive assignment: gain has priority over maintenance.
  maintain_keep <- !gain_keep & !is.na(maintain_dist) & maintain_dist <= 10000L
  safe_class <- if (cl == "Comparable shared") "Comparable_TD_VD" else "VD_favored"
  assignments <- list(
    TAZ_gain = list(keep = gain_keep, label = "TAZ_gain"),
    TAZ_maintain = list(keep = maintain_keep, label = "TAZ_maintain_VAL_only_loss")
  )
  rbindlist(lapply(names(assignments), function(st) {
    a <- assignments[[st]]
    fn <- paste0(safe_class, "__", a$label, "__within10kb.bed")
    fwrite(q[a$keep, .(chr, start, end, peak_id)],
           file.path(state_bed_out, fn), sep = "\t", col.names = FALSE)
    data.table(accessibility_class = cl, H3K27me1_state = st,
               assignment_priority = "TAZ_gain_then_TAZ_maintain",
               peaks = sum(a$keep), bed_file = fn)
  }))
}))
fwrite(state_manifest, file.path(state_bed_out, "BED_manifest.tsv"), sep = "\t")

windows <- c(seq(0L, 1000L, 100L), seq(2000L, 10000L, 1000L))
curve <- rbindlist(lapply(names(groups), function(g) {
  q <- groups[[g]]
  rbindlist(lapply(names(me1), function(st) {
    dst <- distance_vector(q, me1[[st]])
    rbindlist(lapply(windows, function(w) data.table(
      group = g, state = st, window_bp = w,
      near = sum(dst <= w, na.rm = TRUE), total = nrow(q),
      percent = 100 * mean(dst <= w, na.rm = TRUE)
    )))
  }))
}))

group_levels <- c("DAC", "Comparable shared", "VD-favored")
state_levels <- c("TAZ_gain", "VAL_only_loss", "Shared_TAZ_VAL_loss")
state_labels <- c(
  TAZ_gain = "TAZ-specific H3K27me1 gain",
  VAL_only_loss = "TAZ-maintained / VAL-only H3K27me1 loss",
  Shared_TAZ_VAL_loss = "Shared TAZ/VAL H3K27me1 loss"
)
curve[, `:=`(
  group = factor(group, levels = group_levels),
  state = factor(state, levels = state_levels),
  H3K27me1_state = factor(unname(state_labels[as.character(state)]),
                          levels = unname(state_labels[state_levels]))
)]
fwrite(curve, file.path(out, "DAC_comparable_VDfavored_H3K27me1_cumulative_distance.tsv"), sep = "\t")

endpoint <- curve[window_bp == 10000L]
endpoint[, pct_label := sprintf("%.1f%%", percent)]
comparisons <- data.table(
  numerator = c("VD-favored", "VD-favored", "Comparable shared"),
  denominator = c("Comparable shared", "DAC", "DAC"),
  Contrast = c("VD-favored vs comparable", "VD-favored vs DAC", "Comparable vs DAC")
)
stats <- rbindlist(lapply(seq_len(nrow(comparisons)), function(i) {
  cp <- comparisons[i]
  rbindlist(lapply(state_levels, function(st) {
    a <- endpoint[group == cp$numerator & state == st]
    b <- endpoint[group == cp$denominator & state == st]
    ft <- fisher.test(matrix(c(a$near, a$total-a$near,
                               b$near, b$total-b$near), nrow = 2, byrow = TRUE))
    data.table(
      Contrast = cp$Contrast, state = st,
      numerator = cp$numerator, denominator = cp$denominator,
      numerator_near = a$near, numerator_total = a$total,
      denominator_near = b$near, denominator_total = b$total,
      OR = unname(ft$estimate), CI_lower = unname(ft$conf.int[1]),
      CI_upper = unname(ft$conf.int[2]), p = ft$p.value
    )
  }))
}))
contrast_levels <- rev(comparisons$Contrast)
stats[, `:=`(
  Contrast = factor(Contrast, levels = contrast_levels),
  H3K27me1_state = factor(unname(state_labels[state]),
                          levels = unname(state_labels[state_levels])),
  log2OR = log2(OR), lo = log2(CI_lower), hi = log2(CI_upper)
)]
stats[, ptxt := fifelse(p == 0, "<2e-308",
                        fifelse(p < 0.001, format(p, digits = 2, scientific = TRUE), sprintf("%.3f", p)))]
stats[, label := sprintf("log2 OR %.2f; p=%s", log2OR, ptxt)]
stats[, label_x := pmax(log2OR, hi) + 0.06]
fwrite(stats, file.path(out, "DAC_comparable_VDfavored_H3K27me1_10kb_Fisher_tests.tsv"), sep = "\t")

cols <- c(DAC = "#00c7f5", `Comparable shared` = "#921c52", `VD-favored` = "#ff84ff")
contrast_cols <- c(
  "VD-favored vs comparable" = "#ff84ff",
  "VD-favored vs DAC" = "#c84fc8",
  "Comparable vs DAC" = "#921c52"
)
boxed <- theme_bw(base_size = 10.5) + theme(
  panel.border = element_rect(color = "black", fill = NA, linewidth = 0.95),
  panel.grid = element_blank(),
  strip.background = element_rect(fill = "white", color = "black", linewidth = 0.95),
  strip.text = element_text(face = "bold", size = 9.2),
  axis.text = element_text(color = "black"),
  axis.ticks = element_line(color = "black", linewidth = 0.6),
  legend.position = "top"
)

distance_theme <- theme_classic(base_size = 10.5) + theme(
  panel.grid = element_blank(),
  panel.border = element_blank(),
  axis.line = element_line(color = "black", linewidth = 0.8),
  axis.ticks = element_line(color = "black", linewidth = 0.6),
  axis.text = element_text(color = "black"),
  strip.background = element_blank(),
  strip.text = element_text(face = "bold", size = 9.5),
  panel.spacing.x = grid::unit(0.55, "lines"),
  legend.position = "top"
)

p_or <- ggplot(stats, aes(log2OR, Contrast, color = as.character(Contrast))) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray35", linewidth = 0.75) +
  geom_errorbar(aes(xmin = lo, xmax = hi), orientation = "y", width = 0.16, linewidth = 0.9) +
  geom_point(size = 3.4) +
  geom_text(aes(x = label_x, label = label), color = "black", hjust = 0,
            size = 2.3, fontface = "bold") +
  facet_wrap(~H3K27me1_state, nrow = 1, scales = "free_x") +
  scale_color_manual(values = contrast_cols) +
  scale_x_continuous(expand = expansion(mult = c(0.10, 0.52))) +
  labs(x = "log2 odds ratio", y = NULL, color = NULL,
       title = "Relative enrichment within 10 kb",
       subtitle = "Positive values favor the first group named") +
  guides(color = "none") + boxed

legend_labels <- c(
  DAC = "DAC (n = 3,711)",
  `Comparable shared` = "Comparable TD/VD (n = 5,654)",
  `VD-favored` = "VD-favored (n = 10,020)"
)
make_distance <- function(yvar, ylabel, suffix, free_y = FALSE) {
  y_scale <- if (yvar == "near" && !free_y) {
    scale_y_continuous(
      limits = c(0, 3500),
      breaks = c(0, 1000, 2000, 3000, 3500),
      labels = scales::label_comma(),
      expand = expansion(mult = c(0, 0))
    )
  } else {
    scale_y_continuous(
      labels = if (yvar == "percent") function(x) paste0(x, "%") else scales::label_comma(),
      expand = expansion(mult = c(0, 0.05))
    )
  }
  p_dist <- ggplot(curve, aes(window_bp / 1000, .data[[yvar]], color = group, group = group)) +
    geom_line(linewidth = 1.05) + geom_point(size = 1.75) +
    facet_wrap(~H3K27me1_state, nrow = 1,
               scales = if (free_y) "free_y" else "fixed") +
    scale_color_manual(values = cols, labels = legend_labels) +
    scale_x_continuous(breaks = c(0, 1, 3, 5, 7, 10), limits = c(0, 11.25),
                       expand = expansion(mult = c(0, 0))) +
    y_scale +
    labs(x = expression("Accessibility distance to "*Delta*"H3K27me1 (kb)"),
         y = ylabel, color = NULL,
         title = if (yvar == "percent") "Cumulative percentage of accessibility-gain peaks" else "Cumulative raw number of accessibility-gain peaks",
         subtitle = if (free_y) "Y-axis optimized independently for each H3K27me1 state" else NULL) + distance_theme

  if (yvar == "near") {
    p_dist <- p_dist +
      geom_text(
        data = endpoint,
        aes(x = window_bp / 1000, y = near, label = pct_label, color = group),
        inherit.aes = FALSE, hjust = -0.16, size = 3.0,
        fontface = "bold", show.legend = FALSE
      )
  }

  combined <- p_or / p_dist + plot_layout(heights = c(1, 1.25), guides = "collect") +
    plot_annotation(
      title = "DAC, comparable TD/VD, and VD-favored accessibility near H3K27me1 dynamics",
      subtitle = "VD-favored combines VD-only peaks with shared peaks whose VD response is at least 20% stronger than TD",
      tag_levels = "A",
      theme = theme(plot.title = element_text(face = "bold", size = 15),
                    plot.subtitle = element_text(size = 10),
                    plot.tag = element_text(face = "bold", size = 13))
    ) & theme(legend.position = "top")

  stem <- file.path(out, paste0("DAC_COMPARABLE_VS_VD_FAVORED_H3K27me1_OR_", suffix))
  ggsave(paste0(stem, ".png"), combined, width = 16, height = 9.5, dpi = 600, bg = "white")
  ggsave(paste0(stem, ".pdf"), combined, width = 16, height = 9.5, device = pdf)
}

make_distance("near", "Number of accessibility increases", "RAW_DISTANCE")
make_distance("percent", "% Accessibility increases", "PERCENT_DISTANCE")
make_distance("near", "Number of accessibility increases", "RAW_DISTANCE_FREE_Y", free_y = TRUE)

cat("Saved to:", out, "\n")
