#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(data.table))

root <- "RKO_TAZ_v_VAL_project/Final_figure_analysis/"
great_dir <- file.path(root, "outputs", "LOCAL_rGREAT_FINAL_REGIONS")
outdir <- file.path(root, "outputs", "FINAL_rGREAT_RNA_PRIORITY_LINKS")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

great_file <- file.path(great_dir,
  "rGREAT_peak_to_gene_associations_LONG.tsv.gz")
old_great_file <- file.path(root, "outputs", "ATAC_RNA_peak_gene_map_20pct",
  "all_GREAT_peak_gene_links_with_ATAC_stats.tsv")
rna_file <- paste0(
  "RKO_TAZ_v_VAL_project/",
  "PR001414_ROTS_RNA/trimmed_fastq/count_matrix/",
  "RKO_TAZ_VAL_RNA.0.25_norm_counts.stats.ann.txt"
)

# New local GREAT peak-to-gene links for the four final region sets.
links <- fread(great_file)
stopifnot(all(c("region_set", "peak_id", "chr", "start", "end", "gene",
                "distance_to_TSS") %in% names(links)))

links[, `:=`(
  accessibility_class = fifelse(grepl("^VD_favored", region_set),
    "VD-favored", "Comparable TD/VD"),
  H3K27me1_state = fifelse(grepl("TAZ_gain$", region_set),
    "TAZ gain", "TAZ maintain / VAL-only loss"),
  absolute_distance_bp = abs(distance_to_TSS)
)]

# Recover the exact pre-existing GREAT basal-domain class whenever the same
# peak-gene link was present in the earlier analysis. The fallback implements
# the same GREAT basalPlusExt definition: 5 kb upstream and 1 kb downstream.
old <- fread(old_great_file)
old_class <- unique(old[, .(
  peak_id, gene,
  prior_regulatory_class = regulatory_class
)])
links <- merge(links, old_class, by = c("peak_id", "gene"),
  all.x = TRUE, sort = FALSE)
links[, regulatory_class := fifelse(
  !is.na(prior_regulatory_class),
  fifelse(prior_regulatory_class == "TSS_proximal", "TSS", "CRE"),
  fifelse(distance_to_TSS >= -5000 & distance_to_TSS <= 1000, "TSS", "CRE")
)]
links[, regulatory_class_source := fifelse(
  !is.na(prior_regulatory_class), "prior_exact_GREAT_link",
  "local_GREAT_basal_domain_fallback"
)]
links[, prior_regulatory_class := NULL]

# ATAC values are peak-level, so take one unambiguous record per peak from the
# earlier consensus-peak table and attach all five treatment contrasts.
stat_cols <- c(
  "peak_id", "TAZ_vs_Veh_log2FC", "TAZ_vs_Veh_padj",
  "VAL_vs_Veh_log2FC", "VAL_vs_Veh_padj",
  "DAC_vs_Veh_log2FC", "DAC_vs_Veh_padj",
  "TAZ_DAC_vs_Veh_log2FC", "TAZ_DAC_vs_Veh_padj",
  "VAL_DAC_vs_Veh_log2FC", "VAL_DAC_vs_Veh_padj"
)
peak_stats <- unique(old[, ..stat_cols], by = "peak_id")
stopifnot(anyDuplicated(peak_stats$peak_id) == 0L)
links <- merge(links, peak_stats, by = "peak_id", all.x = TRUE, sort = FALSE)

# Preserve all GREAT links before any one-gene prioritization.
setcolorder(links, c(
  "gene", "peak_id", "chr", "start", "end", "region_set",
  "accessibility_class", "H3K27me1_state", "regulatory_class",
  "distance_to_TSS", "absolute_distance_bp", "regulatory_class_source",
  setdiff(names(links), c("gene", "peak_id", "chr", "start", "end",
    "region_set", "accessibility_class", "H3K27me1_state",
    "regulatory_class", "distance_to_TSS", "absolute_distance_bp",
    "regulatory_class_source"))
))
fwrite(links, file.path(outdir,
  "ALL_peak_gene_links_with_final_categories_and_ATAC_stats.tsv.gz"), sep = "\t")

# RNA symbols mapping to multiple rows are exported and excluded rather than
# choosing an arbitrary Ensembl/transcript record.
rna <- as.data.table(read.table(rna_file, header = TRUE, check.names = FALSE,
  stringsAsFactors = FALSE, quote = "\""))
rna_cols <- c(
  "ext_gene", "ENSG", "ens_gene", "gene_biotype",
  "LogFC.DAC", "FDR.DAC", "logCPM.DAC",
  "LogFC.TAZ", "FDR.TAZ", "logCPM.TAZ",
  "LogFC.VAL", "FDR.VAL", "logCPM.VAL",
  "LogFC.TD", "FDR.TD", "logCPM.TD",
  "LogFC.VD", "FDR.VD", "logCPM.VD"
)
stopifnot(all(rna_cols %in% names(rna)))
ambiguous <- rna[, .N, by = ext_gene][N > 1L, ext_gene]
fwrite(rna[ext_gene %in% ambiguous],
  file.path(outdir, "RNA_ambiguous_gene_symbols_excluded.tsv"), sep = "\t")
rna_unique <- rna[!ext_gene %in% ambiguous, ..rna_cols]
stopifnot(anyDuplicated(rna_unique$ext_gene) == 0L)

all_joined <- merge(links, rna_unique, by.x = "gene", by.y = "ext_gene",
  all.x = TRUE, sort = FALSE)
all_joined[, RNA_match_status := fifelse(
  gene %in% ambiguous, "ambiguous_symbol_excluded",
  fifelse(is.na(ENSG), "no_unique_RNA_match", "unique_RNA_match")
)]
fwrite(all_joined, file.path(outdir,
  "ALL_peak_gene_links_with_RNA.tsv.gz"), sep = "\t")

# Expression-blind representative selection, exactly in the requested order:
# TSS > CRE; TAZ gain > maintain; VD-favored > comparable; shortest |TSS
# distance; deterministic peak ID tie-breaker.
all_joined[, `:=`(
  regulatory_priority = fifelse(regulatory_class == "TSS", 1L, 2L),
  H3K27me1_priority = fifelse(H3K27me1_state == "TAZ gain", 1L, 2L),
  accessibility_priority = fifelse(accessibility_class == "VD-favored", 1L, 2L)
)]
setorder(all_joined, gene, regulatory_priority, H3K27me1_priority,
  accessibility_priority, absolute_distance_bp, peak_id)
representative <- all_joined[, .SD[1L], by = gene]
representative[, c("regulatory_priority", "H3K27me1_priority",
  "accessibility_priority") := NULL]
stopifnot(anyDuplicated(representative$gene) == 0L)

fwrite(representative, file.path(outdir,
  "ONE_representative_peak_per_gene_PRIORITY_RESOLVED.tsv"), sep = "\t")
fwrite(representative[RNA_match_status == "unique_RNA_match"],
  file.path(outdir,
    "ONE_representative_peak_per_gene_WITH_RNA_PRIORITY_RESOLVED.tsv"),
  sep = "\t")

# Plot-ready long table for the upcoming Combo/expression scatterplot. Both TD
# and VD values are retained for the identical representative gene/peak set.
rep_rna <- representative[RNA_match_status == "unique_RNA_match"]
plot_long <- rbindlist(list(
  rep_rna[, .(
    gene, ENSG, ens_gene, gene_biotype, peak_id, chr, start, end,
    region_set, accessibility_class, H3K27me1_state, regulatory_class,
    distance_to_TSS, absolute_distance_bp,
    Treatment = "TD",
    ATAC_log2FC = TAZ_DAC_vs_Veh_log2FC,
    ATAC_padj = TAZ_DAC_vs_Veh_padj,
    RNA_log2FC = LogFC.TD, RNA_FDR = FDR.TD, RNA_logCPM = logCPM.TD
  )],
  rep_rna[, .(
    gene, ENSG, ens_gene, gene_biotype, peak_id, chr, start, end,
    region_set, accessibility_class, H3K27me1_state, regulatory_class,
    distance_to_TSS, absolute_distance_bp,
    Treatment = "VD",
    ATAC_log2FC = VAL_DAC_vs_Veh_log2FC,
    ATAC_padj = VAL_DAC_vs_Veh_padj,
    RNA_log2FC = LogFC.VD, RNA_FDR = FDR.VD, RNA_logCPM = logCPM.VD
  )]
))
fwrite(plot_long, file.path(outdir,
  "Combo_ATAC_RNA_scatter_INPUT_LONG.tsv.gz"), sep = "\t")

summary <- rbindlist(list(
  data.table(metric = "all_peak_gene_links", value = nrow(links)),
  data.table(metric = "unique_peaks_in_links", value = uniqueN(links$peak_id)),
  data.table(metric = "unique_GREAT_genes", value = uniqueN(links$gene)),
  data.table(metric = "ambiguous_RNA_symbols_excluded", value = length(ambiguous)),
  data.table(metric = "representative_genes", value = nrow(representative)),
  data.table(metric = "representative_genes_with_unique_RNA_match",
    value = nrow(rep_rna)),
  data.table(metric = "plot_long_rows_TD_plus_VD", value = nrow(plot_long))
))
fwrite(summary, file.path(outdir, "integration_summary.tsv"), sep = "\t")

priority <- data.table(
  rank = 1:5,
  rule = c(
    "TSS over CRE",
    "TAZ gain over TAZ maintain / VAL-only loss",
    "VD-favored over Comparable TD/VD",
    "Smallest absolute GREAT distance to TSS",
    "Lexical peak_id tie-breaker"
  )
)
fwrite(priority, file.path(outdir, "representative_gene_priority_rules.tsv"),
  sep = "\t")

print(summary)
print(representative[, .N,
  by = .(regulatory_class, H3K27me1_state, accessibility_class)][order(
    regulatory_class, H3K27me1_state, accessibility_class)])
