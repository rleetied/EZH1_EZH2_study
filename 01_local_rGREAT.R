#!/usr/bin/env Rscript

# Uses the active R library paths. If needed, set R_LIBS_USER before running.
suppressPackageStartupMessages({
  library(rGREAT)
  library(data.table)
  library(GenomicRanges)
})

root <- "RKO_TAZ_v_VAL_project/Final_figure_analysis/"
bed_dir <- file.path(root, "outputs", "DAC_TD_VD_COMPARABLE_VS_VD_FAVORED_H3K27me1_OR_DISTANCE",
                     "BEDs", "by_H3K27me1_state")
out <- file.path(root, "outputs", "LOCAL_rGREAT_FINAL_REGIONS")
dir.create(out, recursive = TRUE, showWarnings = FALSE)

bed_paths <- c(
  Comparable_TAZ_gain = file.path(bed_dir, "Comparable_TD_VD__TAZ_gain__within10kb.bed"),
  Comparable_TAZ_maintain = file.path(bed_dir, "Comparable_TD_VD__TAZ_maintain_VAL_only_loss__within10kb.bed"),
  VD_favored_TAZ_gain = file.path(bed_dir, "VD_favored__TAZ_gain__within10kb.bed"),
  VD_favored_TAZ_maintain = file.path(bed_dir, "VD_favored__TAZ_maintain_VAL_only_loss__within10kb.bed")
)
stopifnot(all(file.exists(bed_paths)))

collections <- c(
  GO_BP = "GO:BP",
  GO_CC = "GO:CC",
  GO_MF = "GO:MF",
  HALLMARK = "MSigDB:H",
  REACTOME = "MSigDB:C2:CP:REACTOME"
)

read_regions <- function(path) {
  x <- fread(path, header = FALSE)
  setnames(x, c("chr", "start", "end", "peak_id"))
  x[, `:=`(chr = as.character(chr), start = as.integer(start), end = as.integer(end))]
  gr <- GRanges(x$chr, IRanges(x$start + 1L, x$end), peak_id = x$peak_id)
  list(dt = x, gr = gr)
}

association_long <- list()
association_collapsed <- list()
enrichment_all <- list()
run_manifest <- list()

for (set_name in names(bed_paths)) {
  message("Processing: ", set_name)
  reg <- read_regions(bed_paths[[set_name]])
  assoc_written <- FALSE

  for (collection_name in names(collections)) {
    collection <- collections[[collection_name]]
    message("  Collection: ", collection)
    started <- Sys.time()
    ans <- tryCatch({
      res <- great(
        reg$gr,
        gene_sets = collection,
        tss_source = "GREAT:hg38",
        mode = "basalPlusExt",
        basal_upstream = 5000,
        basal_downstream = 1000,
        extension = 1e6,
        cores = 2,
        verbose = FALSE
      )

      if (!assoc_written) {
        a <- getRegionGeneAssociations(res)
        ad <- as.data.table(as.data.frame(a))
        # CharacterList/IntegerList become comma-separated strings in data.frame.
        ad[, region_set := set_name]
        setcolorder(ad, c("region_set", setdiff(names(ad), "region_set")))
        association_collapsed[[set_name]] <- ad

        genes <- mcols(a)$annotated_genes
        distances <- mcols(a)$dist_to_TSS
        long <- rbindlist(lapply(seq_along(a), function(i) {
          if (!length(genes[[i]])) {
            return(data.table(
              region_set = set_name,
              peak_id = mcols(a)$peak_id[i],
              chr = as.character(seqnames(a)[i]),
              start = start(a)[i] - 1L,
              end = end(a)[i],
              gene = NA_character_,
              distance_to_TSS = NA_integer_
            ))
          }
          data.table(
            region_set = set_name,
            peak_id = mcols(a)$peak_id[i],
            chr = as.character(seqnames(a)[i]),
            start = start(a)[i] - 1L,
            end = end(a)[i],
            gene = as.character(genes[[i]]),
            distance_to_TSS = as.integer(distances[[i]])
          )
        }))
        association_long[[set_name]] <- long
        assoc_written <- TRUE
      }

      e <- as.data.table(getEnrichmentTable(res))
      e[, `:=`(region_set = set_name, collection = collection_name)]
      setcolorder(e, c("region_set", "collection", setdiff(names(e), c("region_set", "collection"))))
      enrichment_all[[paste(set_name, collection_name, sep = "__")]] <- e
      list(ok = TRUE, error = NA_character_, terms = nrow(e))
    }, error = function(err) {
      list(ok = FALSE, error = conditionMessage(err), terms = 0L)
    })

    run_manifest[[length(run_manifest) + 1L]] <- data.table(
      region_set = set_name,
      input_bed = bed_paths[[set_name]],
      input_regions = nrow(reg$dt),
      collection = collection_name,
      collection_code = collection,
      success = ans$ok,
      enrichment_terms = ans$terms,
      error = ans$error,
      elapsed_seconds = as.numeric(difftime(Sys.time(), started, units = "secs"))
    )
  }
}

manifest <- rbindlist(run_manifest, fill = TRUE)
fwrite(manifest, file.path(out, "rGREAT_run_manifest.tsv"), sep = "\t")

if (length(association_long)) {
  al <- rbindlist(association_long, use.names = TRUE, fill = TRUE)
  fwrite(al, file.path(out, "rGREAT_peak_to_gene_associations_LONG.tsv.gz"), sep = "\t")
  # A nearest-gene representative is useful for one-gene-per-peak downstream tables.
  nearest <- al[!is.na(gene)][order(abs(distance_to_TSS)), .SD[1], by = .(region_set, peak_id)]
  fwrite(nearest, file.path(out, "rGREAT_peak_to_nearest_associated_gene.tsv.gz"), sep = "\t")
}
if (length(association_collapsed)) {
  fwrite(rbindlist(association_collapsed, use.names = TRUE, fill = TRUE),
         file.path(out, "rGREAT_peak_to_gene_associations_COLLAPSED.tsv.gz"), sep = "\t")
}
if (length(enrichment_all)) {
  ea <- rbindlist(enrichment_all, use.names = TRUE, fill = TRUE)
  fwrite(ea, file.path(out, "rGREAT_enrichment_ALL.tsv.gz"), sep = "\t")
  sig <- ea[p_adjust < 0.05 | p_adjust_hyper < 0.05]
  setorder(sig, region_set, collection, p_adjust, p_adjust_hyper)
  fwrite(sig, file.path(out, "rGREAT_enrichment_FDR05.tsv"), sep = "\t")
  for (nm in unique(ea$collection)) {
    fwrite(ea[collection == nm], file.path(out, paste0("rGREAT_enrichment_", nm, ".tsv.gz")), sep = "\t")
  }
}

settings <- data.table(
  setting = c("genome", "TSS source", "association rule", "basal upstream bp",
              "basal downstream bp", "maximum extension bp", "background", "coordinate conversion"),
  value = c("hg38", "GREAT:hg38", "basalPlusExt", "5000", "1000", "1000000",
            "whole genome excluding assembly gaps (rGREAT default)",
            "BED zero-based start converted to one-based GRanges start")
)
fwrite(settings, file.path(out, "rGREAT_analysis_settings.tsv"), sep = "\t")

cat("Saved local GREAT analysis to:", out, "\n")
print(manifest)
