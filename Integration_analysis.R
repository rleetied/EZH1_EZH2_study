# Multi-omic integration analysis summary
# Final use: Figure 4C-J
# Detailed implementation: ./analysis_code/ (numbered in workflow order).

library(data.table)
library(GenomicRanges)
library(rGREAT)
library(ggplot2)

PROJECT_ROOT <- "RKO_TAZ_v_VAL_project"

# 1. Begin with significant TD/VD accessibility gains from the common ATAC
# consensus universe. Define VD-favored peaks as VD-only peaks plus shared
# peaks with VD log2FC >=1.20 x TD log2FC; all remaining shared peaks are
# comparable TD/VD responses.

# 2. Associate ATAC peaks with the closest H3K27me1 change by edge-to-edge
# distance through 10 kb. Prioritize TAZ H3K27me1 gain over TAZ-maintained /
# VAL-only loss when both are within 10 kb.

# 3. Link regions to genes with rGREAT using hg38 basal-plus-extension domains:
# 5 kb upstream, 1 kb downstream, and extension up to 1 Mb or the neighboring
# gene basal domain. Basal-domain links are labeled TSS; extended links are CRE.

# 4. Join GREAT-linked genes to RNA-seq by gene symbol. Ambiguous duplicated
# symbols are removed. Select one representative peak per gene without using
# RNA significance, with the hierarchy:
# TSS > CRE; H3K27me1 gain > maintenance > loss; VD-favored > comparable;
# shortest absolute TSS distance; peak identifier as deterministic tie-breaker.

# 5. Chromatin-barrier interpretation used in Figure 4:
# - VD-favored + H3K27me1 gain: adaptive H3K27me1 barrier;
# - comparable + gain or VD-favored + maintenance:
#   H3K27me1 plus induced bivalency;
# - comparable + maintenance: induced bivalency.

# 6. Plot VD-versus-VEH accessibility against VD-versus-VEH expression,
# faceted by accessibility class and TSS/CRE assignment. RNA thresholds were
# +1 and -0.5 log2FC; RNA FDR <0.05 was emphasized visually.

# 7. Generate the five-treatment gene-expression heatmap and pathway bubble
# plot from the final unique gene table. Bubble size represents gene ratio and
# color represents -log10(FDR), capped for display where indicated.

