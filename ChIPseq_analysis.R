# ChIP-seq analysis summary
# Final use: Figures 3 and 4E-F
# Detailed implementation: ./analysis_code/ (numbered in workflow order).

library(data.table)
library(GenomicRanges)
library(ggplot2)

PROJECT_ROOT <- "RKO_TAZ_v_VAL_project"

# H3K27me1 signal was averaged across three VEH, TAZ, and VAL replicates.
# Per-region variance and treatment-versus-VEH log2 fold changes were computed.
# Final called regions used variance <=1 and |log2FC| >=0.5.

# TAZ-associated H3K27me1 classes:
# 1. TAZ-specific gain;
# 2. TAZ-maintained / VAL-only loss;
# 3. shared TAZ/VAL loss.
# Gain was prioritized over maintenance for mutually exclusive downstream
# assignments. ChromHMM states were assigned by maximum interval overlap and
# collapsed for display (for example, EnhA1/EnhA2 -> EnhA).

# Final summaries included:
# - numbers and genomic coverage of gain/loss intervals;
# - ChromHMM enrichment relative to the tested region background;
# - average H3K27me1/H3K27me3/H3K27ac signal around region or ATAC centers;
# - accessibility counts within cumulative 0-10-kb edge-to-edge windows.

