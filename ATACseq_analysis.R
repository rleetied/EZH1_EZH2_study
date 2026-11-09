# ATAC-seq analysis summary
# Final use: Figures 3F-I and 4A-D
# Detailed implementation: ./analysis_code/ (numbered in workflow order).

library(data.table)
library(DESeq2)
library(GenomicRanges)

PROJECT_ROOT <- "RKO_TAZ_v_VAL_project"
counts <- fread(file.path(PROJECT_ROOT, "ATAC-seq",
                          "ATAC_consensus_featureCounts.txt"), skip = 1)
metadata <- fread(file.path(PROJECT_ROOT, "metadata", "ATAC_samples.tsv"))

# Consensus peaks were retained with >=10 fragments in >=3 samples.
count_matrix <- as.matrix(counts[, -(1:6)])
rownames(count_matrix) <- counts$Geneid
keep <- rowSums(count_matrix >= 10) >= 3

# Three experimental blocks contained VEH, DAC, TAZ, VAL, TD, and VD.
dds <- DESeqDataSetFromMatrix(
  countData = count_matrix[keep, ],
  colData = as.data.frame(metadata),
  design = ~ block + condition
)
dds <- DESeq(dds)

# Peaks were called increased at FDR <0.05 and log2FC >=1, and decreased at
# FDR <0.05 and log2FC <=-1. Variance-stabilized counts were used for QC/PCA.

# TD/VD increased peaks were classified in a common consensus-peak universe:
# - VD-only;
# - shared but VD log2FC >=1.20 x TD log2FC;
# - comparable TD/VD response.
# VD-only and shared peaks >=20% stronger in VD formed the VD-favored class
# (n=10,020); the comparable class contained 5,654 peaks.

# Edge-to-edge distances to H3K27me1 intervals were calculated with
# GenomicRanges::distanceToNearest(ignore.strand=TRUE) through 10 kb.
# H3K27me1 gain was prioritized when a peak was near both gain and maintenance.

