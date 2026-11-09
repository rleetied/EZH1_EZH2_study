# RNA-seq analysis summary
# Final use: Figure 1G-J
# Detailed implementation: ./analysis_code/ (numbered in workflow order).

library(edgeR)
library(limma)
library(data.table)
library(ggplot2)

PROJECT_ROOT <- "RKO_TAZ_v_VAL_project"

# 1. Import the fourth column of STAR ReadsPerGene.out.tab files.
count_files <- list.files(
  file.path(PROJECT_ROOT, "RNA-seq", "STAR_counts"),
  pattern = "ReadsPerGene.out.tab$", full.names = TRUE
)
sample_table <- fread(file.path(PROJECT_ROOT, "metadata", "RNAseq_samples.tsv"))

# Production code used edgeR::readDGE(columns = c(1, 4)).
dge <- readDGE(count_files, columns = c(1, 4),
               group = sample_table$treatment)

# 2. Filter low-expression genes and apply TMM normalization.
keep <- rowSums(cpm(dge) > 0.25) >= 1
dge <- calcNormFactors(dge[keep, , keep.lib.sizes = FALSE])

# 3. Fit treatment effects with limma-voom quality weights.
treatment <- relevel(factor(sample_table$treatment), ref = "VEH")
design <- model.matrix(~ treatment)
v <- voomWithQualityWeights(dge, design, plot = FALSE)
fit <- eBayes(lmFit(v, design))

# Treatment-versus-vehicle and VD-versus-TD contrasts were extracted with
# limma. Benjamini-Hochberg FDR values were used throughout.

# 4. Cross-cell-line analysis used protein-coding genes shared by RKO and
# SW620. Genes were ranked by the mean of their within-cell-line variances.
# The primary feature space contained the top 1,000 genes; 2,000, 5,000, and
# 10,000 genes were evaluated as sensitivity analyses.

# 5. Final displays:
# - row-z-scored expression heatmaps across all treatments (Figure 1G);
# - shared-scale PCA for RKO and SW620 (Figure 1H);
# - RMS distance from the line-specific reference centroid and relative
#   dual-versus-selective inhibitor displacement (Figure 1I);
# - ranked VD-versus-TD gene-set enrichment using Hallmark, Reactome, KEGG,
#   canonical-pathway, CGP, and Ben-Porath collections (Figure 1J).

