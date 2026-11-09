#load libraries
library(edgeR)
library(ggplot2)
library(biomaRt)
library(reshape2)
library(grid)
library(plyr)
library(dplyr)
library(limma)
library(stringr)
library(data.table)
library(gplots)
library(colorRamps)

#Load in the data matrix
setwd("RKO_TAZ_v_VAL_project/RNA-seq/SW620/")
SW620 <- fread("SW620.DAC_TAZ_VAL_treatments_D3.ann.stats.0.25.txt")

setwd("RKO_TAZ_v_VAL_project/RNA-seq/RKO/")
RKO <- fread("RKO.DAC_TAZ_VAL_treatments_D3.ann.stats.0.25.txt")

######playing around with plottying specific things
Sample <- RKO

Sample <- mutate(Sample, expr_stats.VALvTAZ = ifelse(Sample$LogFC.VAL >= 1 & Sample$FDR.VAL <= 0.01 & Sample$LogFC.TAZ >= 1 & Sample$FDR.TAZ <= 0.01, "common.upregulation",
                                                           ifelse(Sample$LogFC.VAL <= -1 & Sample$FDR.VAL <= 0.01 & Sample$LogFC.TAZ <= -1 & Sample$FDR.TAZ <= 0.01, "common.downregulation",
                                                                  ifelse(Sample$LogFC.VAL >= 1 & Sample$FDR.VAL <= 0.01, "VAL.upregulation",
                                                                         ifelse(Sample$LogFC.VAL <= -1 & Sample$FDR.VAL <= 0.01, "VAL.downregulation",
                                                                                ifelse(Sample$LogFC.TAZ >=1 & Sample$FDR.TAZ <= 0.01, "TAZ.upregulation",
                                                                                       ifelse(Sample$LogFC.TAZ <= -1 & Sample$FDR.TAZ <= 0.01, "TAZ.downregulation","n.s")))))))

Sample <- mutate(Sample, expr_stats.VDvTD = ifelse(Sample$LogFC.VD >= 1 & Sample$FDR.VD <= 0.01 & Sample$LogFC.TD >= 1 & Sample$FDR.TD <= 0.01, "common.upregulation",
                                                         ifelse(Sample$LogFC.VD <= -1 & Sample$FDR.VD <= 0.01 & Sample$LogFC.TD <= -1 & Sample$FDR.TD <= 0.01, "common.downregulation",
                                                                ifelse(Sample$LogFC.VD >= 1 & Sample$FDR.VD <= 0.01, "VD.upregulation",
                                                                       ifelse(Sample$LogFC.VD <= -1 & Sample$FDR.VD <= 0.01, "VD.downregulation",
                                                                              ifelse(Sample$LogFC.TD >=1 & Sample$FDR.TD <= 0.01, "TD.upregulation",
                                                                                     ifelse(Sample$LogFC.TD <= -1 & Sample$FDR.TD <= 0.01, "TD.downregulation","n.s")))))))

#Log2FC heatmaps for upregulation
x <- Sample[,c(51,53,54,1,7,12,17,22)]


GOI <- subset(x, expr_stats.VALvTAZ == "common.upregulation" | expr_stats.VALvTAZ == "TAZ.upregulation" | expr_stats.VALvTAZ == "VAL.upregulation" |
                expr_stats.VDvTD == "common.upregulation" | expr_stats.VDvTD == "TD.upregulation" | expr_stats.VDvTD == "VD.upregulation")

GOI[, expr_stats.VDvTD_priority := fifelse(expr_stats.VDvTD == "common.upregulation", 1L, 2L)]
GOI[, expr_stats.VALvTAZ_priority := fifelse(expr_stats.VALvTAZ == "common.upregulation", 1L, 2L)]


setorder(
  GOI,
  ext_gene,
  expr_stats.VDvTD_priority,
  expr_stats.VALvTAZ_priority
)

GOI_unique <- GOI[, .SD[1L], by = ext_gene]

# Remove temporary columns
GOI_unique[, c("chromatin_barrier_priority","regulatory_priority", "accessibility_priority") := NULL]
# ----------------------------
# Step 1: Prepare heatmap matrix and row annotation
# ----------------------------
# Set rownames for GOI_subset
rownames(GOI_unique) <- GOI_unique$ext_gene

# Numeric matrix for heatmap
#heatmap_matrix <- as.matrix(GOI_subset[, .(Veh1, Veh2, Veh3, DAC1, DAC2, DAC3,
#                                    TAZ1, TAZ2, TAZ3, VAL1, VAL2, VAL3, 
#                                    TD1, TD2, TD3, VD1, VD2, VD3)])


heatmap_matrix <- as.matrix(GOI_unique[, .(LogFC.DAC, LogFC.TAZ, LogFC.VAL, LogFC.TD, LogFC.VD)])
rownames(heatmap_matrix) <- GOI_unique$ext_gene

# ----------------------------
# Step 1: Row annotations
# ----------------------------
# Desired ordering
accessibility_order <- c(
  "VD-favored",
  "Comparable"
)

regulatory_order <- c(
  "TSS",
  "CRE"
)

chromatin_barrier_order <- c(
  "Adaptive H3K27me1 barrier",
  "H3K27me1 + induced bivalency",
  "Induced bivalency"
)
# Row annotations
row_annotation <- data.frame(
  accessibility = factor(
    GOI_unique$accessibility,
    levels = accessibility_order
  ),
  Regulatory_class = factor(
    GOI_unique$regulatory_element,
    levels = regulatory_order
  ),
  chromatin_barrier = factor(
    GOI_unique$chromatin_barrier,
    levels = chromatin_barrier_order),
  row.names = GOI_unique$ext_gene
)

# Sort first by category, then regulatory class, then VD expression
sort_index <- order(
  row_annotation$accessibility,
  row_annotation$chromatin_barrier,
  -heatmap_matrix[, "LogFC.VD"]
)

heatmap_matrix_sorted <- heatmap_matrix[sort_index, , drop = FALSE]
row_annotation_sorted <- row_annotation[sort_index, , drop = FALSE]

# Confirm identical row order
stopifnot(
  identical(
    rownames(heatmap_matrix_sorted),
    rownames(row_annotation_sorted)
  )
)
# ----------------------------
# Step 2: Annotation colors
# ----------------------------

annotation_colors <- list(
  chromatin_barrier = c(
    "Adaptive H3K27me1 barrier" = "#e9bbda",
    "H3K27me1 + induced bivalency" = "#cca84e",
    "Induced bivalency" = "#7bb6d7"
  ),
  Category = c(
    "VD-favored" = "#D682FF",
    "Comparable" = "#931100"
  ),
  Regulatory_class = c(
    "TSS"     = "black",
    "CRE"     = "gray70")
)

# Heatmap color scale
col <- colorRampPalette(
  c("blue", "cornsilk", "red")
)(100)
#c("#D682FF","#931100","#009091","#306791","#00c7f5","#bdbdbd")
# Color scale for heatmap
col <- colorRampPalette(c("blue","cornsilk","red"))(100)
#col <- colorRampPalette(c("black","red","gold"))(50)
#col <- colorRampPalette(c("black","gray50","red"))(50)


#color_limit <- quantile(
# abs(heatmap_matrix),
#  0.99,
#  na.rm = TRUE
#)

#breaks <- seq(
#  -color_limit,
#  color_limit,
#  length.out = length(col) + 1
#)

breaks <- seq(-8, 8, length.out = length(col)+1)

# ----------------------------
# Step 3: Sort rows by Category and LogFC.VD descending
# ----------------------------
#sort_index <- order(row_annotation$Category, -heatmap_matrix[, "LogFC.VD"])
#heatmap_matrix_sorted <- heatmap_matrix[sort_index, ]
#row_annotation_sorted <- row_annotation[sort_index, , drop = FALSE]

# Check rownames match
#stopifnot(all(rownames(heatmap_matrix_sorted) == row_annotation_sorted$Category))

# ----------------------------
# Step 4: Plot heatmap
# ----------------------------
pheatmap(
  heatmap_matrix_sorted,
  cluster_rows = FALSE,          # preserve manual sort
  cluster_cols = FALSE,
  annotation_row = row_annotation_sorted,
  annotation_colors = annotation_colors,
  show_rownames = TRUE,
  fontsize_row = 8,
  main = "Gene LogFC Heatmap by Category",
  color = col,
  scale = "none",
  breaks = breaks
)

