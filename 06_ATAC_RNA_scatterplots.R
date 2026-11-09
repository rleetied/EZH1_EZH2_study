setwd("RKO_TAZ_v_VAL_project/Final_figure_analysis/outputs/ATAC_RNA_scatter_20pct//")

test <- fread("ATAC_RNA_representative_peak_plot_gene_classification.tsv")

TD <- subset(test, Treatment == "TD")
VD <- subset(test, Treatment == "VD")

setwd("~/Desktop/EZH1_EZH2_ms/Post_CR_review_analysis/ATAC/GOI")

Combo_CRE <-subset(VD, regulatory_class == "Distal_candidate_enhancer" & both_significant == "TRUE" & direction_relationship == "Concordant" & peak_response_direction == "Increased" & GREAT_peak_category == "Combo_driven")
Combo_TSS <-subset(VD, regulatory_class == "TSS_proximal" & both_significant == "TRUE" & direction_relationship == "Concordant" & peak_response_direction == "Increased" & GREAT_peak_category == "Combo_driven")

VD_VAL_CRE <- subset(VD, regulatory_class == "Distal_candidate_enhancer" & both_significant == "TRUE" & direction_relationship == "Concordant" & peak_response_direction == "Increased" & GREAT_peak_category %in% c("VD_dominant", "VAL_driven", "VAL_or_TD")) 
VD_VAL_TSS <- subset(VD, regulatory_class == "TSS_proximal" & both_significant == "TRUE" & direction_relationship == "Concordant" & peak_response_direction == "Increased" & GREAT_peak_category %in% c("VD_dominant", "VAL_driven", "VAL_or_TD"))

write.table(Combo_CRE, file = "Combo_CRE_full_annotation.txt", sep = "\t", col.names = T, row.names = F, quote = F)
write.table(Combo_TSS, file = "Combo_TSS_full_annotation.txt", sep = "\t", col.names = T, row.names = F, quote = F)

write.table(VD_VAL_CRE, file = "VD_VAL_CRE_full_annotation.txt", sep = "\t", col.names = T, row.names = F, quote = F)
write.table(VD_VAL_TSS, file = "VD_VAL_TSS_full_annotation.txt", sep = "\t", col.names = T, row.names = F, quote = F)

Combo_CRE_gene <- Combo_CRE[,c(2)]
Combo_TSS_gene <- Combo_TSS[,c(2)]

VD_VAL_CRE_gene <- VD_VAL_CRE[,c(2)]
VD_VAL_TSS_gene <- VD_VAL_TSS[,c(2)]

write.table(Combo_CRE_gene, file = "Combo_CRE_gene.txt", sep = "\t", col.names = F, row.names = F, quote = F)
write.table(Combo_TSS_gene, file = "Combo_TSS_gene.txt", sep = "\t", col.names = F, row.names = F, quote = F)

write.table(VD_VAL_CRE_gene, file = "VD_VAL_CRE_gene.txt", sep = "\t", col.names = F, row.names = F, quote = F)
write.table(VD_VAL_TSS_gene, file = "VD_VAL_TSS_gene.txt", sep = "\t", col.names = F, row.names = F, quote = F)

#Lets make some heatmaps!
setwd("RKO_TAZ_v_VAL_project/Final_figure_analysis/outputs/unbiased_ATAC_H3K27me1_RNA_20pct/")
setwd("~/Desktop/EZH1_EZH2_ms/Post_CR_review_analysis/ATAC/GOI")

ATAC <- fread("VD_upregulated_genes_H3K27me1_accessibility_category.csv", header = T)
colnames(ATAC) <- c("ext_gene","H3K27me1","category","regulatory_class")

ATAC$regulatory_class <- sub(
  "^TSS_CRE$",
  "TSS",
  ATAC$regulatory_class
)

table(ATAC$regulatory_class)

# BRING in RNA-seq
#read in expression data for RKO cells
setwd("RKO_TAZ_v_VAL_project/PR001414_ROTS_RNA/trimmed_fastq/count_matrix/")
final.ann <- fread("RKO_TAZ_VAL_RNA.0.25_norm_counts.stats.ann.txt")

#x <- final.ann[,c(79,43:60)]

x <- final.ann[,c(79,8,13,28,18,33)]
GOI <- merge(x = GREAT_categories_subset, y = x, by = "ext_gene")

GOI[, chromatin_barrier_priority := fifelse(chromatin_barrier == "Adaptive H3K27me1 barrier", 1L, 2L)]
GOI[, regulatory_priority := fifelse(regulatory_element == "TSS", 1L, 2L)]
GOI[, accessibility_priority := fifelse(accessibility == "VD-favored", 1L, 2L)]

setorder(
  GOI,
  ext_gene,
  chromatin_barrier_priority,
  regulatory_priority,
  accessibility_priority
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



# Boxplot: 
library(data.table)
library(ggplot2)

# ----------------------------
# Convert from wide to long
# ----------------------------

plot_data <- melt(
  as.data.table(GOI_unique),
  id.vars = c(
    "ext_gene",
    "category",
    "regulatory_class"
  ),
  measure.vars = c(
    "LogFC.DAC",
    "LogFC.TAZ",
    "LogFC.VAL",
    "LogFC.TD",
    "LogFC.VD"
  ),
  variable.name = "Treatment",
  value.name = "Log2FC"
)

# Clean treatment names
plot_data[, Treatment := sub("^LogFC\\.", "", Treatment)]

# Control plotting order
plot_data[, category := factor(
  category,
  levels = c("VD_VAL", "Combo_driven")
)]

plot_data[, regulatory_class := factor(
  regulatory_class,
  levels = c("TSS", "CRE")
)]

plot_data[, Treatment := factor(
  Treatment,
  levels = c("DAC", "TAZ", "VAL", "TD", "VD")
)]

# ----------------------------
# Colors
# ----------------------------

treatment_colors <- c(
  "DAC" = "#00c7f5",
  "TAZ" = "#ec357a",
  "VAL" = "#008f8f",
  "TD"  = "#931100",
  "VD"  = "#D682FF"
)

# ----------------------------
# Boxplot
# ----------------------------

p <- ggplot(
  plot_data,
  aes(
    x = category,
    y = Log2FC,
    color = Treatment,
    fill = Treatment
  )
) +
  geom_hline(
    yintercept = 0,
    color = "black",
    linewidth = 0.45
  ) +
  geom_hline(
    yintercept = c(-1, 1),
    color = "gray50",
    linetype = "dashed",
    linewidth = 0.4
  ) +
  geom_point(
    aes(group = interaction(category, Treatment)),
    position = position_jitterdodge(
      jitter.width = 0.15,
      jitter.height = 0,
      dodge.width = 0.8,
      seed = 123
    ),
    size = 0.7,
    alpha = 0.16,
    stroke = 0
  ) +
  geom_boxplot(
    aes(group = interaction(category, Treatment)),
    position = position_dodge(width = 0.8),
    width = 0.68,
    outlier.shape = NA,
    alpha = 0.20,
    linewidth = 0.7
  ) +
  facet_wrap(
    ~ regulatory_class,
    nrow = 1,
    drop = FALSE
  ) +
  scale_color_manual(values = treatment_colors) +
  scale_fill_manual(values = treatment_colors) +
  labs(
    x = NULL,
    y = "RNA expression log2 fold change",
    color = "Treatment",
    title = "Expression by accessibility category and regulatory class"
  ) +
  guides(fill = "none") +
  theme_classic(base_size = 12) +
  theme(
    panel.border = element_rect(
      color = "black",
      fill = NA,
      linewidth = 0.7
    ),
    strip.background = element_rect(
      fill = "gray95",
      color = "black"
    ),
    strip.text = element_text(face = "bold"),
    axis.text.x = element_text(
      angle = 30,
      hjust = 1
    ),
    legend.position = "top",
    legend.title = element_text(face = "bold")
  )

p





p <- ggplot(
  plot_data,
  aes(
    x = regulatory_class,
    y = Log2FC,
    color = Treatment,
    fill = Treatment
  )
) +
  geom_hline(
    yintercept = 0,
    color = "black",
    linewidth = 0.45
  ) +
  geom_hline(
    yintercept = c(-1, 1),
    color = "gray20",
    linetype = "dashed",
    linewidth = 0.4
  ) +
  geom_point(
    aes(group = interaction(regulatory_class, Treatment)),
    position = position_jitterdodge(
      jitter.width = 0.3,
      jitter.height = 0,
      dodge.width = 0.8,
      seed = 123
    ),
    size = 0.5,
    alpha = 0.2,
    stroke = 0
  ) +
  geom_boxplot(
    aes(group = interaction(regulatory_class, Treatment)),
    position = position_dodge(width = 0.8),
    width = 0.68,
    outlier.shape = NA,
    alpha = 0.20,
    linewidth = 0.7
  ) +
  facet_wrap(
    ~ category,
    nrow = 1,
    drop = FALSE
  ) +
  scale_color_manual(values = treatment_colors) +
  scale_fill_manual(values = treatment_colors) +
  labs(
    x = NULL,
    y = "RNA expression log2 fold change",
    color = "Treatment",
    title = "Expression by regulatory class and accessibility category"
  ) +
  guides(fill = "none") +
  theme_classic(base_size = 12) +
  theme(
    panel.border = element_rect(
      color = "black",
      fill = NA,
      linewidth = 0.7
    ),
    strip.background = element_rect(
      fill = "gray95",
      color = "black"
    ),
    strip.text = element_text(face = "bold"),
    legend.position = "top",
    legend.title = element_text(face = "bold")
  ) +
coord_cartesian(ylim = c(-4,12))

p






RATIO <- fread("RATIO_increase.txt")
FDR <- fread("FDR_increase.txt")

colnames(RATIO) <- c("Gene_set","VD","Comparable")
colnames(FDR) <- c("Gene_set","VD","Comparable")


melt.RATIO <- melt(RATIO, id.vars = c("Gene_set"),
                   measure.vars = c(2:3),
                   value.name = "RATIO")


melt.FDR <- melt(FDR, id.vars = c("Gene_set"),
                 measure.vars = (c(2:3)),
                 value.name = "FDR")


merged_dt <- merge(melt.RATIO, melt.FDR, by = c("Gene_set", "variable"))

p <- ggplot(merged_dt, aes(x = Gene_set, y = RATIO, fill = -log10(FDR))) +
  geom_col(position = position_dodge(width = 0.8)) +  # side-by-side bars
  facet_wrap(~variable, ncol = 1, scales = "free_y") + # optional: one row per variable
  coord_flip() +
  theme_bw() +
  labs(x = "Gene_set", y = "RATIO", fill = expression(-log[10](FDR))) +
  scale_fill_gradient2(
    low = "black", mid = "#66ffff", high = "red",
    midpoint = 7,
    limits = c(0, 20)
  ) +
  theme(axis.text.y = element_text(size = 7))


p <- ggplot(merged_dt, aes(x = Gene_set, y = -log10(FDR), fill = RATIO)) +
  geom_col(position = position_dodge(width = 0.8)) +  # side-by-side bars
  facet_wrap(~variable, ncol = 2, scales = "free_y") + # optional: one row per variable
  coord_flip() +
  theme_bw() +
  labs(x = "Motif", y = "FDR", fill = "RATIO") +
  scale_fill_gradient2(
    low = "black", mid = "#66ffff", high = "red",
    midpoint = 0.05,
    limits = c(0, 0.1)
  ) +
  theme(axis.text.y = element_text(size = 7))

p <- ggplot(merged_dt, aes(x = Gene_set, y = -log10(FDR), fill = RATIO)) +
  geom_col(aes(group = variable), position = position_dodge(width = 0.9)) +  # side-by-side
  coord_flip() +
  theme_bw() +
  labs(x = "Motif", y = expression(-log[10](FDR)), fill = "RATIO") +
  scale_fill_gradient2(
    low = "black", mid = "#66ffff", high = "red",
    midpoint = 0.05,
    limits = c(0, 0.1)
  ) +
  theme(axis.text.y = element_text(size = 7))


# install.packages("ggbreak")  # if you haven't already
library(ggbreak)

p <- ggplot(merged_dt, aes(x = Gene_set, y = -log10(FDR), fill = RATIO)) +
  geom_col(aes(group = variable), position = position_dodge(width = 0.9)) +  # side-by-side
  coord_flip() +
  theme_bw() +
  labs(x = "Motif", y = expression(-log[10](FDR)), fill = "RATIO") +
  scale_fill_gradient2(
    low = "black", mid = "#66ffff", high = "red",
    midpoint = 0.03,
    limits = c(0, 0.1)
  ) +
  theme(axis.text.y = element_text(size = 7))

merged_dt$variable <- factor(merged_dt$variable, levels = rev(unique(merged_dt$variable)))
merged_dt$Gene_set <- factor(merged_dt$Gene_set, levels = c("REACTOME_RNA_POLYMERASE_II_TRANSCRIPTION",
                                                            "REACTOME_METABOLISM_OF_RNA",
                                                            "REACTOME_POST_TRANSLATIONAL_PROTEIN_MODIFICATION",
                                                            "REACTOME_CELL_CYCLE",
                                                            "HALLMARK_MYC_TARGETS_V1",
                                                            "REACTOME_MEMBRANE_TRAFFICKING",
                                                            "HALLMARK_G2M_CHECKPOINT",
                                                            "HALLMARK_MTORC1_SIGNALING",
                                                            "REACTOME_CELLULAR_RESPONSES_TO_STIMULI",
                                                            "HALLMARK_FATTY_ACID_METABOLISM",
                                                            "HALLMARK_OXIDATIVE_PHOSPHORYLATION"))





library(data.table)
library(ggplot2)
library(scales)

# Replace enrichment_results with your data.table/data.frame name
plot_data <- as.data.table(copy(merged_dt))

# ----------------------------
# Prepare plotting variables
# ----------------------------

plot_data[, significance := -log10(FDR)]

# Clean pathway names
plot_data[, pathway_label := Gene_set]

plot_data[, pathway_label := sub(
  "^(HALLMARK_|KEGG_|REACTOME_)",
  "",
  pathway_label
)]

plot_data[, pathway_label := gsub(
  "_",
  " ",
  pathway_label
)]

# Control x-axis order
plot_data[, variable := factor(
  variable,
  levels = c("VD", "Comparable")
)]

# Sort pathways by their largest enrichment ratio
pathway_order <- plot_data[
  ,
  .(maximum_ratio = max(RATIO, na.rm = TRUE)),
  by = pathway_label
][
  order(maximum_ratio)
]$pathway_label

plot_data[, pathway_label := factor(
  pathway_label,
  levels = pathway_order
)]

# ----------------------------
# Bubble plot
# ----------------------------

p <- ggplot(
  plot_data,
  aes(
    x = variable,
    y = pathway_label
  )
) +
  geom_point(
    aes(
      size = RATIO,
      fill = significance
    ),
    shape = 21,
    color = "black",
    stroke = 0.55,
    alpha = 0.90
  ) +
  scale_size_continuous(
    name = "Gene ratio",
    range = c(3, 11),
    labels = function(x) paste0(round(x * 100), "%")
  ) + 
  scale_fill_gradientn(
    name = expression(-log[10](FDR)),
    colors = c(
      "#F2F2F2",
      "#FFB56B",
      "#E34A33",
      "#8E0152"
    )
  ) +
  labs(
    x = NULL,
    y = NULL,
    title = "Pathway enrichment",
    subtitle = "Bubble size represents gene ratio; color represents statistical significance"
  ) +
  theme_classic(base_size = 12) +
  theme(
    panel.border = element_rect(
      color = "black",
      fill = NA,
      linewidth = 0.75
    ),
    panel.grid.major.y = element_line(
      color = "gray90",
      linewidth = 0.3
    ),
    axis.text.x = element_text(
      face = "bold",
      size = 12
    ),
    axis.text.y = element_text(
      size = 9,
      color = "black"
    ),
    axis.ticks.y = element_blank(),
    legend.position = "right",
    legend.title = element_text(
      face = "bold"
    ),
    plot.title = element_text(
      face = "bold",
      size = 14
    ),
    plot.margin = margin(
      t = 10,
      r = 15,
      b = 10,
      l = 10
    )
  )

p
