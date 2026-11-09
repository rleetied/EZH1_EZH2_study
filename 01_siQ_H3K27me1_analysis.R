#H3K27me1_TAZVAL_processing_and_visualization
library(data.table)
library(ggplot2)
library(dplyr)

##need to remove overlapping annotations --- they are currently annotated based on proximity to the peak response coordinate. Overlap of 129 = fully within annotation; otherwise, chose the annotation that has more of an overlap
setwd("RKO_TAZ_v_VAL_project/ChIP-seq/siQ/siQ-ChIP/preanno/K27me1/tab/chromHMM_annotated/")

test <- fread("VEH1_vs_VEH2_K27me1.chromHMM.bed")
colnames(test) <- c("chr","start","stop","state", "overlap")


test$coordinate <- paste0(test$chr, ":", test$start)
overlaps <- subset(test, overlap != 129)

overlap_max <- overlaps %>%
  group_by(coordinate) %>%
  summarise(max_overlap = max(overlap))

overlap_max <- merge(x = overlap_max, y = overlaps, by = "coordinate")
overlap_max <- mutate(overlap_max, final = ifelse(max_overlap == overlap, TRUE, FALSE))
overlap_max <- subset(overlap_max, final == TRUE)
overlap_max <- overlap_max[,c(3:7,1)]

test_final <- subset(test, overlap == 129)
final <- rbind(test_final, overlap_max)

final <- final[,c(1:5)]
write.table(final, file = "VEH1_vs_VEH2_K27me1.chromHMM.final.bed", row.names = F, col.names = F, quote = F, sep ="\t")

###ChatGPT solution for for loop
# Set working directory
setwd("RKO_TAZ_v_VAL_project/ChIP-seq/siQ/siQ-ChIP/preanno/K27me1/tab/chromHMM_annotated/final_chromHMM/")

# List all .chromHMM.bed files in the directory
files <- list.files(pattern = "\\.chromHMM\\.bed$")

# Loop through each file
for (file in files) {
  
  # Read the current file
  test <- fread(file)
  
  # Assign column names
  colnames(test) <- c("chr", "start", "stop", "state", "overlap")
  
  # Create coordinate column
  test$coordinate <- paste0(test$chr, ":", test$start)
  
  # Filter rows where overlap is not equal to 129
  overlaps <- subset(test, overlap != 129)
  
  # Find maximum overlap for each coordinate
  overlap_max <- overlaps %>%
    group_by(coordinate) %>%
    summarise(max_overlap = max(overlap))
  
  # Merge with original overlaps to keep all information
  overlap_max <- merge(x = overlap_max, y = overlaps, by = "coordinate")
  
  # Mark the maximum overlap
  overlap_max <- mutate(overlap_max, final = ifelse(max_overlap == overlap, TRUE, FALSE))
  
  # Keep only the rows with the maximum overlap
  overlap_max <- subset(overlap_max, final == TRUE)
  
  # Reorder columns
  overlap_max <- overlap_max[, c(3:7, 1)]
  
  # Filter rows where overlap is 129
  test_final <- subset(test, overlap == 129)
  
  # Combine both filtered results
  final <- rbind(test_final, overlap_max)
  
  # Keep only relevant columns
  final <- final[, c(1:5)]
  
  # Write the modified data frame back to a file
  write.table(final, file = file, row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")
  
  # Print progress
  cat("Processed file:", file, "\n")
}

#VEH1
setwd("RKO_TAZ_v_VAL_project/ChIP-seq/siQ/consistent_peak_beds/H3K27me1/VEH/VEH1/")
test <- fread("VEH1_peaks.all_response_values.chromHMM.bed")

colnames(test) <- c("chr","start","stop", "VEH1.res","VEH2.res","VEH3.res",
                    "TAZ1.res","TAZ2.res","TAZ3.res","VAL1.res","VAL2.res","VAL3.res",
                    "max.chr","max.start","max.stop","state","overlap1","overlap2")

test <- test[,c(1:12,16)]

test <- mutate(test, VEH.avg_res = (VEH1.res + VEH2.res + VEH3.res)/3)
test <- mutate(test, TAZ.avg_res = (TAZ1.res + TAZ2.res + TAZ3.res)/3)
test <- mutate(test, VAL.avg_res = (VAL1.res + VAL2.res + VAL3.res)/3)

test$coordinate <- paste0(test$chr, ":", test$start)

test <- mutate(test, VEH.log2FC = log2(VEH.avg_res/VEH.avg_res))
test <- mutate(test, TAZ.log2FC = log2(TAZ.avg_res/VEH.avg_res))
test <- mutate(test, VAL.log2FC = log2(VAL.avg_res/VEH.avg_res))


write.table(test, file = "VEH1_all_response_values_K27me1.chromHMM.log2FC.bed", col.names = T, row.names = F, sep = "\t", quote = F)
write.table(test, file = "VEH1_all_response_values_K27me1.chromHMM.log2FC.noheader.bed", col.names = F, row.names = F, sep = "\t", quote = F)


#VEH2
setwd("RKO_TAZ_v_VAL_project/ChIP-seq/siQ/consistent_peak_beds/H3K27me1/VAL/VAL3")
test <- fread("VAL3_peaks.all_response_values.chromHMM.bed")

colnames(test) <- c("chr","start","stop", "VAL3.res","VEH1.res","VEH2.res",
                    "VEH3.res","TAZ1.res","TAZ2.res","TAZ3.res","VAL1.res","VAL2.res",
                    "max.chr","max.start","max.stop","state","overlap1","overlap2")

test <- test[,c(1:12,16)]

test <- mutate(test, VEH.avg_res = (VEH1.res + VEH2.res + VEH3.res)/3)
test <- mutate(test, TAZ.avg_res = (TAZ1.res + TAZ2.res + TAZ3.res)/3)
test <- mutate(test, VAL.avg_res = (VAL1.res + VAL2.res + VAL3.res)/3)

test$coordinate <- paste0(test$chr, ":", test$start)

test <- mutate(test, VEH.log2FC = log2(VEH.avg_res/VEH.avg_res))
test <- mutate(test, TAZ.log2FC = log2(TAZ.avg_res/VEH.avg_res))
test <- mutate(test, VAL.log2FC = log2(VAL.avg_res/VEH.avg_res))


write.table(test, file = "VAL3_all_response_values_K27me1.chromHMM.log2FC.bed", col.names = T, row.names = F, sep = "\t", quote = F)
write.table(test, file = "VAL3_all_response_values_K27me1.chromHMM.log2FC.noheader.bed", col.names = F, row.names = F, sep = "\t", quote = F)

melt_test.log2FC <- melt(test, id.vars = c("coordinate","state"),
                         measure.vars = c(18:20),
                         value.name = "log2FC")

#Differential_peaks
#Isolating out regions of interest for profiling with H3K27me1/H3K27me3/H3K27ac
setwd("RKO_TAZ_v_VAL_project/ChIP-seq/siQ/consistent_peak_beds/H3K27me1/TAZ/TAZ1/")
TAZ1 <- fread("TAZ1_all_response_values_K27me1.chromHMM.log2FC.bed", header = T)
TAZ1 <- mutate(TAZ1, TAZ.variance = ((TAZ1.res - TAZ.avg_res)^2 + (TAZ2.res - TAZ.avg_res)^2 + (TAZ3.res - TAZ.avg_res)^2)/2)
TAZ1 <- mutate(TAZ1, VEH.variance = ((VEH1.res - VEH.avg_res)^2 + (VEH2.res - VEH.avg_res)^2 + (VEH3.res - VEH.avg_res)^2)/2)
TAZ1 <- mutate(TAZ1, VAL.variance = ((VAL1.res - VAL.avg_res)^2 + (VAL2.res - VAL.avg_res)^2 + (VAL3.res - VAL.avg_res)^2)/2)


TAZ1.low_variance <- subset(TAZ1, TAZ.variance <= 0.25)
TAZ1.gain <- subset(TAZ1.low_variance, TAZ.log2FC >= 1.0)

melt_test.log2FC <- melt(TAZ1.low_variance, id.vars = c("coordinate","state"),
                         measure.vars = c(18:20),
                         value.name = "log2FC")

p <- ggplot(melt_test.log2FC, aes(x=variable, y = log2FC, fill = state)) +
  geom_boxplot(outlier.shape = NA) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90)) +
  geom_hline(yintercept = 0) +
  ylim(-5,5)


#####
setwd("RKO_TAZ_v_VAL_project/ChIP-seq/siQ/consistent_peak_beds/H3K27me1/VEH/VEH1/")
test <- fread("VAL3_all_response_values_K27me1.chromHMM.log2FC.bed", header = T)

#calculate variance
test <- mutate(test, TAZ.variance = ((TAZ1.res - TAZ.avg_res)^2 + (TAZ2.res - TAZ.avg_res)^2 + (TAZ3.res - TAZ.avg_res)^2)/2)
test <- mutate(test, VEH.variance = ((VEH1.res - VEH.avg_res)^2 + (VEH2.res - VEH.avg_res)^2 + (VEH3.res - VEH.avg_res)^2)/2)
test <- mutate(test, VAL.variance = ((VAL1.res - VAL.avg_res)^2 + (VAL2.res - VAL.avg_res)^2 + (VAL3.res - VAL.avg_res)^2)/2)

#threshold by variance
threshold <- 1.0
threshold <- 0.5
threshold <- 0.25
threshold <- 0.1

test.low_variance <- subset(test, VAL.variance <= threshold & VEH.variance <= threshold & TAZ.variance <= threshold)
variance <- test.low_variance[,c(1:3,13:16)]
write.table(variance, file = "VAL3_peaks.variance_threshold_0.1.bed", quote = F, sep = "\t", row.names = F, col.names = F)


#Log2FC 1.0
test.gain <- subset(test.low_variance, TAZ.log2FC >= 1.0)
test.gain <- test.gain[,c(1:3,13:16)]
write.table(test.gain, file = "VAL3_peaks.variance_threshold_0.1.TAZ_logFC_1.gain.bed", quote = F, sep = "\t", row.names = F, col.names = F)

test.lost <- subset(test.low_variance, TAZ.log2FC <= -1.0)
test.lost <- test.lost[,c(1:3,13:16)]
write.table(test.lost, file = "VAL3_peaks.variance_threshold_0.1.TAZ_logFC_1.lost.bed", quote = F, sep = "\t", row.names = F, col.names = F)

#Log2FC 0.5
test.gain <- subset(test.low_variance, TAZ.log2FC >= 0.5)
test.gain <- test.gain[,c(1:3,13:16)]
write.table(test.gain, file = "VAL3_peaks.variance_threshold_0.1.TAZ_logFC_0.5.gain.bed", quote = F, sep = "\t", row.names = F, col.names = F)

test.lost <- subset(test.low_variance, TAZ.log2FC <= -0.5)
test.lost <- test.lost[,c(1:3,13:16)]
write.table(test.lost, file = "VAL3_peaks.variance_threshold_0.1.TAZ_logFC_0.5.lost.bed", quote = F, sep = "\t", row.names = F, col.names = F)

#Log2FC 1.0
test.gain <- subset(test.low_variance, VAL.log2FC >= 1.0)
test.gain <- test.gain[,c(1:3,13:16)]
write.table(test.gain, file = "VAL3_peaks.variance_threshold_0.1.VAL_logFC_1.gain.bed", quote = F, sep = "\t", row.names = F, col.names = F)

test.lost <- subset(test.low_variance, VAL.log2FC <= -1.0)
test.lost <- test.lost[,c(1:3,13:16)]
write.table(test.lost, file = "VAL3_peaks.variance_threshold_0.1.VAL_logFC_1.lost.bed", quote = F, sep = "\t", row.names = F, col.names = F)


#overlap calculations
setwd("RKO_TAZ_v_VAL_project/ChIP-seq/siQ/K27me1_gain_peaks/")

TAZ1 <- fread("TAZ1_peaks.variance_threshold_0.1.TAZ_logFC_0.5.gain.bed")
TAZ2 <- fread("TAZ2_peaks.variance_threshold_0.1.TAZ_logFC_0.5.gain.bed")
TAZ3 <- fread("TAZ3_peaks.variance_threshold_0.1.TAZ_logFC_0.5.gain.bed")


TAZ1.gr <- GRanges(seqnames = TAZ1$V1, ranges = IRanges(start = TAZ1$V2, end = TAZ1$V3))
TAZ2.gr <- GRanges(seqnames = TAZ2$V1, ranges = IRanges(start = TAZ2$V2, end = TAZ2$V3))
TAZ3.gr <- GRanges(seqnames = TAZ3$V1, ranges = IRanges(start = TAZ3$V2, end = TAZ3$V3))

overlaps <- findOverlaps(TAZ1.gr, TAZ2.gr)
 
overlapping.TAZ1 <- data.frame(TAZ1.gr[queryHits(overlaps)])
overlapping.TAZ2 <- data.frame(TAZ2.gr[subjectHits(overlaps)])

write.table(overlapping.TAZ1, file = "TAZ.GRanges_overlap_peaks_TAZ_lost.bed", quote = F, col.names = F, row.names = F, sep = "\t")


setwd("RKO_TAZ_v_VAL_project/ChIP-seq/siQ/K27me1_lost_peaks/")

VEH1 <- fread("VEH1_peaks.variance_threshold_0.1.TAZ_logFC_0.5.lost.bed")
VEH2 <- fread("VEH2_peaks.variance_threshold_0.1.TAZ_logFC_0.5.lost.bed")
VEH3 <- fread("VEH3_peaks.variance_threshold_0.1.TAZ_logFC_0.5.lost.bed")


VEH1.gr <- GRanges(seqnames = VEH1$V1, ranges = IRanges(start = VEH1$V2, end = VEH1$V3))
VEH2.gr <- GRanges(seqnames = VEH2$V1, ranges = IRanges(start = VEH2$V2, end = VEH2$V3))
VEH3.gr <- GRanges(seqnames = VEH3$V1, ranges = IRanges(start = VEH3$V2, end = VEH3$V3))

overlaps <- findOverlaps(VEH1.gr, VEH2.gr)

overlapping.VEH1 <- data.frame(VEH1.gr[queryHits(overlaps)])
overlapping.VEH2 <- data.frame(VEH2.gr[subjectHits(overlaps)])

write.table(overlapping.VEH1, file = "VEH.GRanges_overlap_peaks_TAZ_lost.bed", quote = F, col.names = F, row.names = F, sep = "\t")



setwd("RKO_TAZ_v_VAL_project/ChIP-seq/siQ/K27me1_lost_peaks/")

VAL1 <- fread("VAL1_peaks.variance_threshold_0.1.TAZ_logFC_0.5.lost.bed")
VAL2 <- fread("VAL2_peaks.variance_threshold_0.1.TAZ_logFC_0.5.lost.bed")
VAL3 <- fread("VAL3_peaks.variance_threshold_0.1.TAZ_logFC_0.5.lost.bed")


VAL1.gr <- GRanges(seqnames = VAL1$V1, ranges = IRanges(start = VAL1$V2, end = VAL1$V3))
VAL2.gr <- GRanges(seqnames = VAL2$V1, ranges = IRanges(start = VAL2$V2, end = VAL2$V3))
VAL3.gr <- GRanges(seqnames = VAL3$V1, ranges = IRanges(start = VAL3$V2, end = VAL3$V3))

overlaps <- findOverlaps(VAL2.gr, VAL3.gr)
overlaps <- findOverlaps(VAL1.gr, VAL2.gr)

overlapping.VAL2 <- data.frame(VAL2.gr[queryHits(overlaps)])
overlapping.VAL3 <- data.frame(VAL3.gr[subjectHits(overlaps)])

overlapping.VAL1 <- data.frame(VAL1.gr[queryHits(overlaps)])
overlapping.VAL2 <- data.frame(VAL2.gr[subjectHits(overlaps)])


write.table(overlapping.VAL2, file = "VAL.GRanges_overlap_peaks_TAZ_lost.bed", quote = F, col.names = F, row.names = F, sep = "\t")


setwd("RKO_TAZ_v_VAL_project/ChIP-seq/siQ/K27ac_gain_peaks/")

DAC <- fread("DAC.GRanges_overlap_peaks_K27ac_gain.bed")
TAZ <- fread("TAZ.GRanges_overlap_peaks_K27ac_gain.bed")
VAL <- fread("VAL.GRanges_overlap_peaks_K27ac_gain.bed")
DTAZ <- fread("DTAZ.GRanges_overlap_peaks_K27ac_gain.bed")
DVAL <- fread("DVAL.GRanges_overlap_peaks_K27ac_gain.bed")


DAC.gr <- GRanges(seqnames = DAC$V1, ranges = IRanges(start = DAC$V2, end = DAC$V3))
TAZ.gr <- GRanges(seqnames = TAZ$V1, ranges = IRanges(start = TAZ$V2, end = TAZ$V3))
VAL.gr <- GRanges(seqnames = VAL$V1, ranges = IRanges(start = VAL$V2, end = VAL$V3))
DTAZ.gr <- GRanges(seqnames = DTAZ$V1, ranges = IRanges(start = DTAZ$V2, end = DTAZ$V3))
DVAL.gr <- GRanges(seqnames = DVAL$V1, ranges = IRanges(start = DVAL$V2, end = DVAL$V3))

TAZ_VAL.overlaps <- findOverlaps(TAZ.gr, VAL.gr)

overlapping.TAZ_VAL <- data.frame(TAZ.gr[queryHits(TAZ_VAL.overlaps)])
overlapping.VAL3 <- data.frame(VAL3.gr[subjectHits(overlaps)])

overlapping.VAL1 <- data.frame(VAL1.gr[queryHits(overlaps)])
overlapping.VAL2 <- data.frame(VAL2.gr[subjectHits(overlaps)])


setwd("RKO_TAZ_v_VAL_project/ChIP-seq/siQ/K27me1_gain_peaks/TAZ_gain/")

TAZ <- fread("TAZ.GRanges_overlap_peaks.bed")
VAL <- fread("VAL.GRanges_overlap_peaks.bed")

TAZ.gr <- GRanges(seqnames = TAZ$V1, ranges = IRanges(start = TAZ$V2, end = TAZ$V3))
VAL.gr <- GRanges(seqnames = VAL$V1, ranges = IRanges(start = VAL$V2, end = VAL$V3))

TAZ_VAL.overlaps <- findOverlaps(TAZ.gr, VAL.gr)

overlapping.TAZ_VAL <- data.frame(TAZ.gr[queryHits(TAZ_VAL.overlaps)])


overlapping.VAL1 <- data.frame(VAL1.gr[queryHits(overlaps)])
overlapping.VAL2 <- data.frame(VAL2.gr[subjectHits(overlaps)])

test <- fread("TAZ_gain_K27me1.all_peaks_combined.bed")

test$coordinate <- paste0(test$V1, ":", test$V2)
test.rm_dups <- test[!duplicated(test$coordinate),]

write.table(test.rm_dups, file = "TAZ_gain_K27me1.all_peaks_combined.bed", quote = F, col.names = F, row.names = F, sep = "\t")

setwd("RKO_TAZ_v_VAL_project/ChIP-seq/siQ/K27me1_lost_peaks/VAL_lost/")
test <- fread("VAL_loss_K27me1.cat_peaks.bed")

test$coordinate <- paste0(test$V1, ":", test$V2)
test.rm_dups <- test[!duplicated(test$coordinate),]

write.table(test.rm_dups, file = "VAL_loss_K27me1.all_peaks_combined.bed", quote = F, col.names = F, row.names = F, sep = "\t")

#overlap calculations
setwd("RKO_TAZ_v_VAL_project/ChIP-seq/siQ/consistent_peak_beds/H3K27me1/20250708_H3K27me1_variance_checks/0.5/H3K27me1_lost")

TAZ1 <- fread("TAZ1_peaks.variance_threshold_0.5.VAL_logFC_0.5.lost.bed")
TAZ2 <- fread("TAZ2_peaks.variance_threshold_0.5.VAL_logFC_0.5.lost.bed")
TAZ3 <- fread("TAZ3_peaks.variance_threshold_0.5.VAL_logFC_0.5.lost.bed")


TAZ1.gr <- GRanges(seqnames = TAZ1$V1, ranges = IRanges(start = TAZ1$V2, end = TAZ1$V3))
TAZ2.gr <- GRanges(seqnames = TAZ2$V1, ranges = IRanges(start = TAZ2$V2, end = TAZ2$V3))
TAZ3.gr <- GRanges(seqnames = TAZ3$V1, ranges = IRanges(start = TAZ3$V2, end = TAZ3$V3))

overlaps <- findOverlaps(TAZ1.gr, TAZ2.gr)

overlapping.TAZ1 <- data.frame(TAZ1.gr[queryHits(overlaps)])
overlapping.TAZ2 <- data.frame(TAZ2.gr[subjectHits(overlaps)])

write.table(overlapping.TAZ1, file = "TAZ.GRanges_overlap_peaks_TAZ_lost.bed", quote = F, col.names = F, row.names = F, sep = "\t")


setwd("RKO_TAZ_v_VAL_project/ChIP-seq/siQ/consistent_peak_beds/H3K27me1/20250708_H3K27me1_variance_checks/0.5/H3K27me1_lost")

VEH1 <- fread("VEH1_peaks.variance_threshold_0.5.VAL_logFC_0.5.lost.bed")
VEH2 <- fread("VEH2_peaks.variance_threshold_0.5.VAL_logFC_0.5.lost.bed")
VEH3 <- fread("VEH3_peaks.variance_threshold_0.5.VAL_logFC_0.5.lost.bed")


VEH1.gr <- GRanges(seqnames = VEH1$V1, ranges = IRanges(start = VEH1$V2, end = VEH1$V3))
VEH2.gr <- GRanges(seqnames = VEH2$V1, ranges = IRanges(start = VEH2$V2, end = VEH2$V3))
VEH3.gr <- GRanges(seqnames = VEH3$V1, ranges = IRanges(start = VEH3$V2, end = VEH3$V3))

overlaps <- findOverlaps(VEH1.gr, VEH2.gr)

overlapping.VEH1 <- data.frame(VEH1.gr[queryHits(overlaps)])
overlapping.VEH2 <- data.frame(VEH2.gr[subjectHits(overlaps)])

write.table(overlapping.VEH1, file = "VEH.GRanges_overlap_peaks_TAZ_lost.bed", quote = F, col.names = F, row.names = F, sep = "\t")



setwd("RKO_TAZ_v_VAL_project/ChIP-seq/siQ/K27me1_lost_peaks/")

VAL1 <- fread("VAL1_peaks.variance_threshold_0.1.TAZ_logFC_0.5.lost.bed")
VAL2 <- fread("VAL2_peaks.variance_threshold_0.1.TAZ_logFC_0.5.lost.bed")
VAL3 <- fread("VAL3_peaks.variance_threshold_0.1.TAZ_logFC_0.5.lost.bed")


VAL1.gr <- GRanges(seqnames = VAL1$V1, ranges = IRanges(start = VAL1$V2, end = VAL1$V3))
VAL2.gr <- GRanges(seqnames = VAL2$V1, ranges = IRanges(start = VAL2$V2, end = VAL2$V3))
VAL3.gr <- GRanges(seqnames = VAL3$V1, ranges = IRanges(start = VAL3$V2, end = VAL3$V3))

overlaps <- findOverlaps(VAL2.gr, VAL3.gr)
overlaps <- findOverlaps(VAL1.gr, VAL2.gr)

overlapping.VAL2 <- data.frame(VAL2.gr[queryHits(overlaps)])
overlapping.VAL3 <- data.frame(VAL3.gr[subjectHits(overlaps)])

overlapping.VAL1 <- data.frame(VAL1.gr[queryHits(overlaps)])
overlapping.VAL2 <- data.frame(VAL2.gr[subjectHits(overlaps)])


write.table(overlapping.VAL2, file = "VAL.GRanges_overlap_peaks_TAZ_lost.bed", quote = F, col.names = F, row.names = F, sep = "\t")


setwd("RKO_TAZ_v_VAL_project/ChIP-seq/siQ/K27ac_gain_peaks/")

DAC <- fread("DAC.GRanges_overlap_peaks_K27ac_gain.bed")
TAZ <- fread("TAZ.GRanges_overlap_peaks_K27ac_gain.bed")
VAL <- fread("VAL.GRanges_overlap_peaks_K27ac_gain.bed")
DTAZ <- fread("DTAZ.GRanges_overlap_peaks_K27ac_gain.bed")
DVAL <- fread("DVAL.GRanges_overlap_peaks_K27ac_gain.bed")


DAC.gr <- GRanges(seqnames = DAC$V1, ranges = IRanges(start = DAC$V2, end = DAC$V3))
TAZ.gr <- GRanges(seqnames = TAZ$V1, ranges = IRanges(start = TAZ$V2, end = TAZ$V3))
VAL.gr <- GRanges(seqnames = VAL$V1, ranges = IRanges(start = VAL$V2, end = VAL$V3))
DTAZ.gr <- GRanges(seqnames = DTAZ$V1, ranges = IRanges(start = DTAZ$V2, end = DTAZ$V3))
DVAL.gr <- GRanges(seqnames = DVAL$V1, ranges = IRanges(start = DVAL$V2, end = DVAL$V3))

TAZ_VAL.overlaps <- findOverlaps(TAZ.gr, VAL.gr)

overlapping.TAZ_VAL <- data.frame(TAZ.gr[queryHits(TAZ_VAL.overlaps)])
overlapping.VAL3 <- data.frame(VAL3.gr[subjectHits(overlaps)])

overlapping.VAL1 <- data.frame(VAL1.gr[queryHits(overlaps)])
overlapping.VAL2 <- data.frame(VAL2.gr[subjectHits(overlaps)])


setwd("RKO_TAZ_v_VAL_project/ChIP-seq/siQ/K27me1_gain_peaks/TAZ_gain/")

TAZ <- fread("TAZ.GRanges_overlap_peaks.bed")
VAL <- fread("VAL.GRanges_overlap_peaks.bed")

TAZ.gr <- GRanges(seqnames = TAZ$V1, ranges = IRanges(start = TAZ$V2, end = TAZ$V3))
VAL.gr <- GRanges(seqnames = VAL$V1, ranges = IRanges(start = VAL$V2, end = VAL$V3))

TAZ_VAL.overlaps <- findOverlaps(TAZ.gr, VAL.gr)

overlapping.TAZ_VAL <- data.frame(TAZ.gr[queryHits(TAZ_VAL.overlaps)])


overlapping.VAL1 <- data.frame(VAL1.gr[queryHits(overlaps)])
overlapping.VAL2 <- data.frame(VAL2.gr[subjectHits(overlaps)])

test <- fread("TAZ_gain_K27me1.all_peaks_combined.bed")

test$coordinate <- paste0(test$V1, ":", test$V2)
test.rm_dups <- test[!duplicated(test$coordinate),]

write.table(test.rm_dups, file = "TAZ_gain_K27me1.all_peaks_combined.bed", quote = F, col.names = F, row.names = F, sep = "\t")

setwd("RKO_TAZ_v_VAL_project/ChIP-seq/siQ/K27me1_lost_peaks/VAL_lost/")
test <- fread("VAL_loss_K27me1.cat_peaks.bed")

test$coordinate <- paste0(test$V1, ":", test$V2)
test.rm_dups <- test[!duplicated(test$coordinate),]

write.table(test.rm_dups, file = "VAL_loss_K27me1.all_peaks_combined.bed", quote = F, col.names = F, row.names = F, sep = "\t")


######this worked to narrow it down (look up), but then needed to do merge with bedtools to remove overlapping intervals
setwd("RKO_TAZ_v_VAL_project/ChIP-seq/siQ/20241208_final_called_peaks/")

DAC.loss <- fread("DAC.GRanges_overlap_peaks_K27ac_loss.bed")
TAZ.loss <- fread("TAZ.GRanges_overlap_peaks_K27ac_loss.bed")
VAL.loss <- fread("VAL.GRanges_overlap_peaks_K27ac_loss.bed")
DTAZ.loss <- fread("DTAZ.GRanges_overlap_peaks_K27ac_loss.bed")
DVAL.loss <- fread("DVAL.GRanges_overlap_peaks_K27ac_loss.bed")

DAC.loss$coordinate <- paste0(DAC.loss$V1, ":", DAC.loss$V2)
TAZ.loss$coordinate <- paste0(TAZ.loss$V1, ":", TAZ.loss$V2)
VAL.loss$coordinate <- paste0(VAL.loss$V1, ":", VAL.loss$V2)
DTAZ.loss$coordinate <- paste0(DTAZ.loss$V1, ":", DTAZ.loss$V2)
DVAL.loss$coordinate <- paste0(DVAL.loss$V1, ":", DVAL.loss$V2)


DAC.loss.rm_dups <- DAC.loss[!duplicated(DAC.loss$coordinate),]
TAZ.loss.rm_dups <- TAZ.loss[!duplicated(TAZ.loss$coordinate),]
VAL.loss.rm_dups <- VAL.loss[!duplicated(VAL.loss$coordinate),]
DTAZ.loss.rm_dups <- DTAZ.loss[!duplicated(DTAZ.loss$coordinate),]
DVAL.loss.rm_dups <- DVAL.loss[!duplicated(DVAL.loss$coordinate),]

write.table(DAC.loss.rm_dups, file = "DAC.GRanges_overlap_peaks_K27ac_loss.rm_dups.bed", quote = F, col.names = F, row.names = F, sep = "\t")
write.table(TAZ.loss.rm_dups, file = "TAZ.GRanges_overlap_peaks_K27ac_loss.rm_dups.bed", quote = F, col.names = F, row.names = F, sep = "\t")
write.table(VAL.loss.rm_dups, file = "VAL.GRanges_overlap_peaks_K27ac_loss.rm_dups.bed", quote = F, col.names = F, row.names = F, sep = "\t")
write.table(DTAZ.loss.rm_dups, file = "DTAZ.GRanges_overlap_peaks_K27ac_loss.rm_dups.bed", quote = F, col.names = F, row.names = F, sep = "\t")
write.table(DVAL.loss.rm_dups, file = "DVAL.GRanges_overlap_peaks_K27ac_loss.rm_dups.bed", quote = F, col.names = F, row.names = F, sep = "\t")


DAC.loss.rm_dups_count <- sum(DAC.loss.rm_dups$V3 - DAC.loss.rm_dups$V2)
TAZ.loss.rm_dups_count <- sum(TAZ.loss.rm_dups$V3 - TAZ.loss.rm_dups$V2)
VAL.loss.rm_dups_count <- sum(VAL.loss.rm_dups$V3 - VAL.loss.rm_dups$V2)
DTAZ.loss.rm_dups_count <- sum(DTAZ.loss.rm_dups$V3 - DTAZ.loss.rm_dups$V2)
DVAL.loss.rm_dups_count <- sum(DVAL.loss.rm_dups$V3 - DVAL.loss.rm_dups$V2)


#For Loop
setwd("RKO_TAZ_v_VAL_project/ChIP-seq/siQ/RKO_chromHMM_breakdown/bedtools_intersections/sample_intersections/TAZ")

# List all .bed files in the directory
file_list <- list.files(pattern = "\\.bed$")

# Loop through each file
for (file in file_list) {
  # Read the file
  TAZ <- fread(file)
  
  # Create the coordinate column
  TAZ$coordinate <- paste0(TAZ$V1, ":", TAZ$V2)
  
  # Select and remove duplicate rows
  TAZ.rm_dup <- TAZ[, c(1:4,16)]
  TAZ.rm_dup <- unique(TAZ.rm_dup)
  
  # Generate the output file name
  output_file <- sub("\\.bed$", ".rm_duplicates.bed", file)
  
  # Write the modified data to a new file
  write.table(
    TAZ.rm_dup,
    file = output_file,
    quote = FALSE,
    sep = "\t",
    col.names = FALSE,
    row.names = FALSE
  )
}

#ChromHMM_coverage_with_gained_peaks_and_such
setwd("RKO_TAZ_v_VAL_project/ChIP-seq/siQ/RKO_chromHMM_breakdown/chromHMM_instersection_for_coverage/")

TAZ <- fread("TAZ_loss_K27me1.intersect.chromHMM_ann.bed")

TssA <- subset(TAZ, TAZ$V7 == "TssA")
TssFlnk <- subset(TAZ, TAZ$V7 == "TssFlnk")
TssFlnkD <- subset(TAZ, TAZ$V7 == "TssFlnkD")
TssFlnkU <- subset(TAZ, TAZ$V7 == "TssFlnkU")
TssBiv <- subset(TAZ, TAZ$V7 == "TssBiv")
EnhBiv <- subset(TAZ, TAZ$V7 == "EnhBiv")
EnhA1 <- subset(TAZ, TAZ$V7 == "EnhA1")
EnhA2 <- subset(TAZ, TAZ$V7 == "EnhA2")
EnhWk <- subset(TAZ, TAZ$V7 == "EnhWk")
EnhG1 <- subset(TAZ, TAZ$V7 == "EnhG1")
EnhG2 <- subset(TAZ, TAZ$V7 == "EnhG2")
TxWk <- subset(TAZ, TAZ$V7 == "TxWk")
Tx <- subset(TAZ, TAZ$V7 == "Tx")
ZNF <- subset(TAZ, TAZ$V7 == "ZNF/Rpts")
ReprPCWk <- subset(TAZ, TAZ$V7 == "ReprPCWk")
ReprPC <- subset(TAZ, TAZ$V7 == "ReprPC")
Het <- subset(TAZ, TAZ$V7 == "Het")
Quies <- subset(TAZ, TAZ$V7 == "Quies")

TssA.sum <- sum(TssA$V13)
TssFlnk.sum <- sum(TssFlnk$V13)
TssFlnkD.sum <- sum(TssFlnkD$V13)
TssFlnkU.sum <- sum(TssFlnkU$V13)
TssBiv.sum <- sum(TssBiv$V13)
EnhBiv.sum <- sum(EnhBiv$V13)
EnhA1.sum <- sum(EnhA1$V13)
EnhA2.sum <- sum(EnhA2$V13)
EnhWk.sum <- sum(EnhWk$V13)
EnhG1.sum <- sum(EnhG1$V13)
EnhG2.sum <- sum(EnhG2$V13)
TxWk.sum <- sum(TxWk$V13)
Tx.sum <- sum(Tx$V13)
ZNF.sum <- sum(ZNF$V13)
ReprPCWk.sum <- sum(ReprPCWk$V13)
ReprPC.sum <- sum(ReprPC$V13)
Het.sum <- sum(Het$V13)
Quies.sum <- sum(Quies$V13)

setwd("RKO_TAZ_v_VAL_project/ChIP-seq/siQ/RKO_chromHMM_breakdown/bedtools_intersections/sample_intersections/TAZ")

#H3K27me1 signal
setwd("RKO_TAZ_v_VAL_project/ChIP-seq/siQ/consistent_peak_beds/H3K27me1/VEH/VEH1/variance_threshold_1/")

VEH1 <- fread("VEH1_peaks.variance_threshold_1.VAL_logFC_0.5.lost.bed")
VEH1 <- mutate(VEH1, width = V3-V2)


setwd("RKO_TAZ_v_VAL_project/ChIP-seq/siQ/consistent_peak_beds/H3K27me1/VEH/VEH2/variance_threshold_1/")

VEH2 <- fread("VEH2_peaks.variance_threshold_1.VAL_logFC_0.5.lost.bed")
VEH2 <- mutate(VEH2, width = V3-V2)

setwd("RKO_TAZ_v_VAL_project/ChIP-seq/siQ/consistent_peak_beds/H3K27me1/VEH/VEH3/variance_threshold_1/")

VEH3 <- fread("VEH3_peaks.variance_threshold_1.VAL_logFC_0.5.lost.bed")
VEH3 <- mutate(VEH3, width = V3-V2)


setwd("RKO_TAZ_v_VAL_project/ChIP-seq/siQ/consistent_peak_beds/H3K27me1/TAZ/TAZ1/variance_threshold_1/")

TAZ1 <- fread("TAZ1_peaks.variance_threshold_1.VAL_logFC_0.5.lost.bed")
TAZ1 <- mutate(TAZ1, width = V3-V2)


setwd("RKO_TAZ_v_VAL_project/ChIP-seq/siQ/consistent_peak_beds/H3K27me1/TAZ/TAZ2/variance_threshold_1/")

TAZ2 <- fread("TAZ2_peaks.variance_threshold_1.VAL_logFC_0.5.lost.bed")
TAZ2 <- mutate(TAZ2, width = V3-V2)

setwd("RKO_TAZ_v_VAL_project/ChIP-seq/siQ/consistent_peak_beds/H3K27me1/TAZ/TAZ3/variance_threshold_1/")

TAZ3 <- fread("TAZ3_peaks.variance_threshold_1.VAL_logFC_0.5.lost.bed")
TAZ3 <- mutate(TAZ3, width = V3-V2)

setwd("RKO_TAZ_v_VAL_project/ChIP-seq/siQ/consistent_peak_beds/H3K27me1/VAL/VAL1/variance_threshold_1/")

VAL1 <- fread("VAL1_peaks.variance_threshold_1.VAL_logFC_0.5.lost.bed")
VAL1 <- mutate(VAL1, width = V3-V2)

setwd("RKO_TAZ_v_VAL_project/ChIP-seq/siQ/consistent_peak_beds/H3K27me1/VAL/VAL2/variance_thresold_1/")

VAL2 <- fread("VAL2_peaks.variance_threshold_1.VAL_logFC_0.5.lost.bed")
VAL2 <- mutate(VAL2, width = V3-V2)

setwd("RKO_TAZ_v_VAL_project/ChIP-seq/siQ/consistent_peak_beds/H3K27me1/VAL/VAL3/variance_threshold_1/")

VAL3 <- fread("VAL3_peaks.variance_threshold_1.VAL_logFC_0.5.lost.bed")
VAL3 <- mutate(VAL3, width = V3-V2)

#overlaps of H3K27me1 called peaks
setwd("RKO_TAZ_v_VAL_project/ChIP-seq/siQ/consistent_peak_beds/H3K27me1/20250708_H3K27me1_variance_checks/0.25/H3K27me1_lost/")

#VEH
VEH1 <- fread("VEH1_peaks.variance_threshold_0.25.TAZ_logFC_0.5.lost.bed")
VEH2 <- fread("VEH2_peaks.variance_threshold_0.25.TAZ_logFC_0.5.lost.bed")

VEH1.gr <- GRanges(seqnames = VEH1$V1, ranges = IRanges(start = VEH1$V2, end = VEH1$V3))
VEH2.gr <- GRanges(seqnames = VEH2$V1, ranges = IRanges(start = VEH2$V2, end = VEH2$V3))

hits <- findOverlaps(VEH1.gr, VEH2.gr)
overlap_regions <- data.frame(pintersect(VEH1.gr[queryHits(hits)], VEH2.gr[subjectHits(hits)]))

setwd("RKO_TAZ_v_VAL_project/ChIP-seq/siQ/20250710_H3K27me1_revisited/H3K27me1_lost/")

write.table(overlap_regions, file = "VEH_called.GRanges_H3K27me1_lost.var_0.25.logFC_0.5.TAZ.bed", quote = F, sep = "\t", row.names = F, col.names = F)


#TAZ
setwd("RKO_TAZ_v_VAL_project/ChIP-seq/siQ/consistent_peak_beds/H3K27me1/20250708_H3K27me1_variance_checks/1.0/H3K27me1_lost/")

TAZ1 <- fread("TAZ1_peaks.variance_threshold_1.VAL_logFC_0.5.lost.bed")
TAZ2 <- fread("TAZ2_peaks.variance_threshold_1.VAL_logFC_0.5.lost.bed")

TAZ1.gr <- GRanges(seqnames = TAZ1$V1, ranges = IRanges(start = TAZ1$V2, end = TAZ1$V3))
TAZ2.gr <- GRanges(seqnames = TAZ2$V1, ranges = IRanges(start = TAZ2$V2, end = TAZ2$V3))

hits <- findOverlaps(TAZ1.gr, TAZ2.gr)
overlap_regions <- data.frame(pintersect(TAZ1.gr[queryHits(hits)], TAZ2.gr[subjectHits(hits)]))

setwd("RKO_TAZ_v_VAL_project/ChIP-seq/siQ/20250710_H3K27me1_revisited/H3K27me1_lost/")

write.table(overlap_regions, file = "TAZ_called.GRanges_H3K27me1_lost.var_1.logFC_0.5.VAL.bed", quote = F, sep = "\t", row.names = F, col.names = F)

#TAZ
setwd("RKO_TAZ_v_VAL_project/ChIP-seq/siQ/consistent_peak_beds/H3K27me1/20250708_H3K27me1_variance_checks/0.25/H3K27me1_gain/")

TAZ1 <- fread("TAZ1_peaks.variance_threshold_0.25.TAZ_logFC_0.5.gain.bed")
TAZ2 <- fread("TAZ2_peaks.variance_threshold_0.25.TAZ_logFC_0.5.gain.bed")

TAZ1.gr <- GRanges(seqnames = TAZ1$V1, ranges = IRanges(start = TAZ1$V2, end = TAZ1$V3))
TAZ2.gr <- GRanges(seqnames = TAZ2$V1, ranges = IRanges(start = TAZ2$V2, end = TAZ2$V3))

hits <- findOverlaps(TAZ1.gr, TAZ2.gr)
overlap_regions <- data.frame(pintersect(TAZ1.gr[queryHits(hits)], TAZ2.gr[subjectHits(hits)]))

setwd("RKO_TAZ_v_VAL_project/ChIP-seq/siQ/20250710_H3K27me1_revisited/H3K27me1_gain/")

write.table(overlap_regions, file = "TAZ_called.GRanges_H3K27me1_gain.var_0.25.logFC_0.5.TAZ.bed", quote = F, sep = "\t", row.names = F, col.names = F)


#VAL
setwd("RKO_TAZ_v_VAL_project/ChIP-seq/siQ/consistent_peak_beds/H3K27me1/20250708_H3K27me1_variance_checks/1.0/H3K27me1_gain/")

VAL1 <- fread("VAL1_peaks.variance_threshold_1.VAL_logFC_0.5.gain.bed")
VAL2 <- fread("VAL2_peaks.variance_threshold_1.VAL_logFC_0.5.gain.bed")

VAL1.gr <- GRanges(seqnames = VAL1$V1, ranges = IRanges(start = VAL1$V2, end = VAL1$V3))
VAL2.gr <- GRanges(seqnames = VAL2$V1, ranges = IRanges(start = VAL2$V2, end = VAL2$V3))

hits <- findOverlaps(VAL1.gr, VAL2.gr)
overlap_regions <- data.frame(pintersect(VAL1.gr[queryHits(hits)], VAL2.gr[subjectHits(hits)]))

setwd("RKO_TAZ_v_VAL_project/ChIP-seq/siQ/20250710_H3K27me1_revisited/H3K27me1_gain/")

write.table(overlap_regions, file = "VAL_called.GRanges_H3K27me1_gain.var_1.logFC_0.5.VAL.bed", quote = F, sep = "\t", row.names = F, col.names = F)

#sums of final called peaks
#setwd
setwd("RKO_TAZ_v_VAL_project/ChIP-seq/siQ/20250710_H3K27me1_revisited/H3K27me1_lost/")

VEH_TAZ <- fread("VEH_called.GRanges_H3K27me1_lost.var_1.logFC_0.5.TAZ.bed")
VEH_VAL <- fread("VEH_called.GRanges_H3K27me1_lost.var_1.logFC_0.5.VAL.bed")

#peak center calculation
#original_code
setwd("RKO_TAZ_v_VAL_project/ChIP-seq/siQ/20250710_H3K27me1_revisited/H3K27me1_lost/")

sample <- fread("TAZ_called.GRanges_H3K27me1_lost.var_0.25.logFC_0.5.TAZ.bed", header = F)

sample <- mutate(sample, width = V3-V2)
sample_center <- mutate(sample, sample_center_start = (width/2)+V2)
sample_center <- mutate(sample_center, sample_center_start_round = round(sample_center_start, digits = 0))
sample_center <- mutate(sample_center, sample_center_start_stop = sample_center_start_round + 1)

sample_center <- sample_center[,c(1,9:10,4)]

setwd("RKO_TAZ_v_VAL_project/ChIP-seq/siQ/20250710_H3K27me1_revisited/H3K27me1_lost/")
write.table(sample_center, file = "sample.peak_center.bed", quote = F, row.names = F, col.names = F, sep = "\t")

# List all .bed files you want to process
setwd("RKO_TAZ_v_VAL_project/ChIP-seq/siQ/20250710_H3K27me1_revisited/H3K27me1_gain/")
setwd("RKO_TAZ_v_VAL_project/ChIP-seq/siQ/20250710_H3K27me1_revisited/Variance_1.final_H3K27me1_peaks/")
setwd("RKO_TAZ_v_VAL_project/ATAC-seq/20261008_10kb_K27me1_integration/")

bed_files <- list.files(pattern = "\\.bed$")

# Loop over each BED file
for (file in bed_files) {
  
  # Read file
  sample <- fread(file, header = FALSE)
  
  # Compute center position
  sample <- sample %>%
    mutate(width = V3 - V2) %>%
    mutate(sample_center_start = (width / 2) + V2) %>%
    mutate(sample_center_start_round = round(sample_center_start, digits = 0)) %>%
    mutate(sample_center_start_stop = sample_center_start_round + 1)
  
  # Keep desired columns: chr, center start, center stop, original 4th column
  sample_center <- sample[, c(1, 7:8, 4,5)]
  
  # Generate output filename based on input filename
  outname <- gsub(".bed$", ".peak_center.bed", file)
  
  # Write output
  write.table(sample_center, file = outname, quote = FALSE,
              row.names = FALSE, col.names = FALSE, sep = "\t")
}

#Calculating coverages of different variances
setwd("RKO_TAZ_v_VAL_project/ChIP-seq/siQ/20250710_H3K27me1_revisited/H3K27me1_lost/GMM_clustering/TAZ_called/")

bed_files <- list.files(pattern = "\\.bed$")

# Initialize an empty data frame to store results
results <- data.frame(filename = character(), sum_V4 = numeric(), stringsAsFactors = FALSE)

# Loop over each file
for (file in bed_files) {
  sample <- fread(file, header = FALSE)
  sum_val <- sum(sample$V4, na.rm = TRUE)  # safer with na.rm
  results <- rbind(results, data.frame(filename = file, sum_V4 = sum_val))
}

# Write the results to a file
write.table(results, file = "sum_width_summary.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

#ChromHMM_coverage_for H3K27me1 revised
setwd("RKO_TAZ_v_VAL_project/ChIP-seq/siQ/RKO_chromHMM_breakdown/chromHMM_intersections_for_H3K27ac_cluster_5mCwH3K27me1/")

# List all .bed files
bed_files <- list.files(pattern = "\\chromHMM.bed$")

# Loop through each file
for (file in bed_files) {
  message("Processing ", file)
  
  # Load data
  Sample <- fread(file, header = FALSE)
  
  TssA <- subset(Sample, Sample$V7 == "TssA")
  TssFlnk <- subset(Sample, Sample$V7 == "TssFlnk")
  TssFlnkD <- subset(Sample, Sample$V7 == "TssFlnkD")
  TssFlnkU <- subset(Sample, Sample$V7 == "TssFlnkU")
  TssBiv <- subset(Sample, Sample$V7 == "TssBiv")
  EnhBiv <- subset(Sample, Sample$V7 == "EnhBiv")
  EnhA1 <- subset(Sample, Sample$V7 == "EnhA1")
  EnhA2 <- subset(Sample, Sample$V7 == "EnhA2")
  EnhWk <- subset(Sample, Sample$V7 == "EnhWk")
  EnhG1 <- subset(Sample, Sample$V7 == "EnhG1")
  EnhG2 <- subset(Sample, Sample$V7 == "EnhG2")
  TxWk <- subset(Sample, Sample$V7 == "TxWk")
  Tx <- subset(Sample, Sample$V7 == "Tx")
  ZNF <- subset(Sample, Sample$V7 == "ZNF/Rpts")
  ReprPCWk <- subset(Sample, Sample$V7 == "ReprPCWk")
  ReprPC <- subset(Sample, Sample$V7 == "ReprPC")
  Het <- subset(Sample, Sample$V7 == "Het")
  Quies <- subset(Sample, Sample$V7 == "Quies")
  
  TssA.sum <- sum(TssA$V13)
  TssFlnk.sum <- sum(TssFlnk$V13)
  TssFlnkD.sum <- sum(TssFlnkD$V13)
  TssFlnkU.sum <- sum(TssFlnkU$V13)
  TssBiv.sum <- sum(TssBiv$V13)
  EnhBiv.sum <- sum(EnhBiv$V13)
  EnhA1.sum <- sum(EnhA1$V13)
  EnhA2.sum <- sum(EnhA2$V13)
  EnhWk.sum <- sum(EnhWk$V13)
  EnhG1.sum <- sum(EnhG1$V13)
  EnhG2.sum <- sum(EnhG2$V13)
  TxWk.sum <- sum(TxWk$V13)
  Tx.sum <- sum(Tx$V13)
  ZNF.sum <- sum(ZNF$V13)
  ReprPCWk.sum <- sum(ReprPCWk$V13)
  ReprPC.sum <- sum(ReprPC$V13)
  Het.sum <- sum(Het$V13)
  Quies.sum <- sum(Quies$V13)

# Combine results into a named vector
state_sums <- c(
  TssA = TssA.sum,
  TssFlnk = TssFlnk.sum,
  TssFlnkD = TssFlnkD.sum,
  TssFlnkU = TssFlnkU.sum,
  TssBiv = TssBiv.sum,
  EnhBiv = EnhBiv.sum,
  EnhA1 = EnhA1.sum,
  EnhA2 = EnhA2.sum,
  EnhWk = EnhWk.sum,
  EnhG1 = EnhG1.sum,
  EnhG2 = EnhG2.sum,
  TxWk = TxWk.sum,
  Tx = Tx.sum,
  ZNF = ZNF.sum,
  ReprPCWk = ReprPCWk.sum,
  ReprPC = ReprPC.sum,
  Het = Het.sum,
  Quies = Quies.sum
)

# Convert to data frame
state_sums_df <- data.frame(State = names(state_sums), Sum = as.numeric(state_sums))

# Write to file (TSV)
write.table(state_sums_df, file = "state_sums.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
#save file

base <- sub("\\chromHMM.bed$", "", file)

write.table(state_sums_df, file = paste0(base, ".state_sum.tsv"),
            quote = FALSE, sep = "\t", col.names = FALSE, row.names = FALSE)
}











