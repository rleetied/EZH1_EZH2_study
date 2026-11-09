#!/usr/bin/env bash

# Reconstructed from the production HPC scripts. Replace the example variables
# and paths with the appropriate sample-specific values before execution.

set -euo pipefail

SAMPLE="A1RKOVeh"
R1="${SAMPLE}_L000_R1_001.fastq.gz"
R2="${SAMPLE}_L000_R2_001.fastq.gz"
TRIM_DIR="trimmed_fastq"
HG38_INDEX="/path/to/bowtie_index/hg38.SNAP"
BLACKLIST="/path/to/hg38-blacklist.v2.bed"

# 1. Adapter/quality trimming and post-trim FastQC.
trim_galore \
  --quality 25 \
  --fastqc \
  --output_dir "${TRIM_DIR}" \
  --paired "${R1}" "${R2}"

# 2. Paired-end alignment; duplicate removal; retain primary, properly paired
#    alignments with MAPQ >=30; coordinate sort.
bowtie2 \
  --very-sensitive \
  -I 0 \
  -X 1000 \
  -p 4 \
  -x "${HG38_INDEX}" \
  -1 "${TRIM_DIR}/${SAMPLE}_L000_R1_001_val_1.fq.gz" \
  -2 "${TRIM_DIR}/${SAMPLE}_L000_R2_001_val_2.fq.gz" | \
samblaster --removeDups | \
samtools view -@ 4 -b -q 30 -f 2 -F 2828 | \
samtools sort -@ 4 -o "${SAMPLE}_PE.sorted.dedup.bt2.bam" -

# 3. Remove mitochondrial alignments and index.
samtools view -bh -e 'rname != "chrM"' \
  "${SAMPLE}_PE.sorted.dedup.bt2.bam" \
  > "${SAMPLE}_PE.sorted.dedup.noMT.bt2.bam"
samtools index "${SAMPLE}_PE.sorted.dedup.noMT.bt2.bam"

# 4. Sample-level peak calling. No --call-summits option was used.
macs3 callpeak \
  --treatment "${SAMPLE}_PE.sorted.dedup.noMT.bt2.bam" \
  --format BAMPE \
  --gsize hs \
  --name "${SAMPLE}" \
  --outdir macs3_no_multisummit \
  --qvalue 0.01 \
  --keep-dup all \
  --bdg \
  --SPMR

# 5. Remove blacklist overlaps and retain canonical chromosomes.
bedtools intersect -v \
  -a "macs3_no_multisummit/${SAMPLE}_peaks.narrowPeak" \
  -b "${BLACKLIST}" | \
awk 'BEGIN {OFS="\t"} $1 ~ /^chr([1-9]|1[0-9]|2[0-2]|X|Y)$/' \
  > "${SAMPLE}_peaks.clean.narrowPeak"

# 6. Treatment-level reproducibility (example for three VEH replicates):
#    retain segments represented in >=2/3 replicates, then merge adjacent ones.
bedtools multiinter \
  -i A1RKOVeh_peaks.clean.narrowPeak \
     B7RKOVeh_peaks.clean.narrowPeak \
     C13RKOVeh_peaks.clean.narrowPeak | \
awk 'BEGIN {OFS="\t"} $4 >= 2 {print $1,$2,$3}' | \
bedtools merge -i - \
  > Veh_reproducible.bed

# Repeat the reproducibility step for TAZ, VAL, DAC, TD, and VD. Combine and
# merge the six treatment-level BEDs to create the final consensus universe.
cat Veh_reproducible.bed TAZ_reproducible.bed VAL_reproducible.bed \
    DAC_reproducible.bed TAZ_DAC_reproducible.bed VAL_DAC_reproducible.bed | \
sort -k1,1 -k2,2n | \
bedtools merge -i - \
  > ATAC_consensus_peaks.bed

# Convert the consensus BED to SAF (BED starts are zero-based; SAF is one-based).
awk 'BEGIN {OFS="\t"; print "GeneID","Chr","Start","End","Strand"}
     {printf "ATAC_peak_%06d\t%s\t%d\t%d\t.\n", NR,$1,$2+1,$3}' \
  ATAC_consensus_peaks.bed \
  > ATAC_consensus_peaks.saf

# 7. Count paired-end fragments in the consensus peak universe.
featureCounts \
  -T 8 \
  -F SAF \
  -a ATAC_consensus_peaks.saf \
  -o ATAC_consensus_featureCounts.txt \
  -p --countReadPairs -B -C --primary -Q 30 -s 0 \
  *_PE.sorted.dedup.noMT.bt2.bam
