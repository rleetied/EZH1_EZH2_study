#!/usr/bin/env bash
# Detailed production commands: ./processing_code/.

# ATAC-seq preprocessing summary
# Final use: Figures 3F-I and 4A-D

PROJECT_ROOT="RKO_TAZ_v_VAL_project"
FASTQ_DIR="${PROJECT_ROOT}/ATAC-seq/raw_fastq"
TRIM_DIR="${PROJECT_ROOT}/ATAC-seq/trimmed_fastq"
HG38_INDEX="<path-to-bowtie2-hg38-index>"
BLACKLIST="<path-to-hg38-blacklist.bed>"
SAMPLES=("<18-sample-prefixes>")

for SAMPLE in "${SAMPLES[@]}"; do
  trim_galore --quality 25 --fastqc --paired \
    --output_dir "${TRIM_DIR}" \
    "${FASTQ_DIR}/${SAMPLE}_R1.fastq.gz" \
    "${FASTQ_DIR}/${SAMPLE}_R2.fastq.gz"

  # Align; remove duplicate fragments; retain properly paired, primary
  # alignments with MAPQ >=30; sort; remove chrM.
  bowtie2 --very-sensitive -I 0 -X 1000 -p 4 -x "${HG38_INDEX}" \
    -1 "${TRIM_DIR}/${SAMPLE}_R1_val_1.fq.gz" \
    -2 "${TRIM_DIR}/${SAMPLE}_R2_val_2.fq.gz" | \
    samblaster --removeDups | \
    samtools view -b -q 30 -f 2 -F 2828 | \
    samtools sort -o "${SAMPLE}.sorted.dedup.bam" -

  samtools view -bh -e 'rname != "chrM"' \
    "${SAMPLE}.sorted.dedup.bam" > "${SAMPLE}.final.bam"
  samtools index "${SAMPLE}.final.bam"

  # Sample-level accessible-region calling; multisummit splitting was not used.
  macs3 callpeak -t "${SAMPLE}.final.bam" -f BAMPE -g hs \
    -n "${SAMPLE}" --qvalue 0.01 --keep-dup all --bdg --SPMR

  bedtools intersect -v -a "${SAMPLE}_peaks.narrowPeak" -b "${BLACKLIST}" \
    > "${SAMPLE}.noBlacklist.narrowPeak"
done

# Within each treatment, retain genomic segments supported in >=2 of 3
# replicates using bedtools multiinter, then merge adjacent segments.
# Merge the six treatment-level sets to form the final 88,726-peak consensus.

# Count paired-end fragments in consensus intervals with featureCounts v2.1.1.
featureCounts -T 8 -F SAF -a ATAC_consensus_peaks.saf \
  -o ATAC_consensus_featureCounts.txt \
  -p --countReadPairs -B -C --primary -Q 30 -s 0 *.final.bam

