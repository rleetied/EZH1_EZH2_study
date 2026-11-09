#!/usr/bin/env bash
# Detailed production commands: ./processing_code/.

# H3K27ac ChIP-EM-seq (acEM-seq) processing summary
# Final use: Figure 4F-G

PROJECT_ROOT="RKO_TAZ_v_VAL_project"
FASTQ_DIR="${PROJECT_ROOT}/ChIP_EMseq/raw_fastq"
TRIM_DIR="${PROJECT_ROOT}/ChIP_EMseq/trimmed_fastq"
BISCUIT_REF="<path-to-biscuit-indexed-hg38.fa>"
SAMPLES=(acVEH1 acVEH2 acVAL1 acVAL2)

for SAMPLE in "${SAMPLES[@]}"; do
  trim_galore --quality 25 --fastqc --paired --output_dir "${TRIM_DIR}" \
    "${FASTQ_DIR}/${SAMPLE}_R1_001.fastq.gz" \
    "${FASTQ_DIR}/${SAMPLE}_R2_001.fastq.gz"

  biscuit align -@ 30 -R "@RG\tID:${SAMPLE}\tSM:${SAMPLE}\tPL:Illumina" \
    "${BISCUIT_REF}" \
    "${TRIM_DIR}/${SAMPLE}_R1_001_val_1.fq.gz" \
    "${TRIM_DIR}/${SAMPLE}_R2_001_val_2.fq.gz" | \
    dupsifter "${BISCUIT_REF}" | \
    samtools sort -o "${SAMPLE}.bam" -
  samtools index "${SAMPLE}.bam"

  biscuit pileup -@ 30 -o "${SAMPLE}.vcf" \
    "${BISCUIT_REF}" "${SAMPLE}.bam"
  bgzip "${SAMPLE}.vcf"
  tabix -p vcf "${SAMPLE}.vcf.gz"

  biscuit vcf2bed -t cg "${SAMPLE}.vcf.gz" > "${SAMPLE}.cg.bed"
  biscuit mergecg "${BISCUIT_REF}" "${SAMPLE}.cg.bed" \
    > "${SAMPLE}.mergecg.bed"
done

# CpG calls required at least three reads. Methylation tracks were intersected
# with final ATAC regions and summarized in bins centered on ATAC peak centers.

