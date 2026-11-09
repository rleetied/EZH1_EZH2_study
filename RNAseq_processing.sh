#!/usr/bin/env bash
# Detailed production commands: ./processing_code/.

# RNA-seq preprocessing summary
# Final use: Figure 1G-J

PROJECT_ROOT="RKO_TAZ_v_VAL_project"
FASTQ_DIR="${PROJECT_ROOT}/RNA-seq/raw_fastq"
TRIM_DIR="${PROJECT_ROOT}/RNA-seq/trimmed_fastq"
STAR_INDEX="<path-to-hg38-STAR-index>"

# Samples comprised three replicates of VEH, DAC, TAZ, VAL, TD (DAC+TAZ),
# and VD (DAC+VAL) for RKO and SW620 cells.
SAMPLES=("<sample-prefixes>")

for SAMPLE in "${SAMPLES[@]}"; do
  # Adapter/quality trimming and post-trim FastQC.
  trim_galore \
    --quality 25 \
    --fastqc \
    --paired \
    --output_dir "${TRIM_DIR}" \
    "${FASTQ_DIR}/${SAMPLE}_L000_R1_001.fastq.gz" \
    "${FASTQ_DIR}/${SAMPLE}_L000_R2_001.fastq.gz"

  # Two-pass alignment to hg38 and gene-level counting with STAR.
  STAR \
    --genomeDir "${STAR_INDEX}" \
    --readFilesIn \
      "${TRIM_DIR}/${SAMPLE}_L000_R1_001_val_1.fq.gz" \
      "${TRIM_DIR}/${SAMPLE}_L000_R2_001_val_2.fq.gz" \
    --readFilesCommand zcat \
    --twopassMode Basic \
    --runThreadN 8 \
    --outReadsUnmapped None \
    --sjdbOverhang 50 \
    --quantMode GeneCounts \
    --outFileNamePrefix "${TRIM_DIR}/${SAMPLE}"
done

