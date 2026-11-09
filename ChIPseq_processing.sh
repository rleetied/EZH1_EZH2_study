#!/usr/bin/env bash
# Detailed production commands: ./processing_code/.

# ChIP-seq processing summary for H3K27me1, H3K27me3, and H3K27ac
# Final use: Figures 3 and 4E-F

PROJECT_ROOT="RKO_TAZ_v_VAL_project"
FASTQ_DIR="${PROJECT_ROOT}/ChIP-seq/raw_fastq"
HG38_INDEX="<path-to-bowtie2-hg38-index>"

for R1 in "${FASTQ_DIR}"/*_R1_001.fastq.gz; do
  R2="${R1/_R1_/_R2_}"
  SAMPLE="$(basename "${R1}" _R1_001.fastq.gz)"

  bowtie2 -p 20 -I 0 -X 700 --end-to-end --sensitive \
    -x "${HG38_INDEX}" -1 "${R1}" -2 "${R2}" -S "${SAMPLE}.sam"

  # Retain MAPQ >=20 properly oriented paired fragments and convert to BED.
  awk -v MAQ=20 '($5>=MAQ && $2==99) || ($5>=MAQ && $2==163) \
    {print $3"\t"$4"\t"$4+$9-1}' "${SAMPLE}.sam" | \
    awk '$2<=$3 {print $1"\t"$2"\t"$3}' | \
    sort -k1,1 -k2,2n > "${SAMPLE}.bed"
done

# Spike-in quantitative ChIP-seq signal tracks and called regions were
# generated for each histone modification. Final analyses used H3K27me1
# regions with replicate variance <=1 and absolute log2FC >=0.5.

# deepTools matrices were calculated around called H3K27me1 regions and final
# ATAC peak centers. Representative command:
computeMatrix reference-point \
  -S "<normalized-bigWig-files>" \
  -R "<final-region-BED-files>" \
  --beforeRegionStartLength 10000 \
  --afterRegionStartLength 10000 \
  --skipZeros \
  -o "<output-matrix.gz>"

