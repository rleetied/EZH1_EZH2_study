#!/usr/bin/env bash
# Detailed production commands: ./processing_code/.

# Multi-omic matrix and annotation processing summary
# Final use: Figure 4C-J

PROJECT_ROOT="RKO_TAZ_v_VAL_project"
REGIONS="${PROJECT_ROOT}/Integration/final_region_sets"
OUT="${PROJECT_ROOT}/Integration/final_matrices"

# Region sets represent comparable TD/VD and VD-favored accessibility gains,
# separated by TAZ H3K27me1 gain or TAZ-maintained H3K27me1.

computeMatrix reference-point -R "${REGIONS}"/*.bed \
  -S "${PROJECT_ROOT}/ATAC-seq/bigwigs"/*.bw \
  --beforeRegionStartLength 5000 --afterRegionStartLength 5000 \
  --skipZeros -o "${OUT}/ATAC.matrix.gz"

computeMatrix reference-point -R "${REGIONS}"/*.bed \
  -S "${PROJECT_ROOT}/ChIP-seq/siQ/bigwigs"/*.bw \
  --beforeRegionStartLength 5000 --afterRegionStartLength 5000 \
  --skipZeros -o "${OUT}/histone_marks.matrix.gz"

computeMatrix reference-point -R "${REGIONS}"/*.bed \
  -S "${PROJECT_ROOT}/ChIP_EMseq/bedGraphs"/*.bedGraph \
  --beforeRegionStartLength 1500 --afterRegionStartLength 1500 \
  --binSize 10 --missingDataAsZero \
  -o "${OUT}/acEMseq.matrix.gz"

# Motif/genomic annotation was performed when required with HOMER.
for BED in "${REGIONS}"/*.bed; do
  NAME="$(basename "${BED}" .bed)"
  findMotifsGenome.pl "${BED}" hg38 "${OUT}/${NAME}.motifs" -len 8,10,12 -p 4
  annotatePeaks.pl "${BED}" hg38 > "${OUT}/${NAME}.annotated.txt"
done

