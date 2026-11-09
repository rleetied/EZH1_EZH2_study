#!/bin/bash

# Set the directory where your .bed files live
BED_DIR="/home/rochelle.tiedemann/rothbart-secondary/Rochelle/20240808_RKO_TAZ_v_VAL_project/ATAC-seq/20261008_10kb_K27me1_integration/peak_center"
BW_DIR="/home/rochelle.tiedemann/rothbart-secondary/Rochelle/20240808_RKO_TAZ_v_VAL_project/ChIP_EMseq/bigwigs"

# List all .bed files (without the full path, just filenames)
BED_FILES=($(basename -a "$BED_DIR"/*.bed))

# Loop through each BED file
for i in "${!BED_FILES[@]}"; do
    BED_NAME="${BED_FILES[$i]}"
    BASENAME="${BED_NAME%.bed}"  # Remove .bed extension
    cat > slurm.compute.ac_EMsignal.ATAC_regions.$i <<EOF
#!/bin/bash
#
#SBATCH --job-name=computeMatrix
#SBATCH --ntasks=1
#SBATCH --time=10:00:00

#Change to trimmed data directory

cd /home/rochelle.tiedemann/rothbart-secondary/Rochelle/20240808_RKO_TAZ_v_VAL_project/Integration/

source activate /home/rochelle.tiedemann/rothbart-secondary/Rochelle/miniconda3/envs/Oasis3

#Launch deeptools
computeMatrix reference-point \
-S ${BW_DIR}/acVEH1.CpG.cov_3.bw \
${BW_DIR}/acVEH2.CpG.cov_3.bw \
${BW_DIR}/acVAL1.CpG.cov_3.bw \
${BW_DIR}/acVAL2.CpG.cov_3.bw \
-R ${BED_DIR}/${BED_NAME} \
--beforeRegionStartLength 1500 --afterRegionStartLength 1500 \
--outFileName ${BASENAME}.ATAC.acEMseq_5mC_signal.reference_point.1.5k.tab \
--skipZeros \
--outFileNameMatrix ${BASENAME}.ATAC.acEMseq_5mC_signal.reference_point.1.5k.tab.gz

EOF
done
