#!/bin/bash

# Set the directory where your .bed files live

BED_DIR="/home/rochelle.tiedemann/rothbart-secondary/Rochelle/20240808_RKO_TAZ_v_VAL_project/ChIP-seq/siQ/20250710_H3K27me1_revisited/H3K27me1_lost/peak_center/Variance_1_LogFC_0.5"

# List all .bed files (without the full path, just filenames) 
BED_FILES=($(basename -a "$BED_DIR"/*.bed))

# Loop through each BED file

for i in "${!BED_FILES[@]}"; do
    BED_NAME="${BED_FILES[$i]}"
    BASENAME="${BED_NAME%.bed}"  # Remove .bed extension
    
    cat > slurm.compute.K27me1_signal.$i <<EOF
#!/bin/bash
#
#SBATCH --job-name=computeMatrix
#SBATCH --ntasks=1
#SBATCH --time=128:00:00

#Change to trimmed data directory

cd /home/rochelle.tiedemann/rothbart-secondary/Rochelle/20240808_RKO_TAZ_v_VAL_project/Integration/ChIP_seq_peaks/

source activate /home/rochelle.tiedemann/rothbart-secondary/Rochelle/miniconda3/envs/Oasis3

#Launch deeptools
computeMatrix reference-point \
-S /home/rochelle.tiedemann/rothbart-secondary/Rochelle/20240808_RKO_TAZ_v_VAL_project/ChIP-seq/siQ/siQ-ChIP/bigwigs/VEH1me1_siQ.bw \
/home/rochelle.tiedemann/rothbart-secondary/Rochelle/20240808_RKO_TAZ_v_VAL_project/ChIP-seq/siQ/siQ-ChIP/bigwigs/VEH2me1_siQ.bw \
/home/rochelle.tiedemann/rothbart-secondary/Rochelle/20240808_RKO_TAZ_v_VAL_project/ChIP-seq/siQ/siQ-ChIP/bigwigs/VEH3me1_siQ.bw \
/home/rochelle.tiedemann/rothbart-secondary/Rochelle/20240808_RKO_TAZ_v_VAL_project/ChIP-seq/siQ/siQ-ChIP/bigwigs/TAZ1me1_siQ.bw \
/home/rochelle.tiedemann/rothbart-secondary/Rochelle/20240808_RKO_TAZ_v_VAL_project/ChIP-seq/siQ/siQ-ChIP/bigwigs/TAZ2me1_siQ.bw \
/home/rochelle.tiedemann/rothbart-secondary/Rochelle/20240808_RKO_TAZ_v_VAL_project/ChIP-seq/siQ/siQ-ChIP/bigwigs/TAZ3me1_siQ.bw \
/home/rochelle.tiedemann/rothbart-secondary/Rochelle/20240808_RKO_TAZ_v_VAL_project/ChIP-seq/siQ/siQ-ChIP/bigwigs/VAL1me1_siQ.bw \
/home/rochelle.tiedemann/rothbart-secondary/Rochelle/20240808_RKO_TAZ_v_VAL_project/ChIP-seq/siQ/siQ-ChIP/bigwigs/VAL2me1_siQ.bw \
/home/rochelle.tiedemann/rothbart-secondary/Rochelle/20240808_RKO_TAZ_v_VAL_project/ChIP-seq/siQ/siQ-ChIP/bigwigs/VAL3me1_siQ.bw \
-R ${BED_DIR}/${BED_NAME} \
--beforeRegionStartLength 10000 --afterRegionStartLength 10000 \
--outFileName ${BASENAME}.H3K27me1_signal.10k.tab \
--skipZeros \
--outFileNameMatrix ${BASENAME}.H3K27me1_signal.10k.tab.gz

EOF
done
