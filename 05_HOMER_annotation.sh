#!/bin/bash

# Set the directory where your .bed files live
BED_DIR="/home/rochelle.tiedemann/rothbart-secondary/Rochelle/20240808_RKO_TAZ_v_VAL_project/ATAC-seq/ATAC_VD_dominant"

# List all .bed files (without the full path, just filenames)
BED_FILES=($(basename -a "$BED_DIR"/*.bed))

# Loop through each BED file
for i in "${!BED_FILES[@]}"; do
    BED_NAME="${BED_FILES[$i]}"
    BASENAME="${BED_NAME%.bed}"  # Remove .bed extension

    cat > slurm.HOMER.$i <<EOF
#!/bin/bash
#SBATCH --job-name=HOMER_${BASENAME}
#SBATCH --ntasks=1
#SBATCH --time=10:00:00

# Activate conda environment
source activate /home/rochelle.tiedemann/rothbart-primary/Rochelle/anaconda2/envs/Oasis

# Change to BED file directory
cd $BED_DIR

# Run HOMER motif finding
findMotifsGenome.pl ${BED_NAME} hg38 ${BASENAME}.hg38_motif/ -len 8,10,12 -p 4

# Annotate peaks
annotatePeaks.pl ${BED_NAME} hg38 -CpG -go ${BASENAME}.GO_analysis/ > ${BASENAME}.hg38.bed
EOF

done
