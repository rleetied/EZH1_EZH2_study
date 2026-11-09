#!/bin/bash

# Set the directory where your .bed files live

BED_DIR="/home/rochelle.tiedemann/rothbart-secondary/Rochelle/20240808_RKO_TAZ_v_VAL_project/ATAC-seq/20261008_10kb_K27me1_integration/peak_center"

# List all .bed files (without the full path, just filenames) 
BED_FILES=($(basename -a "$BED_DIR"/*.bed))

# Loop through each BED file

for i in "${!BED_FILES[@]}"; do
    BED_NAME="${BED_FILES[$i]}"
    BASENAME="${BED_NAME%.bed}"  # Remove .bed extension
    
    cat > slurm.compute.atac_signal.$i <<EOF
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
-S /home/rochelle.tiedemann/rothbart-secondary/Rochelle/20240808_RKO_TAZ_v_VAL_project/ATAC-seq/bigwigs/A1RKOVeh_PE.sorted.dedup.noMT.bt2.fragments.cpm.50.bw \
/home/rochelle.tiedemann/rothbart-secondary/Rochelle/20240808_RKO_TAZ_v_VAL_project/ATAC-seq/bigwigs/B7RKOVeh_PE.sorted.dedup.noMT.bt2.fragments.cpm.50.bw \
/home/rochelle.tiedemann/rothbart-secondary/Rochelle/20240808_RKO_TAZ_v_VAL_project/ATAC-seq/bigwigs/C13RKOVeh_PE.sorted.dedup.noMT.bt2.fragments.cpm.50.bw \
/home/rochelle.tiedemann/rothbart-secondary/Rochelle/20240808_RKO_TAZ_v_VAL_project/ATAC-seq/bigwigs/A4RKOD100_PE.sorted.dedup.noMT.bt2.fragments.cpm.50.bw \
/home/rochelle.tiedemann/rothbart-secondary/Rochelle/20240808_RKO_TAZ_v_VAL_project/ATAC-seq/bigwigs/B10RKOD100_PE.sorted.dedup.noMT.bt2.fragments.cpm.50.bw \
/home/rochelle.tiedemann/rothbart-secondary/Rochelle/20240808_RKO_TAZ_v_VAL_project/ATAC-seq/bigwigs/C16RKOD100_PE.sorted.dedup.noMT.bt2.fragments.cpm.50.bw \
/home/rochelle.tiedemann/rothbart-secondary/Rochelle/20240808_RKO_TAZ_v_VAL_project/ATAC-seq/bigwigs/A2RKOTAZ_PE.sorted.dedup.noMT.bt2.fragments.cpm.50.bw \
/home/rochelle.tiedemann/rothbart-secondary/Rochelle/20240808_RKO_TAZ_v_VAL_project/ATAC-seq/bigwigs/B8RKOTAZ_PE.sorted.dedup.noMT.bt2.fragments.cpm.50.bw \
/home/rochelle.tiedemann/rothbart-secondary/Rochelle/20240808_RKO_TAZ_v_VAL_project/ATAC-seq/bigwigs/C14RKOTAZ_PE.sorted.dedup.noMT.bt2.fragments.cpm.50.bw \
/home/rochelle.tiedemann/rothbart-secondary/Rochelle/20240808_RKO_TAZ_v_VAL_project/ATAC-seq/bigwigs/A3RKOVAL_PE.sorted.dedup.noMT.bt2.fragments.cpm.50.bw \
/home/rochelle.tiedemann/rothbart-secondary/Rochelle/20240808_RKO_TAZ_v_VAL_project/ATAC-seq/bigwigs/B9RKOVAL_PE.sorted.dedup.noMT.bt2.fragments.cpm.50.bw \
/home/rochelle.tiedemann/rothbart-secondary/Rochelle/20240808_RKO_TAZ_v_VAL_project/ATAC-seq/bigwigs/C15RKOVAL_PE.sorted.dedup.noMT.bt2.fragments.cpm.50.bw \
/home/rochelle.tiedemann/rothbart-secondary/Rochelle/20240808_RKO_TAZ_v_VAL_project/ATAC-seq/bigwigs/A5RKOTD_PE.sorted.dedup.noMT.bt2.fragments.cpm.50.bw \
/home/rochelle.tiedemann/rothbart-secondary/Rochelle/20240808_RKO_TAZ_v_VAL_project/ATAC-seq/bigwigs/B11RKOTD_PE.sorted.dedup.noMT.bt2.fragments.cpm.50.bw \
/home/rochelle.tiedemann/rothbart-secondary/Rochelle/20240808_RKO_TAZ_v_VAL_project/ATAC-seq/bigwigs/C17RKOTD_PE.sorted.dedup.noMT.bt2.fragments.cpm.50.bw \
/home/rochelle.tiedemann/rothbart-secondary/Rochelle/20240808_RKO_TAZ_v_VAL_project/ATAC-seq/bigwigs/A6RKOVD_PE.sorted.dedup.noMT.bt2.fragments.cpm.50.bw \
/home/rochelle.tiedemann/rothbart-secondary/Rochelle/20240808_RKO_TAZ_v_VAL_project/ATAC-seq/bigwigs/B12RKOVD_PE.sorted.dedup.noMT.bt2.fragments.cpm.50.bw \
/home/rochelle.tiedemann/rothbart-secondary/Rochelle/20240808_RKO_TAZ_v_VAL_project/ATAC-seq/bigwigs/C18RKOVD_PE.sorted.dedup.noMT.bt2.fragments.cpm.50.bw \
-R ${BED_DIR}/${BED_NAME} \
--beforeRegionStartLength 10000 --afterRegionStartLength 10000 \
--outFileName ${BASENAME}.ATAC_signal.10k.tab \
--skipZeros \
--outFileNameMatrix ${BASENAME}.ATAC_signal.10k.tab.gz

EOF
done
