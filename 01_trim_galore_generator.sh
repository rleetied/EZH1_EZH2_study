#!/bin/bash
set -euxo pipefail

#Build up a space sep'd array of sample/file prefix
declare -a SAMPLE=(HCT116DACA HCT116DACB HCT116DACC HCT116VALA HCT116VALB HCT116VALC HCT116VDA HCT116VDB HCT116VDC HCT116VehA HCT116VehB HCT116VehC SW620DACA SW620DACB SW620DACC SW620TAZA SW620TAZB SW620TAZC SW620TDA SW620TDB SW620TDC SW620VALA SW620VALB SW620VALC SW620VDA SW620VDB SW620VDC SW620VehA SW620VehB SW620VehC)

#Adjust the length of the array to iterate over to match the number of samples/files to analyze
for i in {0..29}; do
    cat > slurm.trim_galore.$i <<EOF
#!/bin/bash
#
#SBATCH --job-name=trimgalore
#SBATCH --ntasks=1
#SBATCH --time=168:00:00

source activate /home/rochelle.tiedemann/rothbart-primary/Rochelle/anaconda2/envs/Oasis2

#Change to raw data directory

cd /home/rochelle.tiedemann/rothbart-secondary/Rochelle/20240808_RKO_TAZ_v_VAL_project/PR001634_ROTS_RNA

#Launch Trim Galore!
trim_galore --quality 25 --fastqc --output_dir /home/rochelle.tiedemann/rothbart-secondary/Rochelle/20240808_RKO_TAZ_v_VAL_project/PR001634_ROTS_RNA/trimmed_fastq --paired ${SAMPLE[i]}_L000_R1_001.fastq.gz ${SAMPLE[i]}_L000_R2_001.fastq.gz

source deactivate /home/rochelle.tiedemann/rothbart-primary/Rochelle/anaconda2/envs/Oasis2

EOF
done
