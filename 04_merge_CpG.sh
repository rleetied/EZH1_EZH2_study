#!/bin/bash
set -euxo pipefail

declare -a SAMPLE=(
  acVAL1_combined
  acVAL2_combined
  acVEH1_combined
  acVEH2_combined
)

WORKDIR="/home/rochelle.tiedemann/rothbart-secondary/Rochelle/PR001883_ROTS_EMChIP/trimmed_fastq/bams/H3K27ac"
REF="/home/rochelle.tiedemann/rothbart-secondary/Rochelle/alignment_indexes/hg38_index/biscuit_index/hg38.final.fa"

for i in "${!SAMPLE[@]}"; do
  cat > slurm.biscuitscript.$i <<EOF
#!/bin/bash
#SBATCH --job-name=biscuit_${SAMPLE[i]}
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --time=168:00:00

cd $WORKDIR

module load bbc2/samblaster
module load bbc2/pigz
module load bbc2/parallel
module load bbc2/samtools
module load bbc2/dupsifter

source activate biscuit

biscuit vcf2bed \\
-t cg \\
${SAMPLE[i]}.vcf.gz > ${SAMPLE[i]}.bed

biscuit mergecg $REF ${SAMPLE[i]}.bed > ${SAMPLE[i]}.mergecg.bed

EOF
done
