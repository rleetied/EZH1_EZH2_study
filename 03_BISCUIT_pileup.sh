#!/bin/bash
#SBATCH --job-name=biscuit_${SAMPLE[i]}
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=30
#SBATCH --time=168:00:00

cd /home/rochelle.tiedemann/rothbart-secondary/Rochelle/PR001883_ROTS_EMChIP/trimmed_fastq/bams/H3K27ac

module load bbc2/samblaster
module load bbc2/pigz
module load bbc2/parallel
module load bbc2/samtools
module load bbc2/dupsifter

source activate biscuit

biscuit pileup \\
  -@ 30 \\
-o H3K27ac_EMseq.all_samples.vcf \\
  /home/rochelle.tiedemann/rothbart-secondary/Rochelle/alignment_indexes/hg38_index/biscuit_index/hg38.final.fa \\
  acVEH1_combined.bam acVEH2_combined.bam acVAL1_combined.bam acVAL2_combined.bam

bgzip -@ 30 H3K27ac_EMseq.all_samples.vcf

tabix -p vcf H3K27ac_EMseq.all_samples.vcf.gz

EOF
done
