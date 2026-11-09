#!/bin/bash
set -euxo pipefail

#This script is intended to generate a job script that can be submitted to a PBS manager
#To generate the job scripts
#1) Modify any of the paths/function calls/etc below to meet the analysis needs/tool used
#2) bash your_generator.sh
#Submit to the queue
#Note if the number of jobs is large, you will receive an e-mail for start/stop/fail for each job
#3) Example submission loop
#for i in {0..9}; do; qsub pbs.starscript.$i; sleep 1; done

#Build up a space sep'd array of sample/file prefix
declare -a SAMPLE=(acVAL1_combined acVAL2_combined acVEH1_combined acVEH2_combined me1TAZ1 me1TAZ2 me1VAL2_combined me1VEH1 me1VEH2_combined meVAL1)

#set a tmpdir for parallel
tmpdir='/home/rochelle.tiedemann/rothbart-secondary/Rochelle/PR001883_ROTS_EMChIP/trimmed_fastq'

#Adjust the length of the array to iterate over to match the number of samples/files to analyze
for i in {0..9}; do
    cat > slurm.biscuitscript.$i <<EOF
#!/bin/bash
#
#SBATCH --job-name=biscuit_${SAMPLE[i]}
#SBATCH --ntasks=1
#SBATCH --time=168:00:00

#Change to WGBS directory
cd /home/rochelle.tiedemann/rothbart-secondary/Rochelle/PR001883_ROTS_EMChIP/trimmed_fastq

#load modules
module load bbc2/samblaster
module load bbc2/pigz
module load bbc2/parallel

#these are directional libs
#Launch biscuit
/home/rochelle.tiedemann/rothbart-secondary/Rochelle/biscuit/biscuit align -@ 30 -R '@RG\tLB:hg38\tID:WGBS_${SAMPLE[i]}\tPL:Illumina\tPU:novoseq6000\tSM:${SAMPLE[i]}' \
/home/rochelle.tiedemann/rothbart-secondary/Rochelle/alignment_indexes/hg38_index/biscuit_index/hg38.final.fa \
${SAMPLE[i]}_L000_R1_001_val_1.fq.gz ${SAMPLE[i]}_L000_R2_001_val_2.fq.gz | \
samblaster --addMateTags | parallel --tmpdir ${tmpdir} --pipe --tee {} ::: 'samblaster -a -e -u ${SAMPLE[i]}.clipped.fastq -d ${SAMPLE[i]}.disc.hg38.sam -s ${SAMPLE[i]}.split.hg38.sam -o /dev/null' \
'samtools view -hb | samtools sort -@ 8 -m 5G -o ${SAMPLE[i]}.sorted.markdup.withdisc_split_clip.hg38.bam -O BAM -'

#you can collect the split and discordant bams out of smoove too
#sort and convert to bams
samtools sort -@ 20 -o ${SAMPLE[i]}.disc.hg38.bam -O BAM ${SAMPLE[i]}.disc.hg38.sam

#index
samtools index ${SAMPLE[i]}.disc.hg38.bam

#sort and convert to bams
samtools sort -o ${SAMPLE[i]}.split.hg38.bam -O BAM ${SAMPLE[i]}.split.hg38.sam

#index
samtools index ${SAMPLE[i]}.split.hg38.bam

#compress
pigz -p 40 ${SAMPLE[i]}.clipped.fastq

#index the full bam
samtools index ${SAMPLE[i]}.sorted.markdup.withdisc_split_clip.hg38.bam

#clean up the sams...
if [[ -f ${SAMPLE[i]}.disc.hg38.bam ]] && [[ -f ${SAMPLE[i]}.split.hg38.bam ]]; then rm ${SAMPLE[i]}*.sam; fi

EOF
done
EOF
