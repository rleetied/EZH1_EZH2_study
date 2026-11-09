module load bbc/bowtie2/bowtie2-2.4.1
for w in $1/*_L000_R1_001.fastq.gz; do 
#w=`echo $1`
output=`basename $w _L000_R1_001.fastq.gz`
r1=`echo $w`
r2=`echo $w|sed -e 's/R1/R2/g'`
if [ $output != "Undetermined" ] ; then
echo $output $r1 $r2
bowtie2 -p 20 -I 0 -X 700 --end-to-end --sensitive -x /secondary/projects/rothbart/Rochelle/hg38_fasta_SNAP/bowtie_index/hg38.SNAP -1 $r1 -2 $r2 -S $output.sam
#bowtie2 -p 20 -I 0 -X 700 --end-to-end --sensitive -x dros_bowtie -1 $r1 -2 $r2 -S $output.sam
awk -v MAQ=20 '$5>=MAQ && $2==99 || $5>=MAQ && $2==163 {print $3"\t"$4"\t"$4+$9-1}' $output.sam | awk '$2<=$3 {print $1"\t"$2"\t"$3"\t"$3-$2}' | sort -k1,1 -k2,2n > $output.bed
fi
done
