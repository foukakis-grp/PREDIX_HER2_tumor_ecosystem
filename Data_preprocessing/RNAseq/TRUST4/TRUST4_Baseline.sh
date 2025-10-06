#!/bin/bash -l
#SBATCH -A sens2022005
#SBATCH -p node
#SBATCH -n 4
#SBATCH -t 150:00:00
#SBATCH -J TRUST4_Baseline
module load bioinfo-tools
module load TRUST4/1.0.8
module load R/4.1.1
module load R_packages/4.1.1
cd /home/kangwang/TRUST4/result
for ((i=2; i<=200; i++)); #header(1)
do
Sample=$(awk -v i=$i 'NR==i{print $1}' /home/kangwang/TRUST4/data/samplesheet_rna_Baseline.tsv)
fastq1=$(awk -v i=$i 'NR==i{print $2}' /home/kangwang/TRUST4/data/samplesheet_rna_Baseline.tsv)
fastq2=$(awk -v i=$i 'NR==i{print $3}' /home/kangwang/TRUST4/data/samplesheet_rna_Baseline.tsv)
cd /home/kangwang/TRUST4/result/Baseline
mkdir $Sample
cd /home/kangwang/TRUST4/result/Baseline/$Sample
run-trust4 -t 64 \
-f /home/kangwang/TRUST4/reference/TRUST4_ref/hg38_bcrtcr.fa \
--ref /home/kangwang/TRUST4/reference/TRUST4_ref/human_IMGT+C.fa \
-1 $fastq1 \
-2 $fastq2 \
-o $Sample
sed '1,1s/#count/count/g' ${Sample}_report.tsv > ${Sample}_report.processed.tsv
Rscript /home/kangwang/TRUST4/src/trust4_bcr_process.R --cdr3 ${Sample}_report.processed.tsv \
--sampleid $Sample \
--stat /home/kangwang/TRUST4/data/RNA_markdup.sorted.bam.stats/${Sample}.markdup.sorted.bam.stats.txt \
--outdir /home/kangwang/TRUST4/result/Baseline/$Sample/$Sample
Rscript /home/kangwang/TRUST4/src/trust4_tcr_process.R --cdr3 ${Sample}_report.processed.tsv \
--sampleid $Sample \
--stat /home/kangwang/TRUST4/data/RNA_markdup.sorted.bam.stats/${Sample}.markdup.sorted.bam.stats.txt \
--outdir /home/kangwang/TRUST4/result/Baseline/$Sample/$Sample
done
