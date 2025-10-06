#!/bin/bash -l
#SBATCH -A sens2022005
#SBATCH -p node
#SBATCH -n 4
#SBATCH -t 40:00:00
#SBATCH -J PON_vcf
module load bioinfo-tools
module load GATK/4.2.0.0
cd /proj/sens2022005/nobackup/data/PREDIX_HER2_WES/PureCN/scratch/PON
for ((i=190; i<=190; i++)); 
do 
Sample=$(awk -v i=$i 'NR==i{print $3}' /proj/sens2022005/nobackup/data/PREDIX_HER2_WES/PureCN/scratch/PON/PREDIX_HER2_blood_bam.txt)
bam=$(awk -v i=$i 'NR==i{print $2}' /proj/sens2022005/nobackup/data/PREDIX_HER2_WES/PureCN/scratch/PON/PREDIX_HER2_blood_bam.txt)
echo "parsing $Sample"
gatk Mutect2 --java-options "-Xmx16g -XX:+UseParallelGC -XX:ParallelGCThreads=64" --native-pair-hmm-threads 64 -R /sw/data/uppnex/ToolBox/hg38bundle/Homo_sapiens_assembly38.fasta -I $bam -max-mnp-distance 0 -O /proj/sens2022005/nobackup/data/PREDIX_HER2_WES/PureCN/scratch/PON/VCF/$Sample".vcf.gz"
done
