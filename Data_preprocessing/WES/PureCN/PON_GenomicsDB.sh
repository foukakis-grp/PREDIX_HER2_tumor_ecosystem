#!/bin/bash -l
#SBATCH -A sens2022005
#SBATCH -p node
#SBATCH -n 1
#SBATCH -t 40:00:00
#SBATCH -J PON
module load bioinfo-tools
module load GATK/4.2.0.0 
cd /proj/sens2022005/nobackup/data/PREDIX_HER2_WES/PureCN/scratch/PON/VCF
gatk GenomicsDBImport --java-options "-Xmx8g -XX:+UseParallelGC -XX:ParallelGCThreads=16" \
       --arguments_file /proj/sens2022005/nobackup/data/PREDIX_HER2_WES/PureCN/scratch/PON/VCF/PON_sample \
       --merge-input-intervals \
       -R /sw/data/uppnex/ToolBox/hg38bundle/Homo_sapiens_assembly38.fasta \
       -L /home/kangwang/nf-core/nf-core-sarek-2.7.1/Sarek-data/Twist_Comprehensive_Exome_Covered_Targets_GRCh38.pad.sort.merge.bed \
       --genomicsdb-workspace-path pon_db \
       --batch-size 50 \
       --reader-threads 16 
