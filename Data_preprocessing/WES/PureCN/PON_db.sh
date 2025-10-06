#!/bin/bash -l
#SBATCH -A sens2022005
#SBATCH -p core
#SBATCH -n 8
#SBATCH -t 24:00:00
#SBATCH -J PON
module load bioinfo-tools
module load GATK/4.3.0.0 
cd /proj/sens2022005/nobackup/data/PREDIX_HER2_WES/PureCN/scratch/PON/VCF
gatk CreateSomaticPanelOfNormals -R /sw/data/uppnex/ToolBox/hg38bundle/Homo_sapiens_assembly38.fasta \
     --germline-resource /sw/data/GATK/Mutect2/af-only-gnomad.hg38.vcf.gz \
     --variant gendb://pon_db \
     --output /proj/sens2022005/nobackup/data/PREDIX_HER2_WES/PureCN/scratch/PON/pon_PREDIX_HER2.vcf.gz
