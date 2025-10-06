#!/bin/bash -l
#SBATCH -A sens2022005
#SBATCH -p node
#SBATCH -n 6
#SBATCH -t 12:00:00
#SBATCH -J coverage_baseline
module load bioinfo-tools
module load R/4.2.1
module load R_packages/4.2.1
cd /proj/nobackup/sens2022005/data/PREDIX_HER2_WES/PureCN/scratch/data/Baseline
export PURECN="/sw/apps/R_packages/4.2.1/bianca/PureCN/extdata"
export OUT_REF="/proj/sens2022005/nobackup/data/PREDIX_HER2_WES/PureCN/scratch/reference_files"
export OUT="/proj/nobackup/sens2022005/data/PREDIX_HER2_WES/PureCN/scratch/data/Baseline"
Rscript $PURECN/Coverage.R --out-dir $OUT \
    --bam /proj/nobackup/sens2022005/data/PREDIX_HER2_WES/PureCN/scratch/data/baseline_bam.list \
    --intervals $OUT_REF/Twist_Comprehensive_Exome_baits_hg38_intervals.txt \
    --cores 96 --parallel --force
