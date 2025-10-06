#!/bin/bash -l
#SBATCH -A sens2022005
#SBATCH -p node
#SBATCH -n 2
#SBATCH -t 60:00:00
#SBATCH -J PureCN_baseline
module load bioinfo-tools
module load GATK/4.3.0.0 
module load R/4.2.1
module load R_packages/4.2.1
export PURECN="/sw/apps/R_packages/4.2.1/bianca/PureCN/extdata"
export OUT_REF="/proj/sens2022005/nobackup/data/PREDIX_HER2_WES/PureCN/scratch/reference_files"
export OUT="/proj/sens2022005/nobackup/data/PREDIX_HER2_WES/PureCN/scratch/results/baseline"
for ((i=2; i<=20; i++)); #header(1)  # i=2; i<=194
do
Sample=$(awk -v i=$i 'NR==i{print $2}' /proj/sens2022005/nobackup/data/PREDIX_HER2_WES/PureCN/scratch/data/samplesheet_WESpair_baseline.txt)
mkdir $OUT/$Sample
Rscript $PURECN/PureCN.R --out $OUT/$Sample \
    --tumor /proj/sens2022005/nobackup/data/PREDIX_HER2_WES/PureCN/scratch/data/Baseline/$Sample".recal_coverage_loess.txt.gz" \
    --sampleid $Sample \
    --vcf /proj/sens2022005/nobackup/data/PREDIX_HER2_WES/PureCN/data/$Sample".vcf" \
    --fun-segmentation GATK4 \
    --normaldb $OUT_REF/normalDB_Twist_Comprehensive_Exome_hg38.rds \
    --intervals $OUT_REF/Twist_Comprehensive_Exome_baits_hg38_intervals.txt \
    --snp-blacklist $OUT_REF/hg38-blacklist.v2.bed \
    --max-segments 350 \
    --genome hg38 \
    --min-purity 0.05 --max-copy-number 8 \
    --rds=rda \
    --max-non-clonal 0.4 \
    --model betabin \
    --force --post-optimize --seed 123 \
    --cores 32
done
