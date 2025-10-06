#!/bin/bash -l
#SBATCH -A sens2022005
#SBATCH -p core
#SBATCH -n 1
#SBATCH -t 100:00:00
#SBATCH -J lohhla
module purge
module load bioinfo-tools
module load LOHHLA/20210129-00744c5
for ((i=2; i<=194; i++)); #header(1)
do
sampleID=$(awk -v i=$i 'NR==i{print $1}' /proj/sens2022005/nobackup/lohhla/PREDIX_HER2_baseline/data/lohhla_baseline.txt)
hlatyping=$(awk -v i=$i 'NR==i{print $2}' /proj/sens2022005/nobackup/lohhla/PREDIX_HER2_baseline/data/lohhla_baseline.txt)
patientID=$(awk -v i=$i 'NR==i{print $3}' /proj/sens2022005/nobackup/lohhla/PREDIX_HER2_baseline/data/lohhla_baseline.txt)
tumorID=$(awk -v i=$i 'NR==i{print $12}' /proj/sens2022005/nobackup/lohhla/PREDIX_HER2_baseline/data/lohhla_baseline.txt)
normalID=$(awk -v i=$i 'NR==i{print $14}' /proj/sens2022005/nobackup/lohhla/PREDIX_HER2_baseline/data/lohhla_baseline.txt)
#bam_T=$(awk -v i=$i 'NR==i{print $13}' /proj/sens2022005/nobackup/lohhla/PREDIX_HER2_baseline/data/lohhla_baseline.txt)
#bam_N=$(awk -v i=$i 'NR==i{print $15}' /proj/sens2022005/nobackup/lohhla/PREDIX_HER2_baseline/data/lohhla_baseline.txt)
bam_T=/proj/sens2022005/nobackup/PREDIX_HER2_md_bam/$tumorID/DuplicatesMarked
bam_N=/proj/sens2022005/nobackup/PREDIX_HER2_md_bam/$normalID/DuplicatesMarked
cd /proj/sens2022005/nobackup/lohhla/PREDIX_HER2_baseline/results
mkdir $sampleID
BAMDir=/proj/sens2022005/nobackup/lohhla/PREDIX_HER2_baseline/results/$sampleID
cd $BAMDir
mkdir lohhla
cp $bam_T/$tumorID".bam" .
cp $bam_T/$tumorID".bam.bai" .
cp $bam_N/$normalID".bam" .
cp $bam_N/$normalID".bam.bai" .
outputDir=/proj/sens2022005/nobackup/lohhla/PREDIX_HER2_baseline/results/$sampleID/lohhla
Rscript /proj/sens2022005/nobackup/lohhla/LOHHLAscript_uppmax_GRCh38.R --patientId $sampleID \
--outputDir $outputDir \
--normalBAMfile /proj/sens2022005/nobackup/lohhla/PREDIX_HER2_baseline/results/$sampleID/$normalID".bam" \
--BAMDir /proj/sens2022005/nobackup/lohhla/PREDIX_HER2_baseline/results/$sampleID \
--hlaPath /proj/sens2022005/nobackup/lohhla/PREDIX_HER2_baseline/data/HLAtyping/$hlatyping/hlas \
--HLAfastaLoc /proj/sens2022005/nobackup/lohhla/PREDIX_HER2_baseline/data/abc_complete.fasta \
--CopyNumLoc /proj/sens2022005/nobackup/lohhla/PREDIX_HER2_baseline/data/purity_ploidy_solution/$sampleID.solution.txt \
--mappingStep TRUE --minCoverageFilter 10 --fishingStep TRUE --cleanUp FALSE
rm -r $tumorID".bam"*
rm -r $normalID".bam"*
done
