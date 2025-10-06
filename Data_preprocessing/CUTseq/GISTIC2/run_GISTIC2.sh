#!/bin/bash -l
#SBATCH -A sens2022005
#SBATCH -p core
#SBATCH -n 2
#SBATCH -t 10:00:00
#SBATCH -J GISTIC
module load bioinfo-tools GISTIC
cd /proj/sens2022005/nobackup/GISTIC2
basedir="/proj/sens2022005/nobackup/GISTIC2/results/PREDIX_HER2_baseline"
segfile="/proj/sens2022005/nobackup/GISTIC2/data/PREDIX_100kb_baseline_logratio.txt"
refgenefile="/sw/bioinfo/GISTIC/2.0.23/rackham/refgenefiles/hg38.UCSC.add_miR.160920.refgene.mat"
gistic2 -b $basedir -seg $segfile -refgene $refgenefile -genegistic 1 -ta 0.1 -td -0.1 -smallmem 1 -broad 1 -brlen 0.98 -rx 1 -conf 0.99 -js 4 -arb 1 -maxseg 2000 -armpeel 1 -savegene 1 -gcm extreme
