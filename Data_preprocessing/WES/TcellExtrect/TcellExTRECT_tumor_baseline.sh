#!/bin/bash -l
#SBATCH -A sens2022005
#SBATCH -p node
#SBATCH -n 4
#SBATCH -t 24:00:00
#SBATCH -J TcellExTRECT_tumor
module load bioinfo-tools R_packages/4.1.1 samtools/1.17
R --no-save --quiet < /proj/sens2022005/nobackup/TcellExTRECT/TcellExTRECT_tumor_baseline.R
