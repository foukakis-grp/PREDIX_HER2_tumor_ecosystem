#!/bin/bash
#SBATCH -A sens2022005
#SBATCH -p node
#SBATCH -n 1
#SBATCH -t 10:00:00
#SBATCH -J Bcell
 
module load bioinfo-tools R/4.3.1 R_packages/4.3.1 Seurat/5.2.1
 
cd /proj/sens2022005/Xenium/PREDIX_HER2/code_final_sep25/subcluster

Rscript Bcell.R
