#!/bin/bash -l
 
#SBATCH -A sens2022005
#SBATCH -p core
#SBATCH -n 16
#SBATCH -t 1-12:00:00
#SBATCH -J CUTseq-1e5

module load bioinfo-tools python3/3.7.2 snakemake/6.9.1 GATK/4.1.0.0 R_packages/3.4.0

export LD_PRELOAD=/sw/libs/openblas/0.3.15/bianca/lib/libopenblas.so

snakemake --cores 16 --rerun-incomplete -s /proj/sens2022005/nobackup/CUTseq/snakemake/gatk-cnv/1e5_merged/Snakefile
