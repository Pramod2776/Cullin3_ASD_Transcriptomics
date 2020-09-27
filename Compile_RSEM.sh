#!/usr/bin/env bash
#SBATCH -t 01:00:00
#SBATCH	-p shared
#SBATCH --mem=3G
#SBATCH -J Comile_RSEM_Isoform
#SBATCH --job-name=Compile_RSEM
#SBATCH --output=Compile_RSEM.out

module load R
Rscript Compile_RSEM.R

