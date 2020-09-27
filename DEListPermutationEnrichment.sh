#!/usr/bin/env bash

#SBATCH -p shared
#SBATCH -t 10:00:00
#SBATCH -J DEListPermutation_Redo
#SBATCH -o logs/DEListPermutation_Redo_%A_%a.out
#SBATCH --mem=16G
#SBATCH --cpus-per-task=16
#SBATCH --array=2-100

module load R
Rscript DEListPermutationEnrichment_Redo.R
