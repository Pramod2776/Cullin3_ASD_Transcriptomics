#!/usr/bin/env bash

#SBATCH -t 10:00:00
#SBATCH -p shared
#SBATCH --mem=120G
#SBATCH -J build_index
#SBATCH --output=build_index.out

#source activate rnaseq
STAR --runThreadN 16 --runMode genomeGenerate \
    --genomeDir references/Species \
    --genomeFastaFiles Species_Genome/version.genome.fa \
    --sjdbGTFfile  Species_Genome/gencode.version.annotation.gtf \
    --sjdbOverhang 100


rsem-prepare-reference \
    --gtf Species_Genome/gencode.version.annotation.gtf \
    Species_Genome/version.genome.fa \
    references/RSEM_index/RSEM_Species
