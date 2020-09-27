#!/usr/bin/env bash

#SBATCH -t 05:00:00
#SBATCH	-p shared
#SBATCH --mem=100G
#SBATCH -J Cul3_CX
#SBATCH --job-name=Cul3_CX
#SBATCH --output=Cul3_CX_%A_%a.out
#SBATCH --array=1-2

source activate RNASeq
export LANG=en_US.UTF-8
export LC_ALL=en_US.UTF-8

FILE1=`cat Path/Project/input_filenames/files.txt | awk -v SLURM_ARRAY_TASK_ID=$((${SLURM_ARRAY_TASK_ID} * 2)) 'NR==SLURM_ARRAY_TASK_ID {print}'`
FILE2=`cat PathProject/input_filenames/files.txt | awk -v SLURM_ARRAY_TASK_ID=$(((${SLURM_ARRAY_TASK_ID} * 2)-1)) 'NR==SLURM_ARRAY_TASK_ID {print}'`

FILE_PATH="PATH/Project/folder/"
TEMP_DIR=${SLURM_ARRAY_TASK_ID}
mkdir -p -- "PATH/folder/$TEMP_DIR"
chunky run PEC_DAC_RNAseq_16p_isoform.py \
    --reads  "$FILE_PATH$FILE2":"$FILE_PATH$FILE1" \
    --output PATH/folder/$TEMP_DIR \
    --forward-adapter XXX --reverse-adapter XXX

