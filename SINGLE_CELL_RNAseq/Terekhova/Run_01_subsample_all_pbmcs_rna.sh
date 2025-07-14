#!/bin/bash

#SBATCH --job-name=01_subsample_all_pbmcs_rna
#SBATCH --mail-type=END,FAIL
#SBATCH --cpus-per-task=1
#SBATCH --mem=90gb
#SBATCH --time=1:00:00
#SBATCH --output=./logs/%x_%j.log
#SBATCH -p hpc

module load Anaconda3
source activate celltypist

python3 01_subsample_all_pbmcs_rna.py
