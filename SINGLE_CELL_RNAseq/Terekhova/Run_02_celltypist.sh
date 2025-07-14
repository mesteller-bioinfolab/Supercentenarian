#!/bin/bash

#SBATCH --job-name=02_celltypist
#SBATCH --mail-type=END,FAIL
#SBATCH --cpus-per-task=1
#SBATCH --mem=150gb
#SBATCH --time=24:00:00
#SBATCH --output=./logs/%x_%j.log
####SBATCH --nodelist=c03
#SBATCH -p hpc

module load Anaconda3
source activate celltypist

python3 02_celltypist.py
