#!/bin/bash
#
#SBATCH -p haswell # partition (queue)
#SBATCH -N 1 # number of nodes
#SBATCH -J compa
#SBATCH --mem 32G # memory pool for all cores
#SBATCH -t 0-03:59# time (D-HH:MM)
#SBATCH -o ./real/outputs4/log.%j.out # STDOUT
#SBATCH -e ./real/outputs4/log.%j.err # STDERR
#SBATCH --mail-type=END,FAIL # notifications for job done & fail
#SBATCH --mail-user=claudia.vasallo@upf.edu # send-to address



# Modules
module load R/3.6.0-foss-2018b

# Config

M116_OVER=$1
GENESET=$2
CAT_FILE=$3
OUTPUTS=$4


Rscript real/scripts/confidence.interval.R ${M116_OVER} ${GENESET} ${CAT_FILE} ${OUTPUTS} 
