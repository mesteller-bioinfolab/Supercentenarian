#!/bin/bash
#
#SBATCH -p haswell # partition (queue)
#SBATCH -N 1 # number of nodes
#SBATCH -J genotypes
#SBATCH --mem 32G # memory pool for all cores
#SBATCH -t 9-23:59# time (D-HH:MM)
#SBATCH -o ./real/genotypes/log.%j.out # STDOUT
#SBATCH -e ./real/genotypes/log.%j.err # STDERR
#SBATCH --mail-type=END,FAIL # notifications for job done & fail
#SBATCH --mail-user=claudia.vasallo@upf.edu # send-to address



# Modules
module load R


# Config
ARRAY_LIST=$1
INPUTS_DIR=$2
OUTPUTS_DIR=$3
TISSUED=$4
IDD=$5


# Array Dependent Config
# If there is no array number assigned:
#   - that means that the script is being executed directly. 
#   - it takes the first item of the list as a parameter
# If there is an array number assigned:
#   - that means that the script is being executed through the cluster
#   - it takes the item correspondent to the array number
if [ "${SLURM_ARRAY_TASK_ID}" == "" ]
then
  CHR=$(cat ${ARRAY_LIST} | sed -n 1p)
else
  CHR=$(cat ${ARRAY_LIST} | sed -n ${SLURM_ARRAY_TASK_ID}p)
fi 




Rscript ./real/scripts/add_gt.R ${TISSUED} ${IDD} ${CHR} ${INPUTS_DIR} ${OUTPUTS_DIR}
