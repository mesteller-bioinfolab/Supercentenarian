#!/bin/bash
#
#SBATCH -p haswell # partition (queue)
#SBATCH -N 1 # number of nodes
#SBATCH -J mz
#SBATCH --mem 60G # memory pool for all cores
#SBATCH -t 0-03:00 # time (D-HH:MM)
#SBATCH -o ./real/outputs/log.%j.out # STDOUT
#SBATCH -e ./real/outputs/log.%j.err # STDERR
#SBATCH --mail-type=FAIL # notifications for job done & fail
#SBATCH --mail-user=claudia.vasallo@upf.edu # send-to address

# Modules
module load R


# Config
INPUTS=$1
OUTPUTS=$2
ARRAY_LIST=$3


# Array Dependent Config
# If there is no array number assigned:
#   - that means that the script is being executed directly. 
#   - it takes the first item of the list as a parameter
# If there is an array number assigned:
#   - that means that the script is being executed through the cluster
#   - it takes the item correspondent to the array number
if [ "${SLURM_ARRAY_TASK_ID}" == "" ]
then
  VERSION=$(cat ${ARRAY_LIST} | sed -n 1p)
else
  VERSION=$(cat ${ARRAY_LIST} | sed -n ${SLURM_ARRAY_TASK_ID}p)
fi 




Rscript ./real/scripts/mz.R ${INPUTS} ${OUTPUTS} ${VERSION}

