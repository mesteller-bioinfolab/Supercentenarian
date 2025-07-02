#!/bin/bash
#
#SBATCH -p haswell # partition (queue)
#SBATCH -N 1 # number of nodes
#SBATCH -J resamp
#SBATCH --mem 16G # memory pool for all cores
#SBATCH -t 0-03:59# time (D-HH:MM)
#SBATCH -o ./real/outputs/log.%j.out # STDOUT
#SBATCH -e ./real/outputs/log.%j.err # STDERR

# Modules
module load R/3.6.0-foss-2018b

# Config

INPUTS_DIR=$1
ARRAY_LIST=$2
SET_FILE=$3


# Array Dependent Config
# If there is no array number assigned:
#   - that means that the script is being executed directly. 
#   - it takes the first item of the list as a parameter
# If there is an array number assigned:
#   - that means that the script is being executed through the cluster
#   - it takes the item correspondent to the array number
if [ "${SLURM_ARRAY_TASK_ID}" == "" ]
then
  CAT=$(cat ${ARRAY_LIST} | sed -n 1p)
else
  CAT=$(cat ${ARRAY_LIST} | sed -n ${SLURM_ARRAY_TASK_ID}p)
fi 

SETS=$(cat ${SET_FILE})


for SET in ${SETS}
do
Rscript ./real/scripts/resampling_res_cat.R ${INPUTS_DIR} ${CAT} ${SET}
done
