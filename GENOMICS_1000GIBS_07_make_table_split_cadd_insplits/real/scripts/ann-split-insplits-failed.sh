#!/bin/bash
#
#SBATCH --partition=normal    #haswell  #normal # partition (queue)
#SBATCH -N 1 # number of nodes
#SBATCH -J faster
#SBATCH --mem 32G # memory pool for all cores
#SBATCH -t 0-23:59# time (D-HH:MM)
#SBATCH -o ./real/outputs/log.%j.out # STDOUT
#SBATCH -e ./real/outputs/log.%j.err # STDERR
#SBATCH --mail-type=END,FAIL # notifications for job done & fail
#SBATCH --mail-user=claudia.vasallo@upf.edu # send-to address



# Modules
module load Python/3.7.2-GCCcore-8.2.0

# Config
ARRAY_LIST=$1
INPUTS_DIR=$2
OUTPUTS_DIR=$3
KEY=$4

# Array Dependent Config
# If there is no array number assigned:
#   - that means that the script is being executed directly. 
#   - it takes the first item of the list as a parameter
# If there is an array number assigned:
#   - that means that the script is being executed through the cluster
#   - it takes the item correspondent to the array number
if [ "${SLURM_ARRAY_TASK_ID}" == "" ]
then
  CHR=$(cat ${ARRAY_LIST} | sed -n 1p | cut -d '.' -f4)
  TISSUED=$(cat ${ARRAY_LIST} | sed -n 1p   |  cut -d'/' -f4)
  IDD=$(cat ${ARRAY_LIST} | sed -n 1p | cut -d '.' -f3) 
else
  CHR=$(cat ${ARRAY_LIST} | sed -n ${SLURM_ARRAY_TASK_ID}p | cut -d '.' -f4)
  TISSUED=$(cat ${ARRAY_LIST}  |  sed -n ${SLURM_ARRAY_TASK_ID}p |  cut -d'/' -f4)
  IDD=$(cat ${ARRAY_LIST} | sed -n ${SLURM_ARRAY_TASK_ID}p | cut -d '.' -f3)
fi 


echo ${TISSUED}
echo ${IDD}
echo ${CHR}


# Run
# filtered.table before
TABLE=${INPUTS_DIR}/${TISSUED}/${IDD}/${TISSUED}.${IDD}.${CHR}.header

python3 ./real/scripts/split-VEP-field-CADD_2.py -i ${TABLE} -o ${OUTPUTS_DIR}/${TISSUED}/${IDD}/${TISSUED}.${IDD}.${CHR}.table.split2.insplits -k ${KEY}

