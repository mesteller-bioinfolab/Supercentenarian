#!/bin/sh
#SBATCH -p haswell # partition (queue)
#SBATCH -N 1 # number of nodes
#SBATCH -J exclusi
#SBATCH --mem 32G # memory pool for all cores
#SBATCH -t 0-01:00 # time (D-HH:MM)
#SBATCH -o ./log.%j.out # STDOUT
#SBATCH -e ./log.%j.err # STDERR
#SBATCH --mail-type=FAIL # notifications for job done & fail
#SBATCH --mail-user=claudia.vasallo@upf.edu # send-to address



#CATEGORY=$1

ARRAY_LIST=$1
IBS_DIR=$2
LEVEL=$3

# Array Dependent Config
# If there is no array number assigned:
#   - that means that the script is being executed directly. 
#   - it takes the first item of the list as a parameter
# If there is an array number assigned:
#   - that means that the script is being executed through the cluster
#   - it takes the item correspondent to the array number
if [ "${SLURM_ARRAY_TASK_ID}" == "" ]
then
  CATEGORY=$(cat ${ARRAY_LIST} | sed -n 1p )
else
  CATEGORY=$(cat ${ARRAY_LIST} | sed -n ${SLURM_ARRAY_TASK_ID}p)
fi 



cat $(ls ./real/outputs/*.${LEVEL}.txt | grep "/${CATEGORY}_" )  | grep -v -f ${IBS_DIR}/1000GIBS.${CATEGORY}.${LEVEL} | cut -f1 | sort -u  > M116.${CATEGORY}.${LEVEL}.excl
