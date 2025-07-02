#!/bin/bash
#
#SBATCH --partition=haswell       # normal  # haswell-p normal # partition (queue)
#SBATCH -N 1 # number of nodes
#SBATCH -J tab-summary
#SBATCH --mem 16G # memory pool for all cores
#SBATCH -t 0-00:29# time (D-HH:MM)
#SBATCH -o ./real/outputs/log.%j.out # STDOUT
#SBATCH -e ./real/outputs/log.%j.err # STDERR
#SBATCH --mail-type=END,FAIL # notifications for job done & fail
#SBATCH --mail-user=claudia.vasallo@upf.edu # send-to address


# Modules

#module load R
#module load R/3.6.0-foss-2018b



INPUTS=$1
OUTPUTS=$2
TISSUED=$3
ARRAY_LIST=$4
CAT_FILE=$5


# Array Dependent Config
# If there is no array number assigned:
#   - that means that the script is being executed directly. 
#   - it takes the first item of the list as a parameter
# If there is an array number assigned:
#   - that means that the script is being executed through the cluster
#   - it takes the item correspondent to the array number
if [ "${SLURM_ARRAY_TASK_ID}" == "" ]
then
  IDD=$(cat ${ARRAY_LIST} | sed -n 1p)
else
  IDD=$(cat ${ARRAY_LIST} | sed -n ${SLURM_ARRAY_TASK_ID}p)
fi 










# extra too cats

CATEGORY=non_coding
cat ${INPUTS}/${TISSUED}/${IDD}/*.eur_rare_non_coding.genes > ${OUTPUTS}/${TISSUED}.${IDD}.eur_rare_${CATEGORY}.genes

CATEGORY=coding_with_nc

cat ${INPUTS}/${TISSUED}/${IDD}/*.eur_rare_coding_with_nc.genes > ${OUTPUTS}/${TISSUED}.${IDD}.eur_rare_${CATEGORY}.genes






