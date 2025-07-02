#!/bin/bash
#
#SBATCH --partition=haswell  # haswell-p normal # partition (queue)
#SBATCH -N 1 # number of nodes
#SBATCH -J join
#SBATCH --mem 16G # memory pool for all cores
#SBATCH -t 0-02:59# time (D-HH:MM)
#SBATCH -o ./real/outputs/log.%j.out # STDOUT
#SBATCH -e ./real/outputs/log.%j.err # STDERR
#SBATCH --mail-type=END,FAIL # notifications for job done & fail
#SBATCH --mail-user=claudia.vasallo@upf.edu # send-to address



# Modules

module load R


# Config
CAT=$1
INPUTS=$2
OUTPUTS=$3


# Array Dependent Config
# If there is no array number assigned:
#   - that means that the script is being executed directly. 
#   - it takes the first item of the list as a parameter
# If there is an array number assigned:
#   - that means that the script is being executed through the cluster
#   - it takes the item correspondent to the array number
#if [ "${SLURM_ARRAY_TASK_ID}" == "" ]
#then
#  IDD=$(cat ${ARRAY_LIST} | sed -n 1p)
#else
#  IDD=$(cat ${ARRAY_LIST} | sed -n ${SLURM_ARRAY_TASK_ID}p)
#fi 


#Rscript ./real/scripts/join.chrs.cat.R ${CAT} ${IDD} ${INPUTS} ${OUTPUTS}


cat $(ls ${INPUTS}/*.${CAT}.genes | grep everyones )  > ${OUTPUTS}/everyones.${CAT}.genes
