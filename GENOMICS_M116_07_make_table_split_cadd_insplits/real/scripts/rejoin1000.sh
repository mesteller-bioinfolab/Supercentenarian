#!/bin/bash
#
#SBATCH -p haswell # normal # partition (queue)
#SBATCH -N 1 # number of nodes
#SBATCH -J join1000
#SBATCH --mem 16G # memory pool for all cores
#SBATCH -t 0-00:19# time (D-HH:MM)
#SBATCH -o ./real/tmp/log.%j.out # STDOUT
#SBATCH -e ./real/tmp/log.%j.err # STDERR
#SBATCH --mail-type=END,FAIL # notifications for job done & fail
#SBATCH --mail-user=claudia.vasallo@upf.edu # send-to address


ARRAY_LIST=$1
INPUTS_DIR=$2
OUT_DIR=$3
TISSUES=$4
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




CHR_FILES=(${INPUTS_DIR}/${TISSUES}/${IDD}/${TISSUES}.${IDD}.${CHR}_*.table.split2.insplits)


HEADER=./real/inputs/splitheader


cp ${HEADER} ${OUT_DIR}/${TISSUES}/${IDD}/${TISSUES}.${IDD}.${CHR}.table.split2
cat ${CHR_FILES[@]} | grep  -v -w TRANSCRIPTION_FACTORS >> ${OUT_DIR}/${TISSUES}/${IDD}/${TISSUES}.${IDD}.${CHR}.table.split2




