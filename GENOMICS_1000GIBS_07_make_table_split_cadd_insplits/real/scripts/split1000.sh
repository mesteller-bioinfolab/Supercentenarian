#!/bin/bash
#
#SBATCH -p haswell # partition (queue)
#SBATCH -N 1 # number of nodes
#SBATCH -J split1000
#SBATCH --mem 32G # memory pool for all cores
#SBATCH -t 0-07:59# time (D-HH:MM)
#SBATCH -o ./real/tmp2/log.%j.out # STDOUT
#SBATCH -e ./real/tmp2/log.%j.err # STDERR
#SBATCH --mail-type=END,FAIL # notifications for job done & fail
#SBATCH --mail-user=claudia.vasallo@upf.edu # send-to address


ARRAY_LIST=$1
INPUTS_DIR=$2
TMPS_DIR=$3
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
   CHR=${ARRAY_LIST}
else
  CHR=$(cat ${ARRAY_LIST} | sed -n ${SLURM_ARRAY_TASK_ID}p)
fi 



CHR_FILE=${INPUTS_DIR}/${TISSUES}/${IDD}/${TISSUES}.${IDD}.${CHR}.table.split

 
HEADER=./real/inputs/header

# make splits
cat ${CHR_FILE} | tail -n+2 > ${TMPS_DIR}/${TISSUES}/${IDD}/${TISSUES}.${IDD}.HEADERLESS_${CHR}
split --verbose -l1000 ${TMPS_DIR}/${TISSUES}/${IDD}/${TISSUES}.${IDD}.HEADERLESS_${CHR} ${TMPS_DIR}/${TISSUES}/${IDD}/${TISSUES}.${IDD}.${CHR}_
LCHRS=($(ls ${TMPS_DIR}/${TISSUES}/${IDD}/${TISSUES}.${IDD}.${CHR}_* | grep -v header | grep -v insplits ))
echo ${#LCHRS[@]}


for NSPLIT in `seq 0 $(( ${#LCHRS[@]} - 1 ))`
do
echo ${NSPLIT}
ECHR=$(echo ${LCHRS[${NSPLIT}]} | cut -d '_' -f2)
echo ${ECHR}
cp ${HEADER} ${TMPS_DIR}/${TISSUES}/${IDD}/${TISSUES}.${IDD}.${CHR}_${ECHR}.header
echo ${TMPS_DIR}/${TISSUES}/${IDD}/${TISSUES}.${IDD}.${CHR}_${ECHR}.header
cat ${LCHRS[${NSPLIT}]} >> ${TMPS_DIR}/${TISSUES}/${IDD}/${TISSUES}.${IDD}.${CHR}_${ECHR}.header
done



#ls ${TMP_DIR}/${TISSUES}/${IDD}/${TISSUES}.${IDD}.*.header > ./real/inputs/${TISSUES}.${IDD}.chrs_split.list

ls ${TMP_DIR}/${TISSUES}/${IDD}/${TISSUES}.${IDD}.${CHR}*.header > ./${TISSUES}.${IDD}.${CHR}.chrs_split.list


