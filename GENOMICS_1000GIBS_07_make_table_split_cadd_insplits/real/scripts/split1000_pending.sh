#!/bin/bash
#
#SBATCH -p haswell # normal # partition (queue)
#SBATCH -N 1 # number of nodes
#SBATCH -J split1000
#SBATCH --mem 64G # memory pool for all cores
#SBATCH -t 0-07:59# time (D-HH:MM)
#SBATCH -o ./real/tmp/log.%j.out # STDOUT
#SBATCH -e ./real/tmp/log.%j.err # STDERR
#SBATCH --mail-type=END,FAIL # notifications for job done & fail
#SBATCH --mail-user=claudia.vasallo@upf.edu # send-to address


ARRAY_LIST=$1
INPUTS_DIR=$2
TMP_DIR=$3
TISSUES=$4



# Array Dependent Config
# If there is no array number assigned:
#   - that means that the script is being executed directly. 
#   - it takes the first item of the list as a parameter
# If there is an array number assigned:
#   - that means that the script is being executed through the cluster
#   - it takes the item correspondent to the array number
if [ "${SLURM_ARRAY_TASK_ID}" == "" ]
then
  CHR=$(cat ${ARRAY_LIST} | sed -n 1p | cut -f1)
  IDD=$(cat ${ARRAY_LIST} | sed -n 1p | cut -f2)
else
  CHR=$(cat ${ARRAY_LIST} | sed -n ${SLURM_ARRAY_TASK_ID}p | cut -f1)
  IDD=$(cat ${ARRAY_LIST} | sed -n ${SLURM_ARRAY_TASK_ID}p | cut -f2)
fi 


mkdir ${TMP_DIR}/${TISSUES}
mkdir ${TMP_DIR}/${TISSUES}/${IDD}

CHR_FILE=${INPUTS_DIR}/${TISSUES}/${IDD}/${TISSUES}.${IDD}.${CHR}.table.split

 
HEADER=./real/inputs/header

# make splits
cat ${CHR_FILE} | tail -n+2 > ${TMP_DIR}/${TISSUES}/${IDD}/${TISSUES}.${IDD}.HEADERLESS_${CHR}
#split --verbose -l1000 ${TMP_DIR}/${TISSUES}/${IDD}/${TISSUES}.${IDD}.HEADERLESS_${CHR} ${TMP_DIR}/${TISSUES}/${IDD}/${TISSUES}.${IDD}.${CHR}_
split -l1000 ${TMP_DIR}/${TISSUES}/${IDD}/${TISSUES}.${IDD}.HEADERLESS_${CHR} ${TMP_DIR}/${TISSUES}/${IDD}/${TISSUES}.${IDD}.${CHR}_
LCHRS=($(ls ${TMP_DIR}/${TISSUES}/${IDD}/${TISSUES}.${IDD}.${CHR}_* | grep -v header | grep -v insplits ))

echo ${#LCHRS[@]}


#for NSPLIT in `seq 0 $(( ${#LCHRS[@]} - 1 ))`
for SPLIT in ${LCHRS[@]}
do
echo ${SPLIT}
#ECHR=$(echo ${LCHRS[${NSPLIT}]} | cut -d '_' -f2)
ECHR=$(echo ${SPLIT} | cut -d '_' -f2)
echo ${ECHR}
cat ${HEADER} > ${TMP_DIR}/${TISSUES}/${IDD}/${TISSUES}.${IDD}.${CHR}_${ECHR}.header
echo ${TMP_DIR}/${TISSUES}/${IDD}/${TISSUES}.${IDD}.${CHR}_${ECHR}.header
#cat ${LCHRS[${NSPLIT}]} >> ${TMP_DIR}/${TISSUES}/${IDD}/${TISSUES}.${IDD}.${CHR}_${ECHR}.header
cat ${SPLIT} >> ${TMP_DIR}/${TISSUES}/${IDD}/${TISSUES}.${IDD}.${CHR}_${ECHR}.header
done



ls ${TMP_DIR}/${TISSUES}/${IDD}/${TISSUES}.${IDD}.*.header > ./real/inputs/${TISSUES}.${IDD}.chrs_split.list



