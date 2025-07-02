#!/bin/bash
#
#SBATCH -p haswell # partition (queue)
#SBATCH -N 1 # number of nodes
#SBATCH -J fet
#SBATCH --mem 16G # memory pool for all cores
#SBATCH -t 0-00:29# time (D-HH:MM)
#SBATCH -o ./real/outputs3/log.%j.out # STDOUT
#SBATCH -e ./real/outputs3/log.%j.err # STDERR
#SBATCH --mail-type=END,FAIL # notifications for job done & fail
#SBATCH --mail-user=claudia.vasallo@upf.edu # send-to address


module load R

INPUTS=$1
OUTPUTS=$2
TISSUED=$3
ARRAY_LIST=$4
GENELISTS_FOLDER=$5
CATEGORY=$6
#CATEGORY_BG=$7
BG_TYPE=$7


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
   echo ${IDD}
else
  IDD=$(cat ${ARRAY_LIST} | sed -n ${SLURM_ARRAY_TASK_ID}p)
fi 



module load R



GENELISTS=( $(ls ${GENELISTS_FOLDER} | grep txt |  cut -d '.' -f1  | rev | cut -d'_' -f2,3,4 | rev  ) )




cat ${OUTPUTS}/fet.table.header > ${OUTPUTS}/${TISSUED}.${IDD}.${CATEGORY}.fet



for GENELIST in ${GENELISTS[@]}

do

echo ${GENELIST}

# 2TISSUES.M116.hagr_cellsignatures.overlaps.table
OVER_FILE=${INPUTS}/${TISSUED}.${IDD}.${GENELIST}.overlaps.table

FET=$(Rscript ./real/scripts/fisher.test.R ${TISSUED} ${IDD} ${GENELIST} ${OVER_FILE} ${CATEGORY} ${BG_TYPE})

#echo ${IDD} ${GENELIST} ${CATEGORY} ${FET}   >> ${OUTPUTS}/${TISSUED}.${IDD}.fet

TOT_BG=$(echo ${FET} | cut -d '"' -f2 | cut -d '"' -f1)

# TOT_BG now replaces all, full row taken 


#IN_BG=$(echo ${FET} | cut -d ' ' -f3)
#TOT_CAT=$(echo ${FET} | cut -d ' ' -f4)
#IN_CAT=$(echo ${FET} | cut -d ' ' -f5)

#PVAL=$(echo ${FET} | cut -d ' ' -f6)
#OR=$(echo ${FET} | cut -d ' ' -f7)
#CI_L=$(echo ${FET} | cut -d ' ' -f8)
#CI_U=$(echo ${FET} | cut -d ' ' -f9 )


# mes <- c( pop, m, k, q, f_test$p.value, f_test$estimate, f_test$conf.int[1:2])

echo ${TISSUED} ${IDD} ${GENELIST} ${CATEGORY} ${TOT_BG} ${IN_BG} ${TOT_CAT} ${IN_CAT} ${PVAL} ${OR} ${CI_L} ${CI_U} >> ${OUTPUTS}/${TISSUED}.${IDD}.${CATEGORY}.fet


done
