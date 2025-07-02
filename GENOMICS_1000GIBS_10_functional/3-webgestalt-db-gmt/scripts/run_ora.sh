#!/bin/bash
#
#SBATCH -p haswell # partition (queue)
#SBATCH -N 1 # number of nodes
#SBATCH -J webgesalt
#SBATCH --mem 32G # memory pool for all cores
#SBATCH -t 0-01:00 # time (D-HH:MM)
#SBATCH -o ./real/outputs/log.%j.out # STDOUT
#SBATCH -e ./real/outputs/log.%j.err # STDERR
#SBATCH --mail-type=FAIL # notifications for job done & fail
#SBATCH --mail-user=claudia.vasallo@upf.edu # send-to address

# Modules
module load R


# Config
INPUTS=$1
GMT=$2
OUTPUTS=$3
CAT_FILE=$4
CAT_BG_FOLDER=$5
TISSUED=$6
ARRAY_LIST=$7



# Array Dependent Config
# If there is no array number assigned:
#   - that means that the script is being executed directly. 
#   - it takes the first item of the list as a parameter
# If there is an array number assigned:
#   - that means that the script is being executed through the cluster
#   - it takes the item correspondent to the array number
if [ "${SLURM_ARRAY_TASK_ID}" == "" ]
then
  IDD=$(cat ${ARRAY_LIST} | sed -n 1p )
else
 IDD=$(cat ${ARRAY_LIST} | sed -n ${SLURM_ARRAY_TASK_ID}p)
fi 


CATS=($( cat ${CAT_FILE}  ))


for CAT in ${CATS[@]}

do

echo Working ... ${CAT}



GENEVAR_LS=${INPUTS}/${TISSUED}/${IDD}/*.genevar



#BG_FILE=${CAT_BG_FOLDER}/everyones.${VERSION}.genes


BG_FILE=${CAT_BG_FOLDER}/${TISSUED}.${IDD}.${CAT}.compound.bg


# this is not getting autosomal 
cat $(ls ${GENEVAR_LS} | grep "rare_${CAT}\." ) | cut -f1 | grep -v -w ENSG | grep -f ${BG_FILE} | sort -u > ./real/inputs/${TISSUED}.${IDD}.${CAT}.res.genes

#cat $(ls ${GENEVAR_LS} | grep "rare_${VERSION}\." ) | cut -f1 | grep -v -w ENSG | sort -u > ./real/inputs/${VERSION}.res.genes



RES_FILE=./real/inputs/${TISSUED}.${IDD}.${CAT}.res.genes


# Execution
Rscript ./real/scripts/WebgestaltR.R \
    ${TISSUED} \
    ${IDD} \
    ${BG_FILE} \
    ${RES_FILE} \
    ${CAT} \
    ${GMT} \
    ${OUTPUTS}


done
