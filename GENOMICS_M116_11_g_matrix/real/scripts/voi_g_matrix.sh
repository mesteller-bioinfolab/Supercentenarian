#!/bin/bash
#
#SBATCH --partition=haswell  # haswell-p normal # partition (queue)
#SBATCH -N 1 # number of nodes
#SBATCH -J chrsjoin
#SBATCH --mem 32G # memory pool for all cores
#SBATCH -t 0-12:59# time (D-HH:MM)
#SBATCH -o ./real/tmp/log.%j.out # STDOUT
#SBATCH -e ./real/tmp/log.%j.err # STDERR
#SBATCH --mail-type=END,FAIL # notifications for job done & fail
#SBATCH --mail-user=claudia.vasallo@upf.edu # send-to address



# Modules

module load R


# Config
ARRAY_LIST=$1
INPUTS_DIR1=$2
INPUTS_DIR2=$3
OUTPUTS=$4
CHR_LIST=$5
CAT_FILE=$6


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



if [ "${IDD}" == "M116" ]
then
    INPUTS=${INPUTS_DIR1}
else
    INPUTS=${INPUTS_DIR2}
fi


CHRS=$(cat ${CHR_LIST})



# ----------------------------------------------------------------------------------------------------------------------------

for CHR in ${CHRS}

do

# gene-level matrix
# genotypes are count voi from genevar of specified version

# Get genes to analyse - done in 0-

if [ "${IDD}" == "M116" ]
then
wait 0
else
# IBS individual in list
cat $(ls ../../1000GIBS/08-filter-category-summary/real/outputs/1000GIBS/${IDD}/*.genevar | grep -w ${CHR} |  grep all ) | head -n1 > ./real/inputs/${IDD}.${CHR}.genevar 
cat $(ls ../../1000GIBS/08-filter-category-summary/real/outputs/1000GIBS/${IDD}/*.genevar | grep -w ${CHR} |  grep all ) | grep -v -w ENSG >> ./real/inputs/${IDD}.${CHR}.genevar 
cat ./real/inputs/${IDD}.${CHR}.genevar | cut -f1 | sort -u | grep -v -w ENSG  > ./real/inputs/${IDD}.${CHR}.genes

fi


# Individuals VOI genes
# do not use this, its not optimal
GENEVAR_FILE=./real/inputs/${IDD}.${CHR}.genevar


# For everyones' genes
EVERYONES=./real/inputs/everyones.${CHR}.genes



Rscript ./real/scripts/genes_matrix.R ${INPUTS} ${OUTPUTS} ${IDD} ${CHR} ${EVERYONES} ${CAT_FILE} ${GENEVAR_FILE}




done
