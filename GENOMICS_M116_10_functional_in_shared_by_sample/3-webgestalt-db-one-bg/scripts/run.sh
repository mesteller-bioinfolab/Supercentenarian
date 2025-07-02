#!/bin/bash
#
#SBATCH -p haswell # partition (queue)
#SBATCH -N 1 # number of nodes
#SBATCH -J webgesalt
#SBATCH --mem 60G # memory pool for all cores
#SBATCH -t 0-03:00 # time (D-HH:MM)
#SBATCH -o ./real/outputs/log.%j.out # STDOUT
#SBATCH -e ./real/outputs/log.%j.err # STDERR
#SBATCH --mail-type=FAIL # notifications for job done & fail
#SBATCH --mail-user=claudia.vasallo@upf.edu # send-to address

# Modules
module load R


# Config
INPUTS=$1
OUTPUTS=$2
ARRAY_LIST=$3
CAT_BG_FOLDER=$4




# Array Dependent Config
# If there is no array number assigned:
#   - that means that the script is being executed directly. 
#   - it takes the first item of the list as a parameter
# If there is an array number assigned:
#   - that means that the script is being executed through the cluster
#   - it takes the item correspondent to the array number
if [ "${SLURM_ARRAY_TASK_ID}" == "" ]
then
  VERSION=$(cat ${ARRAY_LIST} | sed -n 1p )
else
  VERSION=$(cat ${ARRAY_LIST} | sed -n ${SLURM_ARRAY_TASK_ID}p)
fi 




echo Working ... ${VERSION}


GENEVAR_LS=${INPUTS}/*/M116/*.genevar
BG_FILE=${CAT_BG_FOLDER}/everyones.${VERSION}.genes


#BG=$(cat $(ls ${GENEVAR_LS} | grep all) | cut -f1 | grep -v -w ENSG | grep -f ${BACKGROUND_GENES} | sort -u) 


cat $(ls ${GENEVAR_LS} | grep ${VERSION}) | cut -f1 | grep -v -w ENSG | grep -f ${BG_FILE} | sort -u > ./real/inputs/${VERSION}.res.genes




#echo $BG | tr ' ' '\n' > ./real/inputs/bg.genes
#echo $RES | tr ' ' '\n' > ./real/inputs/${VERSION}.res.genes
#BG_FILE=./real/inputs/bg.genes
#BG_FILE=./real/inputs/M116_coding.genes.txt





RES_FILE=./real/inputs/${VERSION}.res.genes


# Execution
Rscript ./real/scripts/webgesalt.R \
    ${BG_FILE} \
    ${RES_FILE} \
    ${VERSION}


