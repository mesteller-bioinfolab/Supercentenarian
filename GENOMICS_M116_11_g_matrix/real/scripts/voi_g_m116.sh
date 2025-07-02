#!/bin/bash
#
#SBATCH --partition=haswell  # haswell-p normal # partition (queue)
#SBATCH -N 1 # number of nodes
#SBATCH -J chrsjoin
#SBATCH --mem 32G # memory pool for all cores
#SBATCH -t 0-02:59# time (D-HH:MM)
#SBATCH -o ./real/tmp/log.%j.out # STDOUT
#SBATCH -e ./real/tmp/log.%j.err # STDERR
#SBATCH --mail-type=END,FAIL # notifications for job done & fail
#SBATCH --mail-user=claudia.vasallo@upf.edu # send-to address



# Modules

module load R


# Config
ARRAY_LIST=$1
INPUTS=$2   
OUTPUTS=$3
IDD=$4
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
  CHR=$(cat ${ARRAY_LIST} | sed -n 1p)
else
  CHR=$(cat ${ARRAY_LIST} | sed -n ${SLURM_ARRAY_TASK_ID}p)
fi 



CATEGORIES=$( cat ${CAT_FILE} )


EVERYONES=./real/inputs/everyones.${CHR}.genes


for CAT in ${CATEGORIES[@]}


do

cat header  > ./real/inputs/2TISSUES.M116.${CHR}.eur_rare_${CAT}.fulltable


if [ "${CAT}" == "all" ]
then
cat ../09-summary-in-shared-by-samples/real/outputs/*/M116/*.${CHR}.eur_rare_all.fulltable | sort -u >> ./real/inputs/2TISSUES.M116.${CHR}.eur_rare_${CAT}.fulltable
else
cat ../09-summary-in-shared-by-samples/real/outputs/*/M116/*.${CHR}.eur_rare_all.fulltable | grep  -w "${CAT}" | sort -u   >>  ./real/inputs/2TISSUES.M116.${CHR}.eur_rare_${CAT}.fulltable
fi

FULL=./real/inputs/2TISSUES.M116.${CHR}.eur_rare_${CAT}.fulltable

Rscript ./real/scripts/genes_matrix_new.R  ${INPUTS} ${OUTPUTS} ${IDD} ${CHR} ${EVERYONES} ${CAT} ${FULL}




done
