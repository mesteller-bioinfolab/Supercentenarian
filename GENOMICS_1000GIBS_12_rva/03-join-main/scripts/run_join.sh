#!/bin/bash
#
#SBATCH -p normal # partition (queue)
#SBATCH -N 1 # number of nodes
#SBATCH -J rvajoin
#SBATCH --mem 120G # memory pool for all cores
#SBATCH -t 0-03:00 # time (D-HH:MM)
#SBATCH -o ./real/outputs/log.%j.out # STDOUT
#SBATCH -e ./real/outputs/log.%j.err # STDERR
#SBATCH --mail-type=FAIL # notifications for job done & fail
#SBATCH --mail-user=claudia.vasallo@upf.edu # send-to address






# Config
INPUTS=$1
OUTPUTS=$2
TISSUED=$3
ARRAY_LIST=$4
CATS_FILE=$5


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



CATS=( $(cat ${CATS_FILE} ) )


for CAT in ${CATS[@]}

do


cat ${INPUTS}/01-CAST/real/outputs/a5.${TISSUED}.${IDD}_${CAT} > ./real/tmp/a5_inmain.${TISSUED}.${IDD}_${CAT}
cat ${INPUTS}/02-MZ/real/outputs/rvt1_a5.${TISSUED}.${IDD}_${CAT} >> ./real/tmp/a5_inmain.${TISSUED}.${IDD}_${CAT}
cat ./real/tmp/a5_inmain.${TISSUED}.${IDD}_${CAT} | sort -u > ${OUTPUTS}/a5_${CAT}.${TISSUED}.${IDD}.inmain.genes

cat ${INPUTS}/01-CAST/real/outputs/a5.${TISSUED}.${IDD}_${CAT} | grep -f ${INPUTS}/02-MZ/real/outputs/rvt1_a5.${TISSUED}.${IDD}_${CAT} > ${OUTPUTS}/a5_${CAT}.${TISSUED}.${IDD}.in2.genes


done
