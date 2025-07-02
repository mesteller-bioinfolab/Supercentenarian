#!/bin/bash
#
#SBATCH -p haswell # partition (queue)
#SBATCH -N 1 # number of nodes
#SBATCH -J mz_tot
#SBATCH --mem 60G # memory pool for all cores
#SBATCH -t 0-03:00 # time (D-HH:MM)
#SBATCH -o ./real/tmp/log.%j.out # STDOUT
#SBATCH -e ./real/tmp/log.%j.err # STDERR
#SBATCH --mail-type=FAIL # notifications for job done & fail
#SBATCH --mail-user=claudia.vasallo@upf.edu # send-to address


# Config
INPUTS=$1
OUTPUTS=$2
CATEGORY=$3
ARRAY_LIST=$4



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




ls ../../../Olot/09-summary-in-shared-by-samples/real/outputs/*/M116/*.${CHR}.eur_rare_${CATEGORY}.variants.dput > ${INPUTS}/${CATEGORY}.${CHR}.dputs.list

ls ../../../1000GIBS/08-filter-category-summary/real/outputs/1000GIBS/*/*.${CHR}.eur_rare_${CATEGORY}.variants.dput >> ${INPUTS}/${CATEGORY}.${CHR}.dputs.list



#DPUTS_LIST=$(cat ${INPUTS}/${CATEGORY}.${CHR}.dputs.list)

DPUTS_LIST=${INPUTS}/${CATEGORY}.${CHR}.dputs.list

#DPUT_LIST1=${INPUTS}/${CATEGORY}.${CHR}.dputs.list1
#DPUTS_LIST=${INPUTS}/${CATEGORY}.${CHR}.dputs.list2



module load R

Rscript ./real/scripts/tot_var.R ${CATEGORY} ${CHR} ${OUTPUTS} ${DPUTS_LIST} 

#${DPUTS_LIST2}

