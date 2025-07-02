#!/bin/bash
#
#SBATCH --partition=haswell       # normal  # haswell-p normal # partition (queue)
#SBATCH -N 1 # number of nodes
#SBATCH -J fast-summary
#SBATCH --mem 64G # memory pool for all cores
#SBATCH -t 0-14:59# time (D-HH:MM)
#SBATCH -o ./real/outputs/log.%j.out # STDOUT
#SBATCH -e ./real/outputs/log.%j.err # STDERR
#SBATCH --mail-type=END,FAIL # notifications for job done & fail
#SBATCH --mail-user=claudia.vasallo@upf.edu # send-to address



# Modules

module load R
#module load R/3.6.0-foss-2018b


# Config
ARRAY_LIST=$1
INPUTS_DIR=$2
OUTPUTS_DIR=$3
TISSUED=$4
IDS_FILE=$5
SUFFIX=$6
FILTER=$7

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


IDS=$(cat ${IDS_FILE})


# Run
#TABLE=${INPUTS_DIR}/${TISSUED}/${IDD}/${TISSUED}.${IDD}.${CHR}.filtered.table


#Rscript /gpfs42/projects/lab_anavarro/disease_pleiotropies/supercentenarian/Olot/08-filter-category-summary/real/scripts/summary1b.R ${INPUTS_DIR} ${OUTPUTS_DIR} ${TISSUED} ${IDD} ${CHR}


for IDD in ${IDS}
do

mkdir ${OUTPUTS_DIR}/${TISSUED}/${IDD}

#Rscript ../../Olot/08-filter-category-summary/real/scripts/summary1b.R ${INPUTS_DIR} ${OUTPUTS_DIR} ${TISSUED} ${IDD} ${CHR} ${SUFFIX} ${FILTER}

Rscript ./real/scripts/summary1b.R ${INPUTS_DIR} ${OUTPUTS_DIR} ${TISSUED} ${IDD} ${CHR} ${SUFFIX} ${FILTER}

done

