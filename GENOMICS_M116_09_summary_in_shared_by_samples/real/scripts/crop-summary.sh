#!/bin/bash
#
#SBATCH --partition=haswell  # haswell-p normal # partition (queue)
#SBATCH -N 1 # number of nodes
#SBATCH -J cropsummary
#SBATCH --mem 32G # memory pool for all cores
#SBATCH -t 0-03:59# time (D-HH:MM)
#SBATCH -o ./real/outputs/log.%j.out # STDOUT
#SBATCH -e ./real/outputs/log.%j.err # STDERR
#SBATCH --mail-type=END,FAIL # notifications for job done & fail
#SBATCH --mail-user=claudia.vasallo@upf.edu # send-to address



# Modules

module load R


# Config
ARRAY_LIST=$1
INPUTS_DIR=$2
TMPS_DIR=$3
OUTPUTS_DIR=$4
TISSUED=$5
IDD=$6
SUFFIX=$7
FILTER=$8


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


# Run
#TABLE=${INPUTS_DIR}/${TISSUED}/${IDD}/${TISSUED}.${IDD}.${CHR}.filtered.table
#Rscript ./real/scripts/summary.R ${INPUTS_DIR} ${OUTPUTS_DIR} ${TISSUED} ${IDD} ${CHR}
# crop fullsplit, then run summary on it to get metrics for them instead of full 6-annots 
# cat /gpfs42/projects/lab_anavarro/disease_pleiotropies/supercentenarian/7-summary/real/outputs/blood/M116/blood.M116.chr22.fullsplit_AF



cat ${INPUTS_DIR}/${TISSUED}/${IDD}/${TISSUED}.${IDD}.${CHR}.fullsplit_AF | head -n1 > ${OUTPUTS_DIR}/${TISSUED}/${IDD}/${TISSUED}.${IDD}.${CHR}.fullsplit_AF2 
cat ${INPUTS_DIR}/${TISSUED}/${IDD}/${TISSUED}.${IDD}.${CHR}.fullsplit_AF | grep -f ${TMPS_DIR}/M116.${CHR}.var.rep >> ${OUTPUTS_DIR}/${TISSUED}/${IDD}/${TISSUED}.${IDD}.${CHR}.fullsplit_AF2

# once the cropped version of fullsplit is in real/outputs, inputs dir for R script is outputs dir
#Rscript ./real/scripts/summary2.R ${OUTPUTS_DIR} ${OUTPUTS_DIR} ${TISSUED} ${IDD} ${CHR}

Rscript ../08-filter-category-summary/real/scripts/summary1b.R ${OUTPUTS_DIR} ${OUTPUTS_DIR} ${TISSUED} ${IDD} ${CHR} ${SUFFIX} ${FILTER}

