#!/bin/bash
#
#SBATCH -p normal # partition (queue)
#SBATCH -N 1 # number of nodes
#SBATCH -J vep2table
#SBATCH --mem 32G # memory pool for all cores
#SBATCH -t 9-23:59# time (D-HH:MM)
#SBATCH -o ./real/outputs/log.%j.out # STDOUT
#SBATCH -e ./real/outputs/log.%j.err # STDERR
#SBATCH --mail-type=END,FAIL # notifications for job done & fail
#SBATCH --mail-user=claudia.vasallo@upf.edu # send-to address



# Modules
module load GATK/4.1.6.0-GCCcore-8.3.0-Java-11


# Config
ARRAY_LIST=$1
INPUTS_DIR=$2
OUTPUTS_DIR=$3
TISSUED=$4
IDD=$5


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



# vep 2 to filtered

# Run
#VEP=${INPUTS_DIR}/${TISSUED}/${IDD}/${TISSUED}.${IDD}.${CHR}.vep2
VEP=${INPUTS_DIR}/${TISSUED}/${IDD}/${TISSUED}.${IDD}.${CHR}.vep

# this command was used for CSVS annotated vep files, GT was missing
#gatk VariantsToTable -V $VEP -O ${OUTPUTS_DIR}/${TISSUED}/${IDD}/${TISSUED}.${IDD}.${CHR}.filtered.table -F CHROM -F POS -F ID -F REF -F ALT -F QUAL -F FILTER -F GT -F AF -F DP -F QD -F ANN -F CSVS --split-multi-allelic


# this one is for ANN and CADD annotated vep files
#gatk VariantsToTable -V ${VEP} -O ${OUTPUTS_DIR}/${TISSUED}/${IDD}/${TISSUED}.${IDD}.${CHR}.table -F CHROM -F POS -F ID -F REF -F ALT -F QUAL -F FILTER -F GT -F AF -F DP -F QD -F ANN -F CADD
gatk VariantsToTable -V ${VEP} -O ${OUTPUTS_DIR}/${TISSUED}/${IDD}/${TISSUED}.${IDD}.${CHR}.table -F CHROM -F POS -F ID -F REF -F ALT -F ANN -F CADD -GF GT 
