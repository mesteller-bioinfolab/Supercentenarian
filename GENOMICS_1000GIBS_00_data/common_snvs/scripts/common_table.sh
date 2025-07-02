#!/bin/bash
#
#SBATCH -p haswell # partition (queue)
#SBATCH -N 1 # number of nodes
#SBATCH -J commontable
#SBATCH --mem 16G # memory pool for all cores
#SBATCH -t 0-00:59# time (D-HH:MM)
#SBATCH -o ./real/outputs/log.%j.out # STDOUT
#SBATCH -e ./real/outputs/log.%j.err # STDERR
#SBATCH --mail-type=END,FAIL # notifications for job done & fail
#SBATCH --mail-user=claudia.vasallo@upf.edu # send-to address



# Modules
#module load VEP/102.0-foss-2019a-Perl-5.28.1
module load BCFtools/1.16-GCC-11.2.0


# Config
ARRAY_LIST=$1
INPUTS_DIR=$2
TMPS_DIR=$3
OUTPUTS_DIR=$4
TISSUED=$5
IDD=$6


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

G1000=${INPUTS_DIR}/1000G.${CHR}.common.snvs
# ./real/tmp/1000G.chr6.common.snvs 

echo CHROM POS REF ALT AF AF_EUR AF_EUR_unrel > ${OUTPUTS_DIR}/1000G.${CHR}.eur.common.table
bcftools query -f '%ID %CHROM %POS %REF %ALT %AF %INFO/AF_EUR  %INFO/AF_EUR_unrel\n' ${G1000} >> ${OUTPUTS_DIR}/1000G.${CHR}.eur.common.table
