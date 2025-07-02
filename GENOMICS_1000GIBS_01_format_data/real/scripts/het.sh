#!/bin/bash
#
#SBATCH -p haswell # partition (queue)
#SBATCH -N 1 # number of nodes
#SBATCH -J het
#SBATCH --mem 100G # memory pool for all cores
#SBATCH -t 9-23:59# time (D-HH:MM)
#SBATCH -o ./real/outputs/log.%j.out # STDOUT
#SBATCH -e ./real/outputs/log.%j.err # STDERR
#SBATCH --mail-type=END,FAIL # notifications for job done & fail
#SBATCH --mail-user=claudia.vasallo@upf.edu # send-to address


TMPS_DIR=./real/tmp


# Modules
module load BCFtools/1.16-GCC-11.2.0




# Config
ARRAY_LIST=$1
INPUTS_DIR=$2
OUTPUTS_DIR=$3
TMPS_DIR=$4



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







echo sample nNonRefHom nHets > HET_${CHR}.txt


SPLITDIR=${TMPS_DIR}/splitdir

IDS=($( ls ${SPLITDIR}/${CHR} | cut -d '.' -f1 | cut -d '_' -f1 ))

for IDD in ${IDS[@]}
do
#VCF=real/outputs/1000GIBS/HG01501/HG01501.chr1.vcf.gz

VCF=${INPUTS_DIR}/${IDD}/${IDD}.${CHR}.vcf.gz
echo $(bcftools stats -s - ${VCF} | grep "^PSC" -B1 | tail -n+2 | cut -f3,5,6) >> HET_${CHR}.txt

done


