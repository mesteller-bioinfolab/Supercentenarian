#!/bin/bash
#
#SBATCH -p haswell # partition (queue)
#SBATCH -N 1 # number of nodes
#SBATCH -J qc
#SBATCH --mem 32G # memory pool for all cores
#SBATCH -t 0-03:59# time (D-HH:MM)
#SBATCH -o ./log.%j.out # STDOUT
#SBATCH -e ./log.%j.err # STDERR
#SBATCH --mail-type=END,FAIL # notifications for job done & fail
#SBATCH --mail-user=claudia.vasallo@upf.edu # send-to address



# Modules
module load BCFtools/1.16-GCC-11.2.0




# Config
ARRAY_LIST=$1
INPUTS_DIR=$2
TMPS_DIR=$3
OUTPUTS_DIR=$4


TISSUED=1000GIBS


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




THEDIR=${INPUTS_DIR}/${TISSUED}
IDS=($( ls ${THEDIR}))


for IDD in ${IDS[@]}
do

mkdir ${TMPS_DIR}/${TISSUED}/${IDD}
mkdir ${OUTPUTS_DIR}/${TISSUED}/${IDD}


# M116 FILTER : QD < 2.0 || FS > 60.0 || MQ < 40.0 || HaplotypeScore > 13.0 || MappingQualityRankSum < -12.5 || ReadPosRankSum < -8.0

# exclude
#bcftools view -i 'QD < 2.0 || FS > 60.0 || MQ < 40.0  || MQRankSum < -12.5 || ReadPosRankSum < -8.0' ${INPUTS_DIR}/${TISSUED}/${IDD}/${IDD}.${CHR}.vcf.gz > ${TMPS_DIR}/${TISSUED}/${IDD}/${TISSUED}.${IDD}.${CHR}.to_exclude.vcf

# filter
#bcftools view -e 'QD < 2.0 || FS > 60.0 || MQ < 40.0  || MQRankSum < -12.5 || ReadPosRankSum < -8.0' ${TMPS_DIR}/${TISSUED}/${IDD}/${IDD}.${CHR}.vcf.gz  > ${TISSUED}/${IDD}/${TISSUED}.${IDD}.${CHR}.filtered.vcf.gz

# exclude table

echo CHROM POS REF ALT > ${OUTPUTS_DIR}/${TISSUED}/${IDD}/${TISSUED}.${IDD}.${CHR}.to_exclude.table
bcftools query -f '%CHROM %POS %REF %ALT\n' ${TMPS_DIR}/${TISSUED}/${IDD}/${TISSUED}.${IDD}.${CHR}.to_exclude.vcf >> ${OUTPUTS_DIR}/${TISSUED}/${IDD}/${TISSUED}.${IDD}.${CHR}.to_exclude.table




done


