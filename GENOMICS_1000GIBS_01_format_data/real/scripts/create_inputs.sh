#!/bin/bash
#
#SBATCH -p haswell # partition (queue)
#SBATCH -N 1 # number of nodes
#SBATCH -J inputs
#SBATCH --mem 16 # memory pool for all cores
#SBATCH -t 0-23:59# time (D-HH:MM)
#SBATCH -o ./real/outputs/log.%j.out # STDOUT
#SBATCH -e ./real/outputs/log.%j.err # STDERR
#SBATCH --mail-type=END,FAIL # notifications for job done & fail
#SBATCH --mail-user=claudia.vasallo@upf.edu # send-to address



# Modules
module load BCFtools/1.16-GCC-11.2.0

# Config
ARRAY_LIST=$1
INPUTS_DIR=$2
OUTPUTS_DIR=$3
TMPS_DIR=$4
SAMPLE_FILE=$5


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


# Filter to sample ids: IBS Females
# 20201028_CCDG_14151_B01_GRM_WGS_2020-08-05_chrX.recalibrated_variants.annotated.vcf.gz
#bcftools view -S ${SAMPLE_FILE} ${INPUTS_DIR}/20201028_CCDG_14151_B01_GRM_WGS_2020-08-05_${CHR}.recalibrated_variants.annotated.vcf.gz --force-samples \
# v2 for chrX
bcftools view -S ${SAMPLE_FILE} ${INPUTS_DIR}/CCDG_14151_B01_GRM_WGS_2020-08-05_${CHR}.filtered.eagle2-phased.vcf.gz --force-samples \
> ${TMPS_DIR}/ibs_females_${CHR}.vcf

# Quality filter matching M116 
bcftools view -e 'QD < 2.0 || FS > 60.0 || MQ < 40.0  || MQRankSum < -12.5 || ReadPosRankSum < -8.0' ${TMPS_DIR}/ibs_females_${CHR}.vcf > ${TMPS_DIR}/ibs_females_${CHR}.vcf.filtered

bgzip -c ${TMPS_DIR}/ibs_females_${CHR}.vcf.filtered > ${TMPS_DIR}/ibs_females_${CHR}.vcf.gz
tabix -p vcf ${TMPS_DIR}/ibs_females_${CHR}.vcf.gz


# -----
# Split into single-sample VCFs - one to many

## For one CHR, split into IBS samples
SPLITDIR=${TMPS_DIR}/splitdir

mkdir ${SPLITDIR}/${CHR}
VCF=${TMPS_DIR}/ibs_females_${CHR}.vcf.gz


sed "s/$/\t${CHR}/" ${SAMPLE_FILE} > ${TMPS_DIR}/${CHR}.cols
cat ${TMPS_DIR}/${CHR}.cols | cut -f2 > ${TMPS_DIR}/${CHR}.col
rm -rf ${TMPS_DIR}/${CHR}.cols
paste ${SAMPLE_FILE} ${TMPS_DIR}/${CHR}.col | tr '\t' '_' >	${TMPS_DIR}/${CHR}.col_
paste ${SAMPLE_FILE} ${TMPS_DIR}/${CHR}.col_ > ${SAMPLE_FILE}.${CHR}


# Split into single-sample VCFSs so the M116 pipeline works as it is 
bcftools +split ${VCF} -S ${SAMPLE_FILE}.${CHR} -Oz -o ${SPLITDIR}/${CHR}

# Now, for each sample, create proper folder and save VCF there
IDS=($( ls ${SPLITDIR}/${CHR} | cut -d '.' -f1 | cut -d '_' -f1 ))

for IDD in ${IDS[@]}
do
#---Move to proper destination
mkdir ${OUTPUTS_DIR}/${IDD}
mv ${SPLITDIR}/${CHR}/${IDD}_${CHR}.vcf.gz ${OUTPUTS_DIR}/${IDD}/${IDD}.${CHR}.vcf.gz_
##--Select variants only (exclude genotypes 0|0)
#bcftools view -c1 ${OUTPUTS_DIR}/${IDD}/${IDD}.${CHR}.vcf.gz_ > ${OUTPUTS_DIR}/${IDD}/${IDD}.${CHR}.vcf1
bcftools view -i 'GT[*]="alt"' ${OUTPUTS_DIR}/${IDD}/${IDD}.${CHR}.vcf.gz_  > ${OUTPUTS_DIR}/${IDD}/${IDD}.${CHR}.vcf

bgzip -c ${OUTPUTS_DIR}/${IDD}/${IDD}.${CHR}.vcf > ${OUTPUTS_DIR}/${IDD}/${IDD}.${CHR}.vcf.gz
tabix -p vcf ${OUTPUTS_DIR}/${IDD}/${IDD}.${CHR}.vcf.gz

yes | rm ${OUTPUTS_DIR}/${IDD}/${IDD}.${CHR}.vcf
yes | rm ${OUTPUTS_DIR}/${IDD}/${IDD}.${CHR}.vcf.gz_

done
