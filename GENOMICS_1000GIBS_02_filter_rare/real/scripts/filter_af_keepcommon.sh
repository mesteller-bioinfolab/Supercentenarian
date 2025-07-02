#!/bin/bash
#
#SBATCH -p haswell # partition (queue)
#SBATCH -N 1 # number of nodes
#SBATCH -J vepfilter
#SBATCH --mem 80G # memory pool for all cores
#SBATCH -t 9-23:59# time (D-HH:MM)
#SBATCH -o ./real/outputs_common/log.%j.out # STDOUT
#SBATCH -e ./real/outputs_common/log.%j.err # STDERR
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


#TMP_DIR=./real/tmp
#mkdir ${TMPS_DIR}/${TISSUED}/${IDD}

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




# EUR AF filters



# "EUR_AF < 0.01 or not EUR_AF"
#filter_vep -i ../1-vep-annot/real/outputs/blood/M116/blood.M116.chr1.vep -o filtered.blood.M116.g1000.chr1 --force_overwrite --filter "EUR_AF < 0.01 or not EUR_AF" --format vcf --vcf_info_field ANN --only_matched
# "gnomAD_NFE_AF < 0.01 or not gnomAD_NFE_AF"
#filter_vep -i ../1-vep-annot/real/outputs/blood/M116/blood.M116.chr1.vep -o filtered.blood.M116.gnomad.chr1 --force_overwrite --filter "gnomAD_NFE_AF < 0.01 or not gnomAD_NFE_AF" --format vcf --vcf_info_field ANN --only_matched



# gnomAD_NFE_AF < 0.01 or EUR_AF < 0.01 or (not gnomAD_NFE_AF and not EUR_AF)
#filter_vep -i ${INPUTS_DIR}/${TISSUED}/${IDD}/${TISSUED}.${IDD}.${CHR}.vep -o ${OUTPUTS_DIR}/${TISSUED}/${IDD}/${TISSUED}.${IDD}.${CHR}.filtered --force_overwrite --filter "EUR_AF < 0.01 or gnomAD_NFE_AF < 0.01 or (not gnomAD_NFE_AF and not EUR_AF)" --format vcf --vcf_info_field ANN --only_matched

# on custom annotated vcfs
#filter_vep -i ${INPUTS_DIR}/${TISSUED}/${IDD}/${TISSUED}.${IDD}.${CHR}.vep2 -o ${OUTPUTS_DIR}/${TISSUED}/${IDD}/${TISSUED}.${IDD}.${CHR}.filtered --force_overwrite --filter "EUR_AF < 0.01 or gnomADg_NFE_AF < 0.01 or gnomADe_NFE_AF < 0.01 or (not gnomADg_NFE_AF and not gnomADe_NFE_AF and not EUR_AF)" --format vcf --vcf_info_field ANN --only_matched


# on original 1000G deep coverage
#filter_vep -i ${INPUTS_DIR}/${TISSUED}/${IDD}/${IDD}.${CHR}.vcf.gz -o ${OUTPUTS_DIR}/${TISSUED}/${IDD}/${TISSUED}.${IDD}.${CHR}.filtered --force_overwrite --filter "AF_EUR < 0.015 or not AF_EUR" --format vcf --vcf_info_field INFO --only_matched



bcftools view --types snps ${INPUTS_DIR}/${TISSUED}/${IDD}/${IDD}.${CHR}.vcf.gz > ${TMPS_DIR}/${TISSUED}/${IDD}/${TISSUED}.${IDD}.${CHR}.snvs
bcftools filter -i 'INFO/AF_EUR_unrel>0.015' ${TMPS_DIR}/${TISSUED}/${IDD}/${TISSUED}.${IDD}.${CHR}.snvs > ${OUTPUTS_DIR}/${TISSUED}/${IDD}/${TISSUED}.${IDD}.${CHR}.filtered
gzip ${OUTPUTS_DIR}/${TISSUED}/${IDD}/${TISSUED}.${IDD}.${CHR}.common
