#!/bin/bash
#
#SBATCH -p haswell # partition (queue)
#SBATCH -N 1 # number of nodes
#SBATCH -J listfilter
#SBATCH --mem 16G # memory pool for all cores
#SBATCH -t 0-06:59# time (D-HH:MM)
#SBATCH -o ./real/outputs/log.%j.out # STDOUT
#SBATCH -e ./real/outputs/log.%j.err # STDERR
#SBATCH --mail-type=END,FAIL # notifications for job done & fail
#SBATCH --mail-user=claudia.vasallo@upf.edu # send-to address



# Modules
module load BCFtools/1.16-GCC-11.2.0


# Config
INPUTS=$1
OUTPUTS=$2
ARRAY_LIST=$3
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



# filter -e # other script to not confude INPUTS DIRS


bcftools annotate --set-id '%CHROM:%POS:%REF:%ALT' ${INPUTS}/${TISSUED}/${IDD}/${TISSUED}.${IDD}.${CHR}.vep  > ${INPUTS}/${TISSUED}/${IDD}/${TISSUED}.${IDD}.${CHR}.vep.id 

bcftools view -e 'ID=@snpsexcl1' ${INPUTS}/${TISSUED}/${IDD}/${TISSUED}.${IDD}.${CHR}.vep.id > ${OUTPUTS}/${TISSUED}/${IDD}/${TISSUED}.${IDD}.${CHR}.filtered.vcf

bcftools view --types=snps ${OUTPUTS}/${TISSUED}/${IDD}/${TISSUED}.${IDD}.${CHR}.filtered.vcf > ${OUTPUTS}/${TISSUED}/${IDD}/${TISSUED}.${IDD}.${CHR}.filtered.vcf.snv
