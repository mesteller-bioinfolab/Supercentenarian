#!/bin/bash
#
#SBATCH --partition=haswell       # normal  # haswell-p normal # partition (queue)
#SBATCH -N 1 # number of nodes
#SBATCH -J tab-summary
#SBATCH --mem 8G # memory pool for all cores
#SBATCH -t 0-00:29# time (D-HH:MM)
#SBATCH -o ./real/outputs/log.%j.out # STDOUT
#SBATCH -e ./real/outputs/log.%j.err # STDERR
#SBATCH --mail-type=END,FAIL # notifications for job done & fail
#SBATCH --mail-user=claudia.vasallo@upf.edu # send-to address


# Modules

#module load R
#module load R/3.6.0-foss-2018b



INPUTS=$1
OUTPUTS=$2
TISSUED=$3
ARRAY_LIST=$4
CAT_FILE=$5


# Array Dependent Config
# If there is no array number assigned:
#   - that means that the script is being executed directly. 
#   - it takes the first item of the list as a parameter
# If there is an array number assigned:
#   - that means that the script is being executed through the cluster
#   - it takes the item correspondent to the array number
if [ "${SLURM_ARRAY_TASK_ID}" == "" ]
then
  IDD=$(cat ${ARRAY_LIST} | sed -n 1p)
else
  IDD=$(cat ${ARRAY_LIST} | sed -n ${SLURM_ARRAY_TASK_ID}p)
fi 


CATEGORIES=$(cat ${CAT_FILE})

# Per categories, all chromosomes together



echo Category Variants Genes > ${OUTPUTS}/${TISSUED}.${IDD}.numbers.table.txt
# All eur rare variants
# Variants -----
cat ${INPUTS}/${TISSUED}/${IDD}/*.eur_rare_all.fulltable | cut -f3 | grep -v -w ID | sort -u > ${OUTPUTS}/${TISSUED}.${IDD}.eur_rare.variants


# Genes -----
cat ${INPUTS}/${TISSUED}/${IDD}/*.rare_all.genevar | cut -f1 | grep -v -w ENSG | sort -u > ${OUTPUTS}/${TISSUED}.${IDD}.eur_rare.genes

# CSVS
# Variants -----
#cat $(ls real/outputs/blood/M116/*.higher_healthy.tsv )| cut -f193 | sort -u | grep -v -w ID  > ${OUTPUTS}/higher_healthy_rare.variants
#cat $(ls real/outputs/blood/M116/*.lower_healthy.tsv )| cut -f193 | sort -u | grep -v -w ID  > ${OUTPUTS}/lower_healthy_rare.variants
# Genes -----
#cat $(ls real/outputs/blood/M116/*.higher_healthy.tsv )| cut -f15,14 | sort -u | grep -v -w ID  > ${OUTPUTS}/higher_healthy_rare.genes
#cat $(ls real/outputs/blood/M116/*.lower_healthy.tsv )| cut -f15,14 | sort -u | grep -v -w ID  > ${OUTPUTS}/lower_healthy_rare.genes
#cat ${OUTPUTS}/higher_healthy_rare.genes > ${OUTPUTS}/diff_rare.genes
#cat  ${OUTPUTS}/lower_healthy_rare.genes >> ${OUTPUTS}/diff_rare.genes
#DIFF=${OUTPUTS}/diff_rare.genes
#echo ${CATEGORY} $(cat ${OUTPUTS}/${TISSUED}.${IDD}.rare.variants | wc -l) $(cat ${OUTPUTS}/rare.genes | wc -l) $(cat ${OUTPUTS}/higher_healthy_rare.variants | wc -l ) $(cat ${OUTPUTS}/lower_healthy_rare.variants | wc -l ) $(cat ${DIFF} | sort -u | wc -l) >> ${OUTPUTS}/numbers.table.txt


echo all $(cat ${OUTPUTS}/${TISSUED}.${IDD}.eur_rare.variants | wc -l) $(cat ${OUTPUTS}/${TISSUED}.${IDD}.eur_rare.genes | wc -l) >> ${OUTPUTS}/${TISSUED}.${IDD}.numbers.table.txt


# eur_rare_altering.variants
# Variants of interest

#CATEGORIES=( damaging moderate modifier_cd modifier_nc altering low_only )


for CATEGORY in ${CATEGORIES[@]}
do

# All 
# Variants -----
cat ${INPUTS}/${TISSUED}/${IDD}/*.eur_rare_${CATEGORY}.variants | sort -u > ${OUTPUTS}/${TISSUED}.${IDD}.eur_rare_${CATEGORY}.variants
# Genes -----
cat ${INPUTS}/${TISSUED}/${IDD}/*.rare_${CATEGORY}.genevar | cut -f1 | grep -v -w ENSG | sort -u > ${OUTPUTS}/${TISSUED}.${IDD}.eur_rare_${CATEGORY}.genes

# CSVS diff
# Variants -----
#cat $(ls real/outputs/blood/M116/*.higher_healthy.tsv | grep ${CATEGORY} ) | cut -f193 | sort -u | grep -v -w ID  > ${OUTPUTS}/higher_healthy_rare_${CATEGORY}.variants
#cat $(ls real/outputs/blood/M116/*.lower_healthy.tsv | grep ${CATEGORY} ) | cut -f193 | sort -u | grep -v -w ID  > ${OUTPUTS}/lower_healthy_rare_${CATEGORY}.variants


# Genes -----
#cat $(ls real/outputs/blood/M116/*.higher_healthy.tsv | grep ${CATEGORY} ) | cut -f15,14 | sort -u | grep -v -w ID  > ${OUTPUTS}/higher_healthy_rare_${CATEGORY}.genes
#cat $(ls real/outputs/blood/M116/*.lower_healthy.tsv | grep ${CATEGORY} ) | cut -f15,14 | sort -u | grep -v -w ID  > ${OUTPUTS}/lower_healthy_rare_${CATEGORY}.genes
#cat ${OUTPUTS}/higher_healthy_rare_${CATEGORY}.genes > ${OUTPUTS}/diff_rare_${CATEGORY}.genes
#cat ${OUTPUTS}/lower_healthy_rare_${CATEGORY}.genes >> ${OUTPUTS}/diff_rare_${CATEGORY}.genes
#DIFFCAT=${OUTPUTS}/diff_rare_${CATEGORY}.genes
#echo ${CATEGORY} $(cat ${OUTPUTS}/eur_rare_${CATEGORY}.variants | wc -l) $(cat ${OUTPUTS}/rare_${CATEGORY}.genes | wc -l) $(cat ${OUTPUTS}/higher_healthy_rare_${CATEGORY}.variants | wc -l ) $(cat ${OUTPUTS}/lower_healthy_rare_${CATEGORY}.variants | wc -l ) $(cat ${DIFFCAT} | sort -u | wc -l) >> ${OUTPUTS}/numbers.table.txt

echo ${CATEGORY} $(cat ${OUTPUTS}/${TISSUED}.${IDD}.eur_rare_${CATEGORY}.variants | wc -l) $(cat ${OUTPUTS}/${TISSUED}.${IDD}.eur_rare_${CATEGORY}.genes | wc -l) >> ${OUTPUTS}/${TISSUED}.${IDD}.numbers.table.txt

done


# extra too cats

CATEGORY=non_coding
cat ${INPUTS}/${TISSUED}/${IDD}/*.eur_rare_non_coding.genes > ${OUTPUTS}/${TISSUED}.${IDD}.eur_rare_${CATEGORY}.genes


CATEGORY=coding_with_no

cat ${INPUTS}/${TISSUED}/${IDD}/*.eur_rare_non_coding.genes > ${OUTPUTS}/${TISSUED}.${IDD}.eur_rare_${CATEGORY}.genes

