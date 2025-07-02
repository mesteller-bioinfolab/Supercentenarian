#!/bin/bash
#
#SBATCH --partition=normal   #haswell       # normal  # haswell-p normal # partition (queue)
#SBATCH -N 1 # number of nodes
#SBATCH -J tab-summary
#SBATCH --mem 8G # memory pool for all cores
#SBATCH -t 0-00:29# time (D-HH:MM)
#SBATCH -o ./real/summary_tables/log.%j.out # STDOUT
#SBATCH -e ./real/summary_tables/log.%j.err # STDERR
#SBATCH --mail-type=END,FAIL # notifications for job done & fail
#SBATCH --mail-user=claudia.vasallo@upf.edu # send-to address




INPUTS=$1
OUTPUTS=$2
IDD=$3
TISSUED=$4
CATFILE=$5


CATEGORIES=( $(cat ${CATFILE[@]} ) )


echo ${CATEGORIES}

# Per categories, all chromosomes together



echo Category Variants Genes  > ${OUTPUTS}/${TISSUED}.${IDD}.numbers.table.txt

# All eur rare variants
# Variants -----
cat ${INPUTS}/*/M116/*.eur_rare_all.fulltable | cut -f3 | grep -v -w ID | sort -u > ${OUTPUTS}/${TISSUED}.${IDD}.eur_rare.variants


# Genes -----
cat ${INPUTS}/*/M116/*.rare_all.genevar | cut -f1 | grep -v -w ENSG | sort -u > ${OUTPUTS}/${TISSUED}.${IDD}.eur_rare.genes



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
#echo ${CATEGORY} $(cat ${OUTPUTS}/rare.variants | wc -l) $(cat ${OUTPUTS}/rare.genes | wc -l) $(cat ${OUTPUTS}/higher_healthy_rare.variants | wc -l ) $(cat ${OUTPUTS}/lower_healthy_rare.variants | wc -l ) $(cat ${DIFF} | sort -u | wc -l) >> ${OUTPUTS}/numbers.table.txt





echo "all"  $(cat ${OUTPUTS}/${TISSUED}.${IDD}.eur_rare.variants | wc -l) $(cat ${OUTPUTS}/${TISSUED}.${IDD}.eur_rare.genes | wc -l) >> ${OUTPUTS}/${TISSUED}.${IDD}.numbers.table.txt




# Variants of interest




for CATEGORY in ${CATEGORIES[@]}

do


# All 
# Variants -----
cat ${INPUTS}/*/M116/*.eur_rare_${CATEGORY}.variants | sort -u > ${OUTPUTS}/${TISSUED}.${IDD}.eur_rare_${CATEGORY}.variants


# Genes -----
cat ${INPUTS}/*/M116/*.rare_${CATEGORY}.genevar | cut -f1 | grep -v -w ENSG | sort -u > ${OUTPUTS}/${TISSUED}.${IDD}.eur_rare_${CATEGORY}.genes






# CSVS diff
# Variants -----
#cat $(ls real/outputs/*/M116/*.higher_healthy.tsv | grep ${CATEGORY} ) | cut -f193 | sort -u | grep -v -w ID  > ${OUTPUTS}/higher_healthy_rare_${CATEGORY}.variants
#cat $(ls real/outputs/*/M116/*.lower_healthy.tsv | grep ${CATEGORY} ) | cut -f193 | sort -u | grep -v -w ID  > ${OUTPUTS}/lower_healthy_rare_${CATEGORY}.variants
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
cat ${INPUTS}/*/${IDD}/*.eur_rare_${CATEGORY}.genes | sort -u > ${OUTPUTS}/${TISSUED}.${IDD}.eur_rare_${CATEGORY}.genes


CATEGORY=coding_with_nc

cat ${INPUTS}/*/${IDD}/*.eur_rare_${CATEGORY}.genes | sort -u > ${OUTPUTS}/${TISSUED}.${IDD}.eur_rare_${CATEGORY}.genes

