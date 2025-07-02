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



echo Category Variants Genes  > ${OUTPUTS}/${TISSUED}.${IDD}.autosomal.numbers.table.txt

# All eur rare variants
# Variants -----
cat $( ls ${INPUTS}/*/M116/*.eur_rare_all.fulltable | grep -v chrM | grep -v chrY | grep -v chrX ) | cut -f3 | grep -v -w ID | sort -u > ${OUTPUTS}/${TISSUED}.${IDD}.eur_rare.autosomal.variants


# Genes -----
cat $(ls ${INPUTS}/*/M116/*rare_all.genevar | grep -v chrM | grep -v chrY | grep -v chrX ) | cut -f1 | grep -v -w ENSG | sort -u > ${OUTPUTS}/${TISSUED}.${IDD}.eur_rare.autosomal.genes




echo "all"  $(cat ${OUTPUTS}/${TISSUED}.${IDD}.eur_rare.autosomal.variants | wc -l) $(cat ${OUTPUTS}/${TISSUED}.${IDD}.eur_rare.autosomal.genes | wc -l) >> ${OUTPUTS}/${TISSUED}.${IDD}.autosomal.numbers.table.txt




# Variants of interest




for CATEGORY in ${CATEGORIES[@]}

do


# All 
# Variants -----
cat $(ls ${INPUTS}/*/M116/*.eur_rare_${CATEGORY}.variants |  grep -v chrM | grep -v chrY | grep -v chrX ) | sort -u > ${OUTPUTS}/${TISSUED}.${IDD}.eur_rare_${CATEGORY}.autosomal.variants


# Genes -----
cat $( ls ${INPUTS}/*/M116/*_${CATEGORY}.genevar | grep -v chrM | grep -v chrY | grep -v chrX ) | cut -f1 | grep -v -w ENSG | sort -u > ${OUTPUTS}/${TISSUED}.${IDD}.eur_rare_${CATEGORY}.autosomal.genes





echo ${CATEGORY} $(cat ${OUTPUTS}/${TISSUED}.${IDD}.eur_rare_${CATEGORY}.autosomal.variants | wc -l) $(cat ${OUTPUTS}/${TISSUED}.${IDD}.eur_rare_${CATEGORY}.autosomal.genes | wc -l) >> ${OUTPUTS}/${TISSUED}.${IDD}.autosomal.numbers.table.txt


done


# extra too cats


CATEGORY=non_coding
cat $(ls ${INPUTS}/*/${IDD}/*.eur_rare_${CATEGORY}.genes | grep -v chrM | grep -v chrX | grep -v chrY ) | sort -u > ${OUTPUTS}/${TISSUED}.${IDD}.eur_rare_${CATEGORY}.autosomal.genes


CATEGORY=coding_with_nc

cat $( ls ${INPUTS}/*/${IDD}/*.eur_rare_${CATEGORY}.genes | grep -v chrM | grep -v chrX | grep -v chrY)  | sort -u  > ${OUTPUTS}/${TISSUED}.${IDD}.eur_rare_${CATEGORY}.autosomal.genes

