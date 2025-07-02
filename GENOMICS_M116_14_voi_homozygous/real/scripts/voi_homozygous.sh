#!/bin/bash
#
#SBATCH -p haswell # partition (queue)
#SBATCH -N 1 # number of nodes
#SBATCH -J voihomo
#SBATCH --mem 32G # memory pool for all cores
#SBATCH -t 0-00:59# time (D-HH:MM)
#SBATCH -o ./real/outputs/log.%j.out # STDOUT
#SBATCH -e ./real/outputs/log.%j.err # STDERR
#SBATCH --mail-type=END,FAIL # notifications for job done & fail
#SBATCH --mail-user=claudia.vasallo@upf.edu # send-to address



# Modules
module load R


# Config
TMP=$1
OUTPUTS=$2
TISSUED=$3
ARRAY_LIST=$4



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


#CATEGORIES=( damaging moderate modifier_cd modifier_nc altering low low_only rare )

CATEGORIES=( DAMAGING MODERATE_ MODIFIER_CD MODIFIER_NC ALTERING LOW LOW_ONLY NCNC_CADD15 CADD15 )


#if [[ $IDD=="M116" ]] 
#then
#TISSUED=2TISSUES



INPUTS=../09-summary-in-shared-by-samples/real/outputs/*/M116/*.eur_rare_all.fulltable 

#INPUTS=../09-summary-in-shared-by-samples/real/outputs/*/M116/*.fullsplit_AF


#else
#TISSUED=1000GIBS
#INPUTS=../../1000GIBS/08-filter-category-summary/real/outputs/1000GIBS/${IDD}/1000GIBS/${IDD}/*.eur_rare_all.fulltable
#cat ../../1000GIBS/08-filter-category-summary/real/outputs/1000GIBS/HG01501/1000GIBS.HG01501.chr22.eur_rare_all.fulltable | head -n1 > ${TMPINPUTS=../../1000GIBS/08-filter-category-summary/real/outputs/1000GIBS/${IDD}/*.eur_rare_all.fulltable
#cat ../../1000GIBS/08-filter-category-summary/real/outputs/1000GIBS/HG01501/1000GIBS.HG01501.chr22.eur_rare_all.fulltable | head -n1 > ${TMP}/${TISSUE}.header
#cat ${TMP}/1000GIBS.header | cut -f8 # GT for now, but all NAs till redo
#fi


#real/tmp/2TISSUES.M116.header


# full table of just rare variants in homozygosity
cat ${TMP}/${TISSUED}.${IDD}.header > ${TMP}/${TISSUED}.${IDD}.rare.homozygous
cat ${INPUTS} | awk '{ if ( $9=="1.0" || $9=="1" ) print $0 }' >> ${TMP}/${TISSUED}.${IDD}.rare.homozygous


cat  ${TMP}/${TISSUED}.${IDD}.rare.homozygous | grep -v -w ID | cut -f3 | sort -u > ${OUTPUTS}/${TISSUED}.${IDD}.rare.homozygous.variants
#cat  ${TMP}/${TISSUED}.${IDD}.rare.homozygous | grep -v -w ID | cut -f16 | sort -u > ${OUTPUTS}/${TISSUED}.${IDD}.rare.homozygous.genes
cat  ${TMP}/${TISSUED}.${IDD}.rare.homozygous | grep -v -w ID | grep -f ${OUTPUTS}/${TISSUED}.${IDD}.rare.homozygous.variants | cut -f16 | sort -u > ${OUTPUTS}/${TISSUED}.${IDD}.rare.homozygous.genes


# Filter by CADD>15 
#Rscript ./real/scripts/cadd_fil.R ${TMP}/${IDD}.rare.homozygous ${OUTPUTS} rare


# per category -------------------------

for CATEGORY in ${CATEGORIES[@]}

do

cat ${TMP}/${TISSUED}.${IDD}.header > ${TMP}/${TISSUED}.${IDD}.${CATEGORY}.homozygous
cat ${TMP}/${TISSUED}.${IDD}.rare.homozygous | grep -w ${CATEGORY}  >> ${TMP}/${TISSUED}.${IDD}.${CATEGORY}.homozygous

cat  ${TMP}/${TISSUED}.${IDD}.${CATEGORY}.homozygous | grep -v -w ID | cut -f3 | sort -u > ${OUTPUTS}/${TISSUED}.${IDD}.${CATEGORY}.homozygous.variants

#cat  ${TMP}/${TISSUED}.${IDD}.${CATEGORY}.homozygous | grep -v -w ID | cut -f16 | sort -u > ${OUTPUTS}/${TISSUED}.${IDD}.${CATEGORY}.homozygous.genes


cat  ${TMP}/${TISSUED}.${IDD}.${CATEGORY}.homozygous | grep -v -w ID | grep -f ${OUTPUTS}/${TISSUED}.${IDD}.${CATEGORY}.homozygous.variants | cut -f16 | sort -u  > ${OUTPUTS}/${TISSUED}.${IDD}.${CATEGORY}.homozygous.genes


# not seen in gnomad v4
cat ${OUTPURS}/${TISSUED}.${IDD}.${CATEGORY}.homozygous.variants  | grep -v -f real/inputs/gnomad.homozygous.chr* > ${OUTPUTS}/${TISSUED}.${IDD}.${CATEGORY}.homozygous.variants.novel


# EXTRA CAT
# CADD >15
#cat real/outputs/blood/M116/blood.M116.*.eur_rare_all.fulltable | awk '{ if ( ( $90+0>15)  && ($9=="1.0" || $9=="1") ) print $0 }' | cut -f3 | sort -u 
#no, do it in R, doesnt quite work with 2 digits and stuff
# Filter by CADD > 15, and VOI categories to get VOI categories plus NCNC CADD>15 so far not considered
#Rscript ./real/scripts/cadd_fil.R ${TMP}/${IDD}.${CATEGORY}.homozygous ${OUTPUTS} ${CATEGORY}
# Now he have a category for this

done
