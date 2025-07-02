#!/bin/bash
#
#SBATCH -p normal # partition (queue)
#SBATCH -N 1 # number of nodes
#SBATCH -J compbg
#SBATCH --mem 32G # memory pool for all cores
#SBATCH -t 0-03:59# time (D-HH:MM)
#SBATCH -o ./log.%j.out # STDOUT
#SBATCH -e ./log.%j.err # STDERR
#SBATCH --mail-type=END,FAIL # notifications for job done & fail
#SBATCH --mail-user=claudia.vasallo@upf.edu # send-to address




TISSUED=$1
IDD=$2
INPUTS=$3
OUTPUTS=$4
CAT_BG_DIR=$5



COD_CATEGORIES=( damaging moderate altering modifier_cd cadd15 low_only )


for CATEGORY in ${COD_CATEGORIES[@]}

do

REFERENCE=coding
cat ${INPUTS}/${TISSUED}.${IDD}.eur_rare_${REFERENCE}.autosomal.genes | grep -f ${CAT_BG_DIR}/everyones.${CATEGORY}.genes  > ${OUTPUTS}/${TISSUED}.${IDD}.${CATEGORY}.compound.bg

done



COD_EXTRA_CATEGORIES=( modifier_nc )
for CATEGORY in ${COD_EXTRA_CATEGORIES[@]}
do
REFERENCE=coding_with_nc
cat ${INPUTS}/${TISSUED}.${IDD}.eur_rare_${REFERENCE}.autosomal.genes | grep -f ${CAT_BG_DIR}/everyones.${CATEGORY}.genes  > ${OUTPUTS}/${TISSUED}.${IDD}.${CATEGORY}.compound.bg
done



NC_CATEGORIES=( ncnc_cadd15 )
for CATEGORY in ${NC_CATEGORIES[@]}
do
REFERENCE=non_coding
cat ${INPUTS}/${TISSUED}.${IDD}.eur_rare_${REFERENCE}.autosomal.genes | grep -f ${CAT_BG_DIR}/everyones.${CATEGORY}.genes  > ${OUTPUTS}/${TISSUED}.${IDD}.${CATEGORY}.compound.bg
done

