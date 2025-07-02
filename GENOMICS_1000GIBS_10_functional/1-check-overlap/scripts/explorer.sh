#!/bin/bash
#
#SBATCH -p normal # partition (queue)
#SBATCH -N 1 # number of nodes
#SBATCH -J explorer
#SBATCH --mem 16G # memory pool for all cores
#SBATCH -t 0-00:29# time (D-HH:MM)
#SBATCH -o ./real/outputs2/log.%j.out # STDOUT
#SBATCH -e ./real/outputs2/log.%j.err # STDERR
#SBATCH --mail-type=END,FAIL # notifications for job done & fail
#SBATCH --mail-user=claudia.vasallo@upf.edu # send-to address



INPUTS=$1    
#  ../../09-summary-in-shared-by-samples/real/summary_tables
OUTPUTS=$2
TISSUED=$3    #2TISSUES
ARRAY_LIST=$4         #M116
GENESETS_DIR=$5  
GENESET=$6  
CATEGORY_FILE=$7 
OVERLAP_DIR=$8     # real/outputs


CAT_BG_DIR=/gpfs42/projects/lab_anavarro/disease_pleiotropies/supercentenarian/Olot/10-functional-in-shared-by-sample/x-goof_red_genelists/real/outputs


#ARRAY_LIST=$CAT_FILE

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








CATEGORIES=$( cat ${CATEGORY_FILE} )

cat ${OUTPUTS}/overlaps.table.header  > ${OUTPUTS}/${TISSUED}.${IDD}.${GENESET}.overlaps.table


for CATEGORY in ${CATEGORIES[@]}

do

all="all"
if [[ "$CATEGORY" == "$all" ]]
then
CATEGORYWORD="eur_rare"
else
CATEGORYWORD="eur_rare_${CATEGORY}"
fi

GENESET_SIZE=$(cat ${GENESETS_DIR}/${GENESET}_ensg.txt | wc -l )

#CAT_TOT=$(cat ${INPUTS}/${TISSUED}.${IDD}.eur_rare_${CATEGORY}.autosomal.genes | wc -l)

CAT_TOT=$(cat ${INPUTS}/${TISSUED}.${IDD}.${CATEGORYWORD}.genes | wc -l)


# removed autosomal from file name that is only in M116


#CAT_IN=$(cat ${OVERLAP_DIR}/${TISSUED}.${IDD}.eur_rare_${CATEGORY}.${GENESET} | wc -l)

CAT_IN=$(cat ${OVERLAP_DIR}/${TISSUED}.${IDD}.${CATEGORYWORD}.${GENESET} | wc -l)


#CAT_OUT=$(cat ${INPUTS}/${TISSUED}.${IDD}.eur_rare_${CATEGORY}.autosomal.genes| grep -v -f ${OVERLAP_DIR}/${TISSUED}.${IDD}.eur_rare_${CATEGORY}.${GENESET} | wc -l)

CAT_OUT=$(cat ${INPUTS}/${TISSUED}.${IDD}.${CATEGORYWORD}.genes | grep -v -f ${OVERLAP_DIR}/${TISSUED}.${IDD}.${CATEGORYWORD}.${GENESET} | wc -l)

echo ${GENESET} ${GENESET_SIZE} ${CATEGORY} ${CAT_TOT} ${CAT_IN} ${CAT_OUT} >> ${OUTPUTS}/${TISSUED}.${IDD}.${GENESET}.overlaps.table


done






# CAT REFERENCE CODING

REFERENCE=coding   # added BG of altering
CATEGORY=altering
CAT_TOT=$(cat ${INPUTS}/${TISSUED}.${IDD}.eur_rare_${REFERENCE}.genes | grep -f ${CAT_BG_DIR}/everyones.${CATEGORY}.genes  | wc -l ) 
CAT_IN=$(cat ${INPUTS}/${TISSUED}.${IDD}.eur_rare_${REFERENCE}.genes | grep -f ${CAT_BG_DIR}/everyones.${CATEGORY}.genes  | grep -f ${OVERLAP_DIR}/${TISSUED}.${IDD}.eur_rare_${REFERENCE}.${GENESET} |wc -l)
CAT_OUT=$( cat ${INPUTS}/${TISSUED}.${IDD}.eur_rare_${REFERENCE}.genes | grep -f ${CAT_BG_DIR}/everyones.${CATEGORY}.genes | grep -v -f ${OVERLAP_DIR}/${TISSUED}.${IDD}.eur_rare_${REFERENCE}.${GENESET} | wc -l)
echo ${GENESET} ${GENESET_SIZE} "COD_ALTERING_BG" ${CAT_TOT} ${CAT_IN} ${CAT_OUT} >> ${OUTPUTS}/${TISSUED}.${IDD}.${GENESET}.overlaps.table


REFERENCE=coding
CATEGORY=damaging

CAT_TOT=$(cat ${INPUTS}/${TISSUED}.${IDD}.eur_rare_${REFERENCE}.genes | grep -f ${CAT_BG_DIR}/everyones.${CATEGORY}.genes  | wc -l ) 
CAT_IN=$(cat ${INPUTS}/${TISSUED}.${IDD}.eur_rare_${REFERENCE}.genes | grep -f ${CAT_BG_DIR}/everyones.${CATEGORY}.genes  | grep -f ${OVERLAP_DIR}/${TISSUED}.${IDD}.eur_rare_${REFERENCE}.${GENESET} |wc -l)
CAT_OUT=$( cat ${INPUTS}/${TISSUED}.${IDD}.eur_rare_${REFERENCE}.genes | grep -f ${CAT_BG_DIR}/everyones.${CATEGORY}.genes | grep -v -f ${OVERLAP_DIR}/${TISSUED}.${IDD}.eur_rare_${REFERENCE}.${GENESET} | wc -l)
echo ${GENESET} ${GENESET_SIZE} "COD_DAMAGING_BG" ${CAT_TOT} ${CAT_IN} ${CAT_OUT} >> ${OUTPUTS}/${TISSUED}.${IDD}.${GENESET}.overlaps.table






REFERENCE=coding
CATEGORY=moderate

CAT_TOT=$(cat ${INPUTS}/${TISSUED}.${IDD}.eur_rare_${REFERENCE}.genes | grep -f ${CAT_BG_DIR}/everyones.${CATEGORY}.genes  | wc -l ) 
CAT_IN=$(cat ${INPUTS}/${TISSUED}.${IDD}.eur_rare_${REFERENCE}.genes | grep -f ${CAT_BG_DIR}/everyones.${CATEGORY}.genes  | grep -f ${OVERLAP_DIR}/${TISSUED}.${IDD}.eur_rare_${REFERENCE}.${GENESET} |wc -l)
CAT_OUT=$( cat ${INPUTS}/${TISSUED}.${IDD}.eur_rare_${REFERENCE}.genes | grep -f ${CAT_BG_DIR}/everyones.${CATEGORY}.genes | grep -v -f ${OVERLAP_DIR}/${TISSUED}.${IDD}.eur_rare_${REFERENCE}.${GENESET} | wc -l)
echo ${GENESET} ${GENESET_SIZE} "COD_MODERATE_BG" ${CAT_TOT} ${CAT_IN} ${CAT_OUT} >> ${OUTPUTS}/${TISSUED}.${IDD}.${GENESET}.overlaps.table


REFERENCE=coding
CATEGORY=modifier_cd

CAT_TOT=$(cat ${INPUTS}/${TISSUED}.${IDD}.eur_rare_${REFERENCE}.genes | grep -f ${CAT_BG_DIR}/everyones.${CATEGORY}.genes  | wc -l ) 
CAT_IN=$(cat ${INPUTS}/${TISSUED}.${IDD}.eur_rare_${REFERENCE}.genes | grep -f ${CAT_BG_DIR}/everyones.${CATEGORY}.genes  | grep -f ${OVERLAP_DIR}/${TISSUED}.${IDD}.eur_rare_${REFERENCE}.${GENESET} |wc -l)
CAT_OUT=$( cat ${INPUTS}/${TISSUED}.${IDD}.eur_rare_${REFERENCE}.genes | grep -f ${CAT_BG_DIR}/everyones.${CATEGORY}.genes | grep -v -f ${OVERLAP_DIR}/${TISSUED}.${IDD}.eur_rare_${REFERENCE}.${GENESET} | wc -l)
echo ${GENESET} ${GENESET_SIZE} "COD_MODIFIER_CD_BG" ${CAT_TOT} ${CAT_IN} ${CAT_OUT} >> ${OUTPUTS}/${TISSUED}.${IDD}.${GENESET}.overlaps.table



REFERENCE=coding
CATEGORY=modifier_nc

CAT_TOT=$(cat ${INPUTS}/${TISSUED}.${IDD}.eur_rare_${REFERENCE}.genes | grep -f ${CAT_BG_DIR}/everyones.${CATEGORY}.genes  | wc -l ) 
CAT_IN=$(cat ${INPUTS}/${TISSUED}.${IDD}.eur_rare_${REFERENCE}.genes | grep -f ${CAT_BG_DIR}/everyones.${CATEGORY}.genes  | grep -f ${OVERLAP_DIR}/${TISSUED}.${IDD}.eur_rare_${REFERENCE}.${GENESET} |wc -l)
CAT_OUT=$( cat ${INPUTS}/${TISSUED}.${IDD}.eur_rare_${REFERENCE}.genes | grep -f ${CAT_BG_DIR}/everyones.${CATEGORY}.genes | grep -v -f ${OVERLAP_DIR}/${TISSUED}.${IDD}.eur_rare_${REFERENCE}.${GENESET} | wc -l)
echo ${GENESET} ${GENESET_SIZE} "COD_MODIFIER_NC_BG" ${CAT_TOT} ${CAT_IN} ${CAT_OUT} >> ${OUTPUTS}/${TISSUED}.${IDD}.${GENESET}.overlaps.table



REFERENCE=non_coding
CATEGORY=ncnc_cadd15

CAT_TOT=$(cat ${INPUTS}/${TISSUED}.${IDD}.eur_rare_${REFERENCE}.genes | grep -f ${CAT_BG_DIR}/everyones.${CATEGORY}.genes  | wc -l ) 
CAT_IN=$(cat ${INPUTS}/${TISSUED}.${IDD}.eur_rare_${REFERENCE}.genes | grep -f ${CAT_BG_DIR}/everyones.${CATEGORY}.genes  | grep -f ${OVERLAP_DIR}/${TISSUED}.${IDD}.eur_rare_${REFERENCE}.${GENESET} |wc -l)
CAT_OUT=$( cat ${INPUTS}/${TISSUED}.${IDD}.eur_rare_${REFERENCE}.genes | grep -f ${CAT_BG_DIR}/everyones.${CATEGORY}.genes | grep -v -f ${OVERLAP_DIR}/${TISSUED}.${IDD}.eur_rare_${REFERENCE}.${GENESET} | wc -l)
echo ${GENESET} ${GENESET_SIZE} "NC_NCNC_CADD15_BG" ${CAT_TOT} ${CAT_IN} ${CAT_OUT} >> ${OUTPUTS}/${TISSUED}.${IDD}.${GENESET}.overlaps.table


REFERENCE=coding
CATEGORY=low_only

CAT_TOT=$(cat ${INPUTS}/${TISSUED}.${IDD}.eur_rare_${REFERENCE}.genes | grep -f ${CAT_BG_DIR}/everyones.${CATEGORY}.genes  | wc -l ) 
CAT_IN=$(cat ${INPUTS}/${TISSUED}.${IDD}.eur_rare_${REFERENCE}.genes | grep -f ${CAT_BG_DIR}/everyones.${CATEGORY}.genes  | grep -f ${OVERLAP_DIR}/${TISSUED}.${IDD}.eur_rare_${REFERENCE}.${GENESET} |wc -l)
CAT_OUT=$( cat ${INPUTS}/${TISSUED}.${IDD}.eur_rare_${REFERENCE}.genes | grep -f ${CAT_BG_DIR}/everyones.${CATEGORY}.genes | grep -v -f ${OVERLAP_DIR}/${TISSUED}.${IDD}.eur_rare_${REFERENCE}.${GENESET} | wc -l)
echo ${GENESET} ${GENESET_SIZE} "COD_LOW_ONLY_BG" ${CAT_TOT} ${CAT_IN} ${CAT_OUT} >> ${OUTPUTS}/${TISSUED}.${IDD}.${GENESET}.overlaps.table



REFERENCE=coding_with_nc
CATEGORY=modifier_nc

CAT_TOT=$(cat ${INPUTS}/${TISSUED}.${IDD}.eur_rare_${REFERENCE}.genes | grep -f ${CAT_BG_DIR}/everyones.${CATEGORY}.genes  | wc -l ) 
CAT_IN=$(cat ${INPUTS}/${TISSUED}.${IDD}.eur_rare_${REFERENCE}.genes | grep -f ${CAT_BG_DIR}/everyones.${CATEGORY}.genes  | grep -f ${OVERLAP_DIR}/${TISSUED}.${IDD}.eur_rare_${REFERENCE}.${GENESET} | wc -l)
CAT_OUT=$( cat ${INPUTS}/${TISSUED}.${IDD}.eur_rare_${REFERENCE}.genes | grep -f ${CAT_BG_DIR}/everyones.${CATEGORY}.genes | grep -v -f ${OVERLAP_DIR}/${TISSUED}.${IDD}.eur_rare_${REFERENCE}.${GENESET} | wc -l)
echo ${GENESET} ${GENESET_SIZE} "CDWNC_MODIFIER_NC_BG" ${CAT_TOT} ${CAT_IN} ${CAT_OUT} >> ${OUTPUTS}/${TISSUED}.${IDD}.${GENESET}.overlaps.table


# all as reference for cadd15, but anyway annotations are in coding
REFERENCE=all
CATEGORY=cadd15

CAT_TOT=$(cat ${INPUTS}/${TISSUED}.${IDD}.eur_rare.genes | grep -f ${CAT_BG_DIR}/everyones.${CATEGORY}.genes  | wc -l )
CAT_IN=$(cat ${INPUTS}/${TISSUED}.${IDD}.eur_rare.genes | grep -f ${CAT_BG_DIR}/everyones.${CATEGORY}.genes  | grep -f ${OVERLAP_DIR}/${TISSUED}.${IDD}.eur_rare.${GENESET} | wc -l)
CAT_OUT=$( cat ${INPUTS}/${TISSUED}.${IDD}.eur_rare.genes | grep -f ${CAT_BG_DIR}/everyones.${CATEGORY}.genes | grep -v -f ${OVERLAP_DIR}/${TISSUED}.${IDD}.eur_rare.${GENESET} | wc -l)
echo ${GENESET} ${GENESET_SIZE} "ALL_CADD15_BG" ${CAT_TOT} ${CAT_IN} ${CAT_OUT} >> ${OUTPUTS}/${TISSUED}.${IDD}.${GENESET}.overlaps.table



REFERENCE=coding
CATEGORY=cadd15

CAT_TOT=$(cat ${INPUTS}/${TISSUED}.${IDD}.eur_rare_${REFERENCE}.genes | grep -f ${CAT_BG_DIR}/everyones.${CATEGORY}.genes  | wc -l )
CAT_IN=$(cat ${INPUTS}/${TISSUED}.${IDD}.eur_rare_${REFERENCE}.genes | grep -f ${CAT_BG_DIR}/everyones.${CATEGORY}.genes  | grep -f ${OVERLAP_DIR}/${TISSUED}.${IDD}.eur_rare_${REFERENCE}.${GENESET} | wc -l)
CAT_OUT=$( cat ${INPUTS}/${TISSUED}.${IDD}.eur_rare_${REFERENCE}.genes | grep -f ${CAT_BG_DIR}/everyones.${CATEGORY}.genes | grep -v -f ${OVERLAP_DIR}/${TISSUED}.${IDD}.eur_rare_${REFERENCE}.${GENESET} | wc -l)
echo ${GENESET} ${GENESET_SIZE} "COD_CADD15_BG" ${CAT_TOT} ${CAT_IN} ${CAT_OUT} >> ${OUTPUTS}/${TISSUED}.${IDD}.${GENESET}.overlaps.table

