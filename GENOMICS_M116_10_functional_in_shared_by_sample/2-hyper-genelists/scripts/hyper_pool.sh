#!/bin/bash
#
#SBATCH -p normal # partition (queue)
#SBATCH -N 1 # number of nodes
#SBATCH -J hyperpool
#SBATCH --mem 16G # memory pool for all cores
#SBATCH -t 0-00:29# time (D-HH:MM)
#SBATCH -o ./real/outputs/log.%j.out # STDOUT
#SBATCH -e ./real/outputs/log.%j.err # STDERR
#SBATCH --mail-type=END,FAIL # notifications for job done & fail
#SBATCH --mail-user=claudia.vasallo@upf.edu # send-to address


# For all chromosomes together
TABLE=$1       
INPUTS=$2      
OUTPUTS=$3      
GENELISTS_FOLDER=$4
CAT_BG_FOLDER=$5     
TISSUED=$6
ARRAY_LIST=$7   # IDD
CAT_FILE=$8
TYPE=$8



# Array Dependent Config
#If there is no array number assigned:
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




CATEGORIES=$( cat ${CAT_FILE} )


# Results per category will go here:

#echo GENELIST RAREGENES LONGEVITY_GENES GENELIST LONGEVITY_IN_GENELIST PVAL OR > ${OUTPUTS}/hyper.results.txt
echo CATEGORY BG LONGEVITY_IN_BG CAT LONGEVITY_IN_CAT PVAL OR > ${OUTPUTS}/hyper.results.txt

# All longevity lists together (Manel gene lists) -------------------------------------------------------------------------------------------


#TISSUED=TWO
#IDD=M116

#if [ "$IDD" == "M116" ]
#then
#LEVEL=${TABLE}/*/M116/*.M116.*.eur_rare_all.fulltable
#else
#LEVEL=${TABLE}/blood/M116/blood.M116.*.eur_rare_all.fulltable
#fi


# then BG will be 10866 genes for M116 instead of the 10786 when seleting blood cropped dataset only
# BG per category

#GENELISTS_FOLDER=$1

#CATEGORIES=( damaging moderate modifier_cd modifier_nc altering low_only )


#---------------------------------
cat ${GENELISTS_FOLDER}/*.txt | sort -u > ./real/inputs/all_longv_lists.txt
ALL_LONGV_GENELISTS=./real/inputs/all_longv_lists.txt


# POPULATION (All Genes in M116 sample in LEVEL: here shared in 2 or 3 samples)
# autosomal genes only to compare with 1000GIBS
#cat ${LEVEL} | grep -v -w chrX | grep -v -w chrY | grep -v -w chrM | grep Transcript | grep protein_coding | grep -v -w ID | cut -f16 | sort -u | grep -v -w "" > ./real/inputs/BG_${TISSUED}_${IDD}.genes
#cat ../../09-summary-in-shared-by-samples/real/outputs/*/M116/*.M116.*.eur_rare_all.fulltable | grep -v -w chrX | grep -v -w chrY | grep -v -w chrM | grep Transcript  | grep CODING-GENE | grep -v -w ID | cut -f16 | sort -u | grep -v -w "" > ./real/inputs/BG_M116.coding.genes
#cat ../../../1000GIBS/08-filter-category-summary/real/outputs/1000GBIB/*/*.eur_rare_all.fulltable | grep -v -w chrX | grep -v -w chrY | grep -v -w chrM | grep Transcript  | grep CODING-GENE | grep -v -w ID | cut -f16 | sort -u | grep -v -w "" > ./real/inputs/BG_1000GIBS.coding.genes

#CODING_GENES=./real/inputs/BG_everyones.coding.genes
CODING_GENES=Homo_sapiens.GRCh38.104.protein_coding.txt


# SUCCESES IN POPULATION (Longevity genes in ALLGENES in M116 sample)
#cat ${ALL_LONGV_GENELISTS} | grep -f  ${BG_GENES} | sort -u > ./real/inputs/LONGEVITY_IN_BG_${TISSUED}_${IDD}.genes
#LONGEVITY_IN_BG_GENES=./real/inputs/LONGEVITY_IN_BG_${TISSUED}_${IDD}.genes
#LONGEVITY_IN_BG=$(cat ${LONGEVITY_IN_BG_GENES}  | grep -v -i -w ensg | cut -f2 | sort -u | wc -l)



for CATEGORY in ${CATEGORIES[@]}
do

#CATEGORY=otra
#all="all"
#if [[ "$CATEGORY" == "$all" ]]
#then
#	#CATEGORY=""
#	echo "all"
#else
#	CATEGORY="_"$CATEGORY
#fi
#echo $CATEGORY


CODE="CODING"
if [[ "$TYPE" == "$CODE" ]]
then
BG_GENES_ORI=${CAT_BG_FOLDER}/everyones.${CATEGORY}.genes
cat ${BG_GENES_ORI} | grep -f ${CODING_GENES} > ./real/inputs/everyones.${CATEGORY}.coding.genes
BG_GENES=./real/inputs/everyones.${CATEGORY}.coding.genes
else
BG_GENES=${CAT_BG_FOLDER}/everyones.${CATEGORY}.genes
fi


BG=$(cat ${BG_GENES} | wc -l)



# SUCCESES IN POPULATION (Longevity genes in ALLGENES in M116 sample)
cat ${ALL_LONGV_GENELISTS} | grep -f  ${BG_GENES} | sort -u > ./real/inputs/LONGEVITY_IN_BG_${TISSUED}_${IDD}.genes
LONGEVITY_IN_BG_GENES=./real/inputs/LONGEVITY_IN_BG_${TISSUED}_${IDD}.genes
LONGEVITY_IN_BG=$(cat ${LONGEVITY_IN_BG_GENES}  | grep -v -i -w ensg | cut -f2 | sort -u | wc -l)



#---------
all="all"

if [[ "$CATEGORY" == "$all" ]]
then
	CATEGORY_GENES_ORI=${INPUTS}/eur_rare.genes
        echo "all"
else
    	CATEGORY_GENES_ORI=${INPUTS}/eur_rare_${CATEGORY}.genes
fi



#CATEGORY_GENES_ORI=${INPUTS}/eur_rare${CATEGORY}.genes
cat ${CATEGORY_GENES_ORI} | grep -f ${BG_GENES} > ./real/inputs/CAT_${CATEGORY}.${TISSUED}.${IDD}.genes
CATEGORY_GENES=./real/inputs/CAT_${CATEGORY}.${TISSUED}.${IDD}.genes
CAT=$(cat ${CATEGORY_GENES} | sort -u | wc -l)

# SUCCESES IN SAMPLE (Longevity genes in Genes with rare and altering/damaging in M116 sample)
#LONGEVITY_IN_LIST=$(cat ../1-check-overlap/real/outputs/blood/M116/blood.M116.*.all.any.genes  | wc -l)
#LONGEVITY_IN_LIST=$( cat ${CATEGORY_GENES} | grep -f ${ALL_LONGV_GENELISTS} | wc -l)
# more explicitly > with new explicite variable LONGEVITY_IN_BG as well

cat ${CATEGORY_GENES} | grep -f ${LONGEVITY_IN_BG_GENES} | sort -u  > ./real/inputs/LONGEVITY_IN_CAT_${CATEGORY}.${TISSUED}_${IDD}.genes
LONGEVITY_IN_CAT_GENES=./real/inputs/LONGEVITY_IN_CAT_${CATEGORY}.${TISSUED}_${IDD}.genes 
LONGEVITY_IN_CAT=$( cat ${LONGEVITY_IN_CAT_GENES} | wc -l ) 

module load R

HYPER=$(Rscript ./real/scripts/hyper.R ${BG} ${LONGEVITY_IN_BG} ${CAT} ${LONGEVITY_IN_CAT})
# 4.615873e-34

#echo $GENELIST $(echo $HYPER | cut -d ' ' -f2) >> ${OUTPUTS}/hyper.results.txt

# now pvalue and odds ratio
PVAL=$(echo $HYPER | cut -d ' ' -f2)
OR=$(echo $HYPER | cut -d ' ' -f3)

echo ${CATEGORY} ${BG} ${LONGEVITY_IN_BG} ${CAT} ${LONGEVITY_IN_CAT} $PVAL $OR >> ${OUTPUTS}/hyper.results.txt



done

#phyper(q, m, n, k, lower.tail = TRUE, log.p = FALSE)
#x, q vector of quantiles representing the number of white balls drawn
#without replacement from an urn which contains both black and white
#balls.
#m the number of white balls in the urn.
#n the number of black balls in the urn.
#k the number of balls drawn from the urn.

