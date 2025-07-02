#!/bin/bash
#
#SBATCH -p normal # partition (queue)
#SBATCH -N 1 # number of nodes
#SBATCH -J overlong
#SBATCH --mem 64G # memory pool for all cores
#SBATCH -t 9-23:59# time (D-HH:MM)
#SBATCH -o ./real/outputs/log.%j.out # STDOUT
#SBATCH -e ./real/outputs/log.%j.err # STDERR
#SBATCH --mail-type=END,FAIL # notifications for job done & fail
#SBATCH --mail-user=claudia.vasallo@upf.edu # send-to address



# Modules

#module load R
module load R/3.6.0-foss-2018b

# Config
ARRAY_LIST=$1
INPUTS=$2
OUTPUTS=$3
TISSUED=$4
IDD=$5



SUFFIX=eur_rare_all.fulltable

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


# Genelists
GENE_LISTS=($(ls ../0-gene-lists/real/outputs/ | grep txt | grep -v all | cut -d '.' -f1))



# ALL GENES WITH RARE DAMAGING/ALTERING/MODIFIER VARIANTS IN CHR 
# GENES WITH VOI, ANY



CATEGORIES=( $(cat ./real/inputs/cats.txt) )


SC="M116"

for CAT in ${CATEGORIES[@]}

do

if [[ "$IDD" ==  "$SC" ]]
then

# Problem with MODERATE which is an IMPACT value and category name: FIXME
#cat ${INPUTS_DIR}/*/${IDD}/*.${IDD}.${CHR}.${SUFFIX} | grep -w ${CAT} | cut -f16  | sort -u > ${OUTPUTS_DIR}/${TISSUED}/${IDD}/${TISSUED}.${IDD}.${CHR}.${CAT}.sample
# USE genevar instead
#cat $(ls ${INPUTS}/*/${IDD}/*.${IDD}.${CHR}.* |  grep genevar | grep -i  "rare_""${CAT}" ) | cut -f16  | sort -u 

else


#cat $(ls ${INPUTS}/${TISSUED}/${IDD}/${TISSUED}.${IDD}.${CHR}.* 

fi

#cat ${INPUTS_DIR}/${TISSUED}/${IDD}/${TISSUED}.${IDD}.${CHR}.${SUFFIX} | grep MODERATE | cut -f16  | sort -u > ${OUTPUTS_DIR}/${TISSUED}/${IDD}/${TISSUED}.${ID}.${CHR}.damaging.sample
#cat ${INPUTS_DIR}/${TISSUED}/${IDD}/${TISSUED}.${IDD}.${CHR}.${SUFFIX} | grep MODIFIER_CD | cut -f16  | sort -u > ${OUTPUTS_DIR}/${TISSUED}/${IDD}/${TISSUED}.${ID}.${CHR}.moderate.sample
#cat ${INPUTS_DIR}/${TISSUED}/${IDD}/${TISSUED}.${IDD}.${CHR}.${SUFFIX} | grep "ALTERING\|DAMAGING\|MODERATE" | cut -f15  | sort -u > ${OUTPUTS_DIR}/${TISSUED}/${IDD}/${TISSUED}.${IDD}.${CHR}.any.sample

done




for GENE_LIST in ${GENE_LISTS[@]}
do

echo ${GENE_LIST}


for CAT in ${CATEGORIES[@]}

do
cat ${INPUTS_DIR}/${TISSUED}/${IDD}/${TISSUED}.${IDD}.${CHR}.${SUFFIX} | head -n1 > ${OUTPUTS_DIR}/${TISSUED}/${IDD}/${TISSUED}.${IDD}.${CHR}.${GENE_LIST}.${CAT}.hits
cat ${INPUTS_DIR}/${TISSUED}/${IDD}/${TISSUED}.${IDD}.${CHR}.${SUFFIX} | grep ${CAT} | grep -f ../0-gene-lists/real/outputs/${GENE_LIST}.txt >> ${OUTPUTS_DIR}/${TISSUED}/${IDD}/${TISSUED}.${IDD}.${CHR}.${GENE_LIST}.${CAT}.hits
done

done



# LONGEVITY GENES WITH RARE DAMAGING/ALTERING VARIANTS IN CHR
# different positions per chr, genes


for GENE_LIST in ${GENE_LISTS[@]}

do
echo ${GENE_LIST}

for CAT in ${CATEGORIES[@]}

do
echo ${GENE_LIST} $(cat ${OUTPUTS_DIR}/${TISSUED}/${IDD}/${TISSUED}.${IDD}.${CHR}.${GENE_LIST}.${CAT}.hits |  grep -v SYMBOL | cut -f2  | sort -u | wc -l) $(cat ${OUTPUTS_DIR}/${TISSUED}/${IDD}/${TISSUED}.${IDD}.${CHR}.${GENE_LIST}.${CAT}.hits  | cut -f12  | sort -u | wc -l) > ${OUTPUTS_DIR}/${TISSUED}/${IDD}/${TISSUED}.${IDD}.${CHR}.${CAT}.variants_genes

cat ${OUTPUTS_DIR}/${TISSUED}/${IDD}/${TISSUED}.${IDD}.${CHR}.${GENE_LIST}.${CAT}.hits | grep -v SYMBOL | cut -f16  | sort -u > ${OUTPUTS_DIR}/${TISSUED}/${IDD}/${TISSUED}.${IDD}.${CHR}.${GENE_LIST}.${CAT}.genes

done


# DAMAGING
#echo ${GENE_LIST} $(cat ${OUTPUTS_DIR}/${TISSUED}/${IDD}/${TISSUED}.${IDD}.${CHR}.${GENE_LIST}.damaging.hits |  grep -v SYMBOL | cut -f2  | sort -u | wc -l) $(cat ${OUTPUTS_DIR}/${TISSUED}/${IDD}/${TISSUED}.${IDD}.${CHR}.${GENE_LIST}.hits  | cut -f12  | sort -u | wc -l) > ${OUTPUTS_DIR}/${TISSUED}/${IDD}/${TISSUED}.${IDD}.${CHR}.damaging.source_variants_genes
#cat ${OUTPUTS_DIR}/${TISSUED}/${IDD}/${TISSUED}.${IDD}.${CHR}.${GENE_LIST}.damaging.hits | grep -v SYMBOL | cut -f14  | sort -u > ${OUTPUTS_DIR}/${TISSUED}/${IDD}/${TISSUED}.${IDD}.${CHR}.${GENE_LIST}.damaging.genes


# ALTERING
#echo ${GENE_LIST} $(cat ${OUTPUTS_DIR}/${TISSUED}/${IDD}/${TISSUED}.${IDD}.${CHR}.${GENE_LIST}.damaging.hits |  grep -v SYMBOL | cut -f2  | sort -u | wc -l) $(cat ${OUTPUTS_DIR}/${TISSUED}/${IDD}/${TISSUED}.${IDD}.${CHR}.${GENE_LIST}.altering.hits  | cut -f12  | sort -u | wc -l) > ${OUTPUTS_DIR}/${TISSUED}/${IDD}/${TISSUED}.${IDD}.${CHR}.altering.source_variants_genes
#cat ${OUTPUTS_DIR}/${TISSUED}/${IDD}/${TISSUED}.${IDD}.${CHR}.${GENE_LIST}.altering.hits | grep -v SYMBOL | cut -f14  | sort -u > ${OUTPUTS_DIR}/${TISSUED}/${IDD}/${TISSUED}.${IDD}.${CHR}.${GENE_LIST}.altering.genes


# MODERATE
#echo ${GENE_LIST} $(cat ${OUTPUTS_DIR}/${TISSUED}/${IDD}/${TISSUED}.${IDD}.${CHR}.${GENE_LIST}.moderate.hits |  grep -v SYMBOL | cut -f2  | sort -u | wc -l) $(cat ${OUTPUTS_DIR}/${TISSUED}/${IDD}/${TISSUED}.${IDD}.${CHR}.${GENE_LIST}.moderate.hits  | cut -f12  | sort -u | wc -l) > ${OUTPUTS_DIR}/${TISSUED}/${IDD}/${TISSUED}.${IDD}.${CHR}.moderate.source_variants_genes
#cat ${OUTPUTS_DIR}/${TISSUED}/${IDD}/${TISSUED}.${IDD}.${CHR}.${GENE_LIST}.moderate.hits | grep -v SYMBOL | cut -f14  | sort -u > ${OUTPUTS_DIR}/${TISSUED}/${IDD}/${TISSUED}.${IDD}.${CHR}.${GENE_LIST}.moderate.genes



# or -f15 for ENSG and then id_convert
#Rscript ./real/scripts/id_convert.R ${OUTPUTS_DIR}/${TISSUED}/${IDD}/${TISSUED}.${IDD}.${CHR}.${GENE_LIST}.genes


done




for CAT in ${CATEGORIES[@]}

do


cat ${OUTPUTS_DIR}/${TISSUED}/${IDD}/${TISSUED}.${IDD}.${CHR}.*.${CAT}.genes  | sort -u > ${OUTPUTS_DIR}/${TISSUED}/${IDD}/${TISSUED}.${IDD}.${CHR}.all.${CAT}.genes
#cat ${OUTPUTS_DIR}/${TISSUED}/${IDD}/${TISSUED}.${IDD}.${CHR}.*.damaging.genes  | sort -u > ${OUTPUTS_DIR}/${TISSUED}/${IDD}/${TISSUED}.${IDD}.${CHR}.all.damaging.genes
#cat ${OUTPUTS_DIR}/${TISSUED}/${IDD}/${TISSUED}.${IDD}.${CHR}.*.altering.genes  | sort -u > ${OUTPUTS_DIR}/${TISSUED}/${IDD}/${TISSUED}.${IDD}.${CHR}.all.altering.genes
#cat ${OUTPUTS_DIR}/${TISSUED}/${IDD}/${TISSUED}.${IDD}.${CHR}.*.moderate.genes  | sort -u > ${OUTPUTS_DIR}/${TISSUED}/${IDD}/${TISSUED}.${IDD}.${CHR}.all.moderate.genes

done
