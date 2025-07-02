#!/bin/bash
#
#SBATCH --partition=normal   #haswell       # normal  # haswell-p normal # partition (queue)
#SBATCH -N 1 # number of nodes
#SBATCH -J autotab-summary
#SBATCH --mem 64G # memory pool for all cores
#SBATCH -t 0-00:29# time (D-HH:MM)
#SBATCH -o ./real/summary_tables/log.%j.out # STDOUT
#SBATCH -e ./real/summary_tables/log.%j.err # STDERR
#SBATCH --mail-type=END,FAIL # notifications for job done & fail
#SBATCH --mail-user=claudia.vasallo@upf.edu # send-to address



INPUTS=./real/summary_tables
TISSUED=2TISSUES
IDD=M116



#AUTO_GENES=$1
#SEX_GENES=$2
#MT_GENES=$3
#AUTO_VARIANTS=$4
#SEX_VARIANTS=$5
#MT_VARIANTS=$6


ARRAY_LIST=$1


if [[ "${SLURM_ARRAY_TASK_ID}" == "" ]]
then
	CAT=$(cat ${ARRAY_LIST} | sed -n 1p)
else
	CAT=$(cat ${ARRAY_LIST} | sed -n ${SLURM_ARRAY_TASK_ID}p)
fi


echo ${CAT}


# autosomal genes of individual
#cat real/outputs/*/${IDD}/*.${IDD}.chr*autosomal.genes | sort -u  | grep  ENSG > ${INPUTS}/${TISSUED}.${IDD}.autosomal.genes
#AUTO_GENES=./real/summary_tables/${TISSUED}.${IDD}.autosomal.genes
#cat ${INPUTS}/${TISSUED}.${IDD}.eur_rare.variants | grep -v chrX | grep -v Y | grep -v chrM > ${INPUTS}/${TISSUED}.${IDD}.eur_rare.autosomal.variants
#cat ${INPUTS}/${TISSUED}.${IDD}.eur_rare.variants | grep -v -f ${INPUTS}/${TISSUED}.${IDD}.eur_rare.autosomal.variants  > ${INPUTS}/${TISSUED}.${IDD}.eur_rare.sex-mt.variants
#CATEGORIES=( $(cat ./real/inputs/cats.txt) )
#for CAT in ${CATEGORIES[@]}
#do
# VARIANTS
#cat ${INPUTS}/${TISSUED}.${IDD}.eur_rare_${CAT}.variants | grep -v chrX | grep -v Y  | grep -v chrM >  ${INPUTS}/${TISSUED}.${IDD}.eur_rare_${CAT}.autosomal.variants
#cat ${INPUTS}/${TISSUED}.${IDD}.eur_rare_${CAT}.variants | grep -v -f ${INPUTS}/${TISSUED}.${IDD}.eur_rare.autosomal.variants  > ${INPUTS}/${TISSUED}.${IDD}.eur_rare_${CAT}.sex-mt.variants



#2TISSUES.M116.eur_rare_ncnc_cadd15.variants
# VARIANTS
#cat ${INPUTS}/${TISSUED}.${IDD}.eur_rare_${CAT}.variants | grep -f ${AUTO_VARIANTS} > ${INPUTS}/${TISSUED}.${IDD}.eur_rare_${CAT}.autosomal.variants
#cat ${INPUTS}/${TISSUED}.${IDD}.eur_rare_${CAT}.variants | grep -f ${SEX_VARIANTS} > ${INPUTS}/${TISSUED}.${IDD}.eur_rare_${CAT}.sex.variants
#cat ${INPUTS}/${TISSUED}.${IDD}.eur_rare_${CAT}.variants | grep -f ${MT_VARIANTS} > ${INPUTS}/${TISSUED}.${IDD}.eur_rare_${CAT}.mt.genes

cat $(ls  real/outputs/*/${IDD}/*.eur_rare_${CAT}.variants | grep -v chrX | grep -v chrY | grep -v chrM) | sort -u > ${INPUTS}/${TISSUED}.${IDD}.eur_rare_${CAT}.autosomal.variants
cat $(ls  real/outputs/*/${IDD}/*.eur_rare_${CAT}.variants | grep chrX )  | sort -u > ${INPUTS}/${TISSUED}.${IDD}.eur_rare_${CAT}.sex.variants
cat $(ls  real/outputs/*/${IDD}/*.eur_rare_${CAT}.variants | grep chrM ) | sort -u > ${INPUTS}/${TISSUED}.${IDD}.eur_rare_${CAT}.mt.variants

# GENES

#cat ${INPUTS}/${TISSUED}.${IDD}.eur_rare_${CAT}.genes | grep -f ${AUTO_GENES} > ${INPUTS}/${TISSUED}.${IDD}.eur_rare_${CAT}.autosomal.genes
#cat ${INPUTS}/${TISSUED}.${IDD}.eur_rare_${CAT}.genes | grep -v -f ${SEX_GENES} > ${INPUTS}/${TISSUED}.${IDD}.eur_rare_${CAT}.sex.genes
#cat ${INPUTS}/${TISSUED}.${IDD}.eur_rare_${CAT}.genes | grep -v -f ${MT_GENES} > ${INPUTS}/${TISSUED}.${IDD}.eur_rare_${CAT}.mt.genes

#cat $(ls  real/outputs/*/${IDD}/*.eur_rare_${CAT}.genes | grep -v chrX | grep -v chrY | grep -v chrM ) > ${INPUTS}/${TISSUED}.${IDD}.eur_rare_${CAT}.autosomal.genes
#cat $(ls  real/outputs/*/${IDD}/*.eur_rare_${CAT}.genes | grep chrX ) > ${INPUTS}/${TISSUED}.${IDD}.eur_rare_${CAT}.sex.genes 
#cat $(ls  real/outputs/*/${IDD}/*.eur_rare_${CAT}.genes	| grep chrM ) > ${INPUTS}/${TISSUED}.${IDD}.eur_rare_${CAT}.mt.genes

cat $(ls real/outputs/*/${IDD}/*.${IDD}.*.genevar | grep ${CAT} | grep -v chrM | grep -v chrY | grep -v chrX) | cut -f1  | sort -u | grep  ENSG | grep -v -w ENSG > ${INPUTS}/${TISSUED}.${IDD}.eur_rare_${CAT}.autosomal.genes
cat $(ls real/outputs/*/${IDD}/*.${IDD}.*.genevar | grep ${CAT} | grep chrX) | cut -f1  | sort -u | grep  ENSG | grep -v -w ENSG > ${INPUTS}/${TISSUED}.${IDD}.eur_rare_${CAT}.sex.genes
cat $(ls real/outputs/*/${IDD}/*.${IDD}.*.genevar | grep ${CAT} | grep chrM) | cut -f1  | sort -u | grep  ENSG | grep -v -w ENSG > ${INPUTS}/${TISSUED}.${IDD}.eur_rare_${CAT}.mt.genes 


#echo all $(cat ${INPUTS}/${TISSUED}.${IDD}.eur_rare.autosomal.variants | wc -l ) $(cat ${INPUTS}/${TISSUED}.${IDD}.eur_rare.autosomal.genes | wc -l ) >> ${INPUTS}/${TISSUED}.${IDD}.autosomal.numbers.txt


echo ${CAT} $(cat ${INPUTS}/${TISSUED}.${IDD}.eur_rare_${CAT}.autosomal.variants | wc -l) $(cat ${INPUTS}/${TISSUED}.${IDD}.eur_rare_${CAT}.autosomal.genes | wc -l ) >> ${INPUTS}/${TISSUED}.${IDD}.autosomal.numbers.txt


#echo all $(cat ${INPUTS}/${TISSUED}.${IDD}.eur_rare.sex.variants | wc -l ) $(cat ${INPUTS}/${TISSUED}.${IDD}.eur_rare.sex.genes | wc -l ) >> ${INPUTS}/${TISSUED}.${IDD}.sex.numbers.txt
#echo all $(cat ${INPUTS}/${TISSUED}.${IDD}.eur_rare.mt.variants | wc -l ) $(cat ${INPUTS}/${TISSUED}.${IDD}.eur_rare.mt.genes | wc -l ) >> ${INPUTS}/${TISSUED}.${IDD}.mt.numbers.txt


echo ${CAT} $(cat ${INPUTS}/${TISSUED}.${IDD}.eur_rare_${CAT}.sex.variants | wc -l) $(cat ${INPUTS}/${TISSUED}.${IDD}.eur_rare_${CAT}.sex.genes | wc -l ) >> ${INPUTS}/${TISSUED}.${IDD}.sex.numbers.txt


echo ${CAT} $(cat ${INPUTS}/${TISSUED}.${IDD}.eur_rare_${CAT}.mt.variants | wc -l) $(cat ${INPUTS}/${TISSUED}.${IDD}.eur_rare_${CAT}.mt.genes | wc -l ) >> ${INPUTS}/${TISSUED}.${IDD}.mt.numbers.txt



