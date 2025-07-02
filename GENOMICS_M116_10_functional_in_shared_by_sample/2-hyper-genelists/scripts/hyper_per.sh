#!/bin/bash
#
#SBATCH -p normal # partition (queue)
#SBATCH -N 1 # number of nodes
#SBATCH -J hyperper
#SBATCH --mem 16G # memory pool for all cores
#SBATCH -t 0-00:59# time (D-HH:MM)
#SBATCH -o ./real/outputs/log.%j.out # STDOUT
#SBATCH -e ./real/outputs/log.%j.err # STDERR
#SBATCH --mail-type=END,FAIL # notifications for job done & fail
#SBATCH --mail-user=claudia.vasallo@upf.edu # send-to address



#GENELISTS_FOLDER=$1



# For all chromosomes together

#TABLE=$1   #../../09-summary-in-shared-by-samples/real/outputs
INPUTS=$2  #../../09-summary-in-shared-by-samples/real/summary_tables
OUTPUTS=$3  #./real/outputs
GENELISTS_FOLDER=$4
CAT_BG_FOLDER=$5     #/gpfs42/projects/lab_anavarro/disease_pleiotropies/supercentenarian/Olot/10-functional-in-shared-by-sample/x-goof_red_genelists/real/outputs
CAT_FILE=$6
TYPE=$7


#CATEGORIES=( damaging  moderate modifier_cd modifier_nc altering low_only )

CATEGORIES=$( cat ${CAT_FILE} )

CODING_GENES=Homo_sapiens.GRCh38.104.protein_coding.txt


# Results per category will go here: script hyper_pool.R



# All longevity lists together (Manel gene lists) -------------------------------------------------------------------------------------------




# All the genelists to test
GENELISTS=$(ls ${GENELISTS_FOLDER}/*.txt | cut -d '/' -f6 | cut -d '.' -f1)


# POPULATION (All Genes in M116 sample in LEVEL: here shared in 2 or 3 samples)
#cat ${LEVEL} | grep Transcript | grep protein_coding | cut -f16 | sort -u | grep -v -w "" > ./real/inputs/ALLGENES_LIST
#ALLGENES=$(cat ./real/inputs/ALLGENES_LIST | wc -l)


#----


#for CATEGORY in ${CATEGORIES[@]}
#do
#BG_GENES=${CAT_BG_FOLDER}/everyones.${CATEGORY}.genes
#BG=$(cat ${BG_GENES} | wc -l)


 
for GENELIST in ${GENELISTS[@]}

do


echo CATEGORY BG GENELIST_IN_BG CAT GENELIST_IN_CAT PVAL OR   > ${OUTPUTS}/hyper.results_${GENELIST}.txt


for CATEGORY in ${CATEGORIES[@]}

do


all="all"

if [[ "$CATEGORY" == "$all" ]]
then
    	CATEGORY_GENES_ORI=${INPUTS}/eur_rare.genes
        echo "all"
else
    	CATEGORY_GENES_ORI=${INPUTS}/eur_rare_${CATEGORY}.genes
fi

echo $CATEGORY



CODE="CODING"
if [[ $TYPE == "$CODE" ]]
then
	BG_GENES_ORI=${CAT_BG_FOLDER}/everyones.${CATEGORY}.genes
	cat ${BG_GENES_ORI} | grep -f ${CODING_GENES} >	./real/inputs/everyones.${CATEGORY}.coding.genes.per
	BG_GENES=./real/inputs/everyones.${CATEGORY}.coding.genes.per
else
	BG_GENES=${CAT_BG_FOLDER}/everyones.${CATEGORY}.genes
fi



BG=$(cat ${BG_GENES} | wc -l)


#CATEGORY_GENES=${INPUTS}/eur_rare_${CATEGORY}.genes

cat ${CATEGORY_GENES_ORI} | grep -f ${BG_GENES} > ./real/inputs/CAT_${CATEGORY}.${GENELIST}.genes
CATEGORY_GENES=./real/inputs/CAT_${CATEGORY}.${GENELIST}.genes
CAT=$(cat ${CATEGORY_GENES} | sort -u | wc -l)


# SUCCESES IN SAMPLE (Longevity genes in Genes with rare and altering/damaging in M116 sample)
GENELIST_GENES=${GENELISTS_FOLDER}/${GENELIST}.txt
GENELIST_IN_CAT=$( cat ${CATEGORY_GENES} | grep -f ${GENELIST_GENES} | wc -l)


# SUCCESSES IN POPULATION -GENES iN POPULATION belonging IN GENELIST
cat ${GENELIST_GENES} | grep -f ${BG_GENES} | sort -u >  ./real/inputs/${GENELIST}_IN_BG.genes
GENELIST_IN_BG_GENES=./real/inputs/${GENELIST}_IN_BG.genes
GENELIST_IN_BG=$(cat ${GENELIST_IN_BG_GENES}  | grep -v -i -w ensg | cut -f2 | sort -u | wc -l)


module load R

HYPER=$(Rscript ./real/scripts/hyper.R ${BG} ${GENELIST_IN_BG} ${CAT} ${GENELIST_IN_CAT})


PVAL=$(echo $HYPER | cut -d ' ' -f2)
OR=$(echo $HYPER | cut -d ' ' -f3)


echo $CATEGORY ${BG} ${GENELIST_IN_BG} ${CAT} ${GENELIST_IN_CAT} ${PVAL} ${OR} >> ${OUTPUTS}/hyper.results_${GENELIST}.txt


done  # done CATEGORIES


done  # done GENELISTS




#phyper(q, m, n, k, lower.tail = TRUE, log.p = FALSE)
#x, q vector of quantiles representing the number of white balls drawn
#without replacement from an urn which contains both black and white
#balls.
#m the number of white balls in the urn.
#n the number of black balls in the urn.
#k the number of balls drawn from the urn.
