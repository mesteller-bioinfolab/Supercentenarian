rm -rf real/outputs/*


TYPE="CATBG"


# Compund bg of same category genes intersect with coding genes list downoaded from Ensembl 18K
# forces the categories into coding
#TYPE="CODING" 

TABLE_DIR=../../09-summary-in-shared-by-samples/real/outputs
INPUT_DIR=../../09-summary-in-shared-by-samples/real/summary_tables
OUTPUT_DIR=./real/outputs
GENELISTS_DIR=../0-gene-lists/real/outputs/manel
CAT_BG_DIR=/gpfs42/projects/lab_anavarro/disease_pleiotropies/supercentenarian/Olot/10-functional-in-shared-by-sample/x-goof_red_genelists/real/outputs
CATS=./real/inputs/versions.txt


#cat ../../09-summary-in-shared-by-samples/real/outputs/*/M116/*.M116.*.eur_rare_all.fulltable | grep -v -w chrX | grep -v -w chrY | grep -v -w chrM | grep Transcript  | grep CODING-GENE | grep -v -w ID | cut -f16 | sort -u | grep -v -w "" > ./real/inputs/BG_M116.coding.genes
#cat ../../../1000GIBS/08-filter-category-summary/real/outputs/1000GBIB/*/*.eur_rare_all.fulltable | grep -v -w chrX | grep -v -w chrY | grep -v -w chrM | grep Transcript  | grep CODING-GENE | grep -v -w ID | cut -f16 | sort -u | grep -v -w "" > ./real/inputs/BG_1000GIBS.coding.genes
#

#cat ./real/inputs/BG_M116.coding.genes > ./real/inputs/BG_everyones.coding.genes
#cat ./real/inputs/BG_1000GIBS.coding.genes >> ./real/inputs/BG_everyones.coding.genes


# enrichment overall longevity genelists
sbatch ./real/scripts/hyper_pool.sh ${TABLE_DIR} ${INPUT_DIR} ${OUTPUT_DIR} ${GENELISTS_DIR} ${CAT_BG_DIR} ${CATS} ${TYPE}

# enrichment per genelist
sbatch ./real/scripts/hyper_per.sh ${TABLE_DIR} ${INPUT_DIR} ${OUTPUT_DIR} ${GENELISTS_DIR} ${CAT_BG_DIR} ${CATS} ${TYPE}

