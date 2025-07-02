rm -rf real/outputs/*

TYPE="CATBG"


TABLE_DIR=../../09-summary-in-shared-by-samples/real/outputs
INPUT_DIR=../../09-summary-in-shared-by-samples/real/summary_tables
OUTPUT_DIR=./real/outputs
GENELISTS_DIR=../0-gene-lists/real/outputs/papers
CAT_BG_DIR=/gpfs42/projects/lab_anavarro/disease_pleiotropies/supercentenarian/Olot/10-functional-in-shared-by-sample/x-goof_red_genelists/real/outputs
CATS=real/inputs/versions.txt



sbatch ../2-hyper-genelists/real/scripts/hyper_pool.sh ${TABLE_DIR} ${INPUT_DIR} ${OUTPUT_DIR} ${GENELISTS_DIR} ${CAT_BG_DIR} ${CATS} ${TYPE}


sbatch ../2-hyper-genelists/real/scripts/hyper_per.sh ${TABLE_DIR} ${INPUT_DIR} ${OUTPUT_DIR} ${GENELISTS_DIR} ${CAT_BG_DIR} ${CATS} ${TYPE}

