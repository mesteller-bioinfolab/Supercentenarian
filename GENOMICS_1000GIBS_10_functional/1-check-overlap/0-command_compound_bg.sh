#!/bin/bash

# Main Code --------------------------------------------------------

main() {

mkdir real/outputs/compound_bg  

rm -rf real/outputs/compound_bg/*


  GENELISTS_DIR=./real/tmp
  TISSUE=1000GIBS
  
 

  INPUT_DIR=/gpfs42/projects/lab_anavarro/disease_pleiotropies/supercentenarian/1000GIBS/09-summary/real/outputs

  #ls -l ../../09-summary/real/outputs/ | grep txt | cut -d '.' -f2 | sort -u > ./real/inputs/id.list
  ID_LIST=./real/inputs/id.list

  CAT_FILE=./real/inputs/cats.txt
  CAT_BG_DIR=/gpfs42/projects/lab_anavarro/disease_pleiotropies/supercentenarian/Olot/10-functional-in-shared-by-sample/x-goof_red_genelists/real/outputs


OUTPUT_DIR=./real/outputs/compound_bg
mkdir ${OUTPUT_DIR}


  # COMMAND
  COMMAND="./real/scripts/compound_bg.sh \
      ${TISSUE} \
      ${ID_LIST} \
      ${INPUT_DIR} \
      ${OUTPUT_DIR} \
      ${CAT_BG_DIR}
    "

  JOBS_COUNT=$(cat ${ID_LIST} | wc -l)
  eval sbatch --array=1-${JOBS_COUNT} ${COMMAND}
 
  #sbatch ${COMMAND}

exit




} 


# FUNCTIONS ==========================================================

main


