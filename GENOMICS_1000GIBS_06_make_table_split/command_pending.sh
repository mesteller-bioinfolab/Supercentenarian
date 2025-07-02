#!/bin/bash

# Main Code --------------------------------------------------------

main() {
  #rm -fR ./real/outputs/*
 
  ANN_KEY=./real/inputs/CSQ_KEY_
  TISSUE=1000GIBS

  # Configuration 
  CHR_ID_LIST=../02-filter-rare/pending

  INPUT_DIR=/gpfs42/projects/lab_anavarro/disease_pleiotropies/supercentenarian/1000GIBS/05-make-table/real/outputs
  OUTPUT_DIR=./real/outputs
  
 

  # COMMAND
  COMMAND=" \
    ./real/scripts/ann-split_pending.sh \
      ${CHR_ID_LIST} \
      ${INPUT_DIR} \
      ${OUTPUT_DIR} \
      ${TISSUE} \
      ${ANN_KEY}
  "

  # Cluster array execution
  JOBS_COUNT=$(cat ${CHR_ID_LIST} | wc -l)
  eval sbatch --array=1-${JOBS_COUNT} ${COMMAND}

 
 exit


  # Direct execution
  #eval bash ${COMMAND}
  #exit 

  # Cluster execution
  #eval sbatch ${COMMAND}
  #exit             

} 


# FUNCTIONS ==========================================================

main
