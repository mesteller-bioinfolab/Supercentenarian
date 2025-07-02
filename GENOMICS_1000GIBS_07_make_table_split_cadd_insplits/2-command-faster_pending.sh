#!/bin/bash

# Main Code --------------------------------------------------------

main() {
  #rm -fR ./real/outputs/*
 
  ANN_KEY=./real/inputs/CSQ_KEY_

  TISSUE=1000GIBS

  # Configuration
  INPUT_DIR=./real/tmp
  OUTPUT_DIR=./real/tmp
  
  
  #CHR_ID_LIST=../02-filter-rare/pending
  #CHR_ID_LIST=herepending
  CHR_ID_LIST=pendingsplit  

  # COMMAND
  COMMAND=" \
    ./real/scripts/ann-split-insplits-lite_pending.sh \
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
