#!/bin/bash

# Main Code --------------------------------------------------------

main() {
 
  ANN_KEY=./real/inputs/CSQ_KEY_
  TISSUE=1000GIBS
  INPUT_DIR=./real/tmp
  OUTPUT_DIR=./real/tmp
  #CHR_SPLIT_LIST=HG01697.chr2.pending.chrs_split.list
  #CHR_SPLIT_LIST=HG02237.chr4.pending.chrs_split.list
  #CHR_SPLIT_LIST=HG02225.chr21.pending.chrs_split.list 
  #CHR_SPLIT_LIST=HG01762.chr2.pending.chrs_split.list
  #CHR_SPLIT_LIST=HG01501.chr1.pending.chrs_split.list
  CHR_SPLIT_LIST=HG01618.chr21.pending.chrs_split.list 


  # COMMAND
  COMMAND=" \
    ./real/scripts/ann-split-insplits-failed.sh \
      ${CHR_SPLIT_LIST} \
      ${INPUT_DIR} \
      ${OUTPUT_DIR} \
      ${ANN_KEY}
  "

  # Cluster array execution
  JOBS_COUNT=$(cat ${CHR_SPLIT_LIST} | wc -l)
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
