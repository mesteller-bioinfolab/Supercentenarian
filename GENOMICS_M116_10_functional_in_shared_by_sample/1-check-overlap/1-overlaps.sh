#!/bin/bash

# Main Code --------------------------------------------------------

main() {
  rm -rf real/tmp/*
  #rm -fR ./real/outputs/*
  
  # Genelists
 
  cp $( ls ../0-gene-lists/real/outputs/manel/*) ./real/tmp  
  cp $( ls ../0-gene-lists/real/outputs/papers/* | grep -v chinese ) ./real/tmp


  cat $( ls ./real/tmp/*_ensg.txt | grep -v paper | grep -v all )| sort -u > ./real/tmp/all_longevity_ensg.txt  


  GENELISTS_DIR=./real/tmp
  TISSUE="2TISSUES"



  CHR_LIST=../../02-vep-annot/real/inputs/chrs.list

  INPUT_DIR=/gpfs42/projects/lab_anavarro/disease_pleiotropies/supercentenarian/Olot/09-summary-in-shared-by-samples/real/summary_tables

 
  OUTPUT_DIR=./real/outputs
    
  ID=M116
 
  CAT_FILE=./real/inputs/cats.txt

 
  # Execution 
  mkdir ${OUTPUT_DIR}


  # COMMAND
  COMMAND=" \
    ./real/scripts/overlap.sh \
      ${INPUT_DIR} \
      ${OUTPUT_DIR} \
      ${TISSUE} \
      ${ID} \
      ${GENELISTS_DIR} \
      ${CAT_FILE} \
      ${AUTOSOMAL_GENES_FILE}
  "

  # Cluster array execution
  #JOBS_COUNT=$(cat ${ID_LIST} | wc -l)
  #eval sbatch --array=1-${JOBS_COUNT} ${COMMAND}
  #exit


  # Cluster execution
  eval sbatch ${COMMAND}
  exit             

} 


# FUNCTIONS ==========================================================

main
