#!/bin/bash

# Main Code --------------------------------------------------------

main() {

  #rm -rf ./real/outputs/*

  TISSUE=urine

  # Configuration
  
 

  CHR_LIST=../02-vep-annot/real/inputs/chrs.list
 

  INPUT_DIR=/gpfs42/projects/lab_anavarro/disease_pleiotropies/supercentenarian/Olot/02-vep-annot/real/outputs
  
  OUTPUT_DIR=./real/outputs
  
  
 
   IDS=$(ls ${INPUT_DIR}/${TISSUE})

 
  # Execution 
  mkdir ${OUTPUT_DIR}/${TISSUE}


  # one job per ID in tissue
  
for ID in ${IDS}

do

echo $ID

mkdir ${OUTPUT_DIR}/${TISSUE}/${ID}


  # COMMAND
  COMMAND=" \
    ./real/scripts/list_filter.sh \
      ${INPUT_DIR} \
      ${OUTPUT_DIR} \
      ${CHR_LIST} \
      ${TISSUE} \
      ${ID}
  "

  # Cluster array execution
  JOBS_COUNT=$(cat ${CHR_LIST} | wc -l)
  eval sbatch --array=1-${JOBS_COUNT} ${COMMAND}

done
 
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
