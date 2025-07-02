#!/bin/bash

# Main Code --------------------------------------------------------

main() {
  
 

  # get samples id
  TISSUE=1000GIBS

  # Configuration

  CHR_LIST=./real/inputs/chrs.list

  #INPUT_DIR=/gpfs42/projects/lab_anavarro/disease_pleiotropies/supercentenarian/1000GIBS/0-data/CCDG_14151_B01_GRM_WGS_2020-08-05_unnannotated
  
  INPUT_DIR=./real/outputs
  TMP_DIR=./real/tmp
  OUTPUT_DIR=./real/outputs
 
   #IDS=$(ls ${INPUT_DIR} | grep -v tbi)
  # Execution   
  INPUT_DIR_THEN=${INPUT_DIR}/${TISSUE}

  # COMMAND
  COMMAND=" \
    ./real/scripts/het.sh \
      ${CHR_LIST} \
      ${INPUT_DIR_THEN} \
      ${OUTPUT_DIR} \
      ${TMP_DIR} 
  "


  # Cluster array execution
  JOBS_COUNT=$(cat ${CHR_LIST} | wc -l)
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
