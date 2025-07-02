#!/bin/bash

# Main Code --------------------------------------------------------

main() {
  #rm -fR ./real/outputs/*


  TISSUE=1000GIBS

  # Configuration

  CHR_ID_LIST=../02-filter-rare/pending 

  INPUT_DIR=/gpfs42/projects/lab_anavarro/disease_pleiotropies/supercentenarian/1000GIBS/04-cadd-annot/real/outputs
  OUTPUT_DIR=./real/outputs
  

  # one job per ID in tissue
  

  # COMMAND
  COMMAND=" \
    ./real/scripts/vep2table_pending.sh \
      ${CHR_ID_LIST} \
      ${INPUT_DIR} \
      ${OUTPUT_DIR} \
      ${TISSUE}
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
