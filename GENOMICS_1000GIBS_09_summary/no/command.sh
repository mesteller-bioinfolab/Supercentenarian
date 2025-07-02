#!/bin/bash

# Main Code --------------------------------------------------------

main() {
  rm -fR ./real/outputs/*


  TISSUE=1000GIBS

  # Configuration
  
   
  INPUT_DIR=/gpfs42/projects/lab_anavarro/disease_pleiotropies/supercentenarian/1000GIBS/08-filter-category-summary/real/outputs
  OUTPUT_DIR=./real/outputs
  
  
  
   #later>IDS=$(ls ${INPUT_DIR}/${TISSUE} | grep -v -f ./real/inputs/ids_already.txt)
   #IDS=$(cat ./real/inputs/ids_already2.txt)
   #IDS=$(ls ${INPUT_DIR}/${TISSUE} )



    ls ${INPUT_DIR}/${TISSUE} > ./real/inputs/id.list
    ID_LIST=$(cat ./real/inputs/id.list)

  # Execution 
  #mkdir ${OUTPUT_DIR}/${TISSUE}


  # COMMAND
  COMMAND=" \
    ./real/scripts/summary-lists.sh \
      ${INPUT_DIR} \
      ${OUTPUT_DIR} \
      ${TISSUE} \
      ${ID_LIST}
  "

  # Cluster array execution
  JOBS_COUNT=$(cat ${ID_LIST} | wc -l)
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
