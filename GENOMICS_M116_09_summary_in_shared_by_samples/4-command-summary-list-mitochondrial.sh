#!/bin/bash



# Main Code --------------------------------------------------------

main() {  
  
 
  # Configuration
   ls ../09-summary-in-shared-by-samples/real/outputs/blood/M116/blood.M116.chr1.*.dput | cut -d '.' -f6 | cut -d '_' -f3,4 | grep -v all > ./real/inputs/cats.txt
   CATS=./real/inputs/cats.txt
   TISSUE="2TISSUES"
  
   INPUT_DIR=./real/outputs
   OUTPUT_DIR=./real/summary_tables
  
   #rm -rf ${OUTPUT_DIR}/* 
  
 
    #IDS=$(ls ${INPUT_DIR}/${TISSUE} )
  

  ID=M116


  # Execution 
  #mkdir ${OUTPUT_DIR}/${TISSUE}
  #rm -fR ./real/outputs/${TISSUE}


  # one job per ID in tissue
  



  # COMMAND
  COMMAND=" \
    ./real/scripts/summary-lists-mitochondrial.sh \
      ${INPUT_DIR} \
      ${OUTPUT_DIR} \
      ${ID} \
      ${TISSUE} \
      ${CATS}
  "

  # Cluster array execution
  #JOBS_COUNT=$(cat ${CHR_LIST} | wc -l)
  #eval sbatch --array=1-${JOBS_COUNT} ${COMMAND}
  sbatch ${COMMAND}
  

 
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
