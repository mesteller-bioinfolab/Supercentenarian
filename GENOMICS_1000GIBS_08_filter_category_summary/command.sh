#!/bin/bash

# Main Code --------------------------------------------------------

main() {
  rm -fR ./real/outputs/*


  TISSUE=1000GIBS

  # Configuration
  
  CHR_LIST=./real/inputs/chrs.list1
  
  INPUT_DIR=/gpfs42/projects/lab_anavarro/disease_pleiotropies/supercentenarian/1000GIBS/07-make-table-split_cadd-insplits/real/outputs
  OUTPUT_DIR=./real/outputs
  
  SUFF=table.split2
  
   # this filter is redundant here but for the sake of homogeneity and use of the same code 
   FILTER_DIR=../00-data/common_snvs/real/outputs


   #later>IDS=$(ls ${INPUT_DIR}/${TISSUE} | grep -v -f ./real/inputs/ids_already.txt)
   #IDS=$(cat ./real/inputs/ids_already2.txt)
 
  # do not use HG01507 which is empty of rare variants
  # IDS=$(ls ${INPUT_DIR}/${TISSUE} | grep -v HG01507)

   ID_LIST=./real/inputs/ids.list

  # Execution 
  mkdir ${OUTPUT_DIR}/${TISSUE}
  #rm -fR ./real/outputs/${TISSUE}




  # COMMAND
  COMMAND=" \
    ./real/scripts/run_summary.sh \
      ${CHR_LIST} \
      ${INPUT_DIR} \
      ${OUTPUT_DIR} \
      ${TISSUE} \
      ${ID_LIST} \
      ${SUFF} \
      ${FILTER_DIR}
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
