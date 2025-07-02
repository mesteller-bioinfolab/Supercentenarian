#!/bin/bash

# Main Code --------------------------------------------------------

main() {
  mkdir ./real/tmp2
  rm -rf ./real/tmp2/*


  # Configuration

    INPUT_DIR=./real/tmp  

    OUTPUT_DIR=./real/tmp2  
 
    echo M116 > ./real/inputs/ids_list.txt
    ls -l ../../1000GIBS/09-summary/real/outputs/  | grep genes | cut -d '.' -f2 | sort -u   >> ./real/inputs/ids_list.txt

    IDS_LIST=./real/inputs/ids_list.txt
  
  # Execution

  cat ./real/inputs/cats.txt > ./real/inputs/cats_plus.txt
  echo rare >>  ./real/inputs/cats_plus.txt

  CAT_LIST=./real/inputs/cats_plus.txt

  # COMMAND
  COMMAND=" \
    ./real/scripts/join.chrs.cat.sh \
      ${CAT_LIST} \
      ${INPUT_DIR} \
      ${OUTPUT_DIR} \
      ${IDS_LIST}
  "


  # Cluster array execution
  JOBS_COUNT=$(cat ${IDS_LIST} | wc -l)
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
