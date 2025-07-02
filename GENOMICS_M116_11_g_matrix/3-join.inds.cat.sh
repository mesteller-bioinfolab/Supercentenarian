#!/bin/bash

# Main Code --------------------------------------------------------

main() {
 
  rm -rf ./real/outputs/*

  # Configuration

  # All CHRs
    
   INPUT_DIR=./real/tmp2
   OUTPUT_DIR=./real/outputs
  
    #echo M116 > ./real/inputs/ids_list.txt
    #ls -l ../../1000GIBS/09-summary/real/outputs/ | grep txt | grep -v -w 116 | cut -d ' ' -f14 | cut -d '.' -f2 >> ./real/inputs/ids_list.txt 
    #IDS=$(cat ./real/inputs/ids_list.txt)
    #IDS_LIST=./real/inputs/ids_list.txt
  

  # Execution 
   
  # one job per CAT
 
 
#CATEGORIES=( damaging moderate modifier_cd modifier_nc altering low rare )

CAT_LIST=./real/inputs/cats_plus.txt

#CATEGORIES=$( cat ./real/inputs/cats.txt)

  # COMMAND
  COMMAND=" \
    ./real/scripts/join.inds.cat.sh \
      ${CAT_LIST} \
      ${INPUT_DIR} \
      ${OUTPUT_DIR}
  "



  # Cluster array execution
  JOBS_COUNT=$(cat ${CAT_LIST} | wc -l)
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
